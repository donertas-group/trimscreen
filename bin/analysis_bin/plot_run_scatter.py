#!/usr/bin/env python3
# example usage:

# ./plot_run_scatter.py -D ulas -C Genus_pasv

import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Plot run-level statistics with coordinate extraction from run IDs.")
    parser.add_argument("-D", "--dataset", required=True, type=str, help="Dataset directory name (e.g. mock16)")
    parser.add_argument("-C", "--column", required=True, type=str, choices=["retained_reads_percent", "Phylum_pasv", "Family_pasv", "Genus_pasv"],
                        help="Summary metric to visualize (e.g., Genus_pasv)")
    return parser.parse_args()


def read_files(dataset):
    result_dir = "/scratch/shire/ssd/pipeline/16s_nf_pipeline"  # Modify if needed
    base_path = os.path.join(result_dir, dataset, "output/compare_runs")

    full_file = os.path.join(base_path, "full_table.csv")
    filtered_file = os.path.join(base_path, "filtered_table.csv")

    full_table = pd.read_csv(full_file)
    filtered_table = pd.read_csv(filtered_file)

    return full_table, filtered_table, result_dir


def extract_xy(run_id):
    """Extract x, y coordinates from run_###### pattern."""
    match = re.match(r"run_(\d{6})", run_id)
    if match:
        digits = match.group(1)
        return int(digits[:3]), int(digits[3:])
    else:
        return None, None


def plot_results(full_table, filtered_table, column, dataset, result_dir):
    # Compute average per run
    full_avg = full_table.groupby("run")[column].mean().reset_index()
    filt_avg = filtered_table.groupby("run")[column].mean().reset_index()

    # Determine which runs are filtered
    filt_runs = set(filt_avg["run"])
    full_avg["is_filtered"] = full_avg["run"].isin(filt_runs)

    # Extract coordinates
    full_avg[["x", "y"]] = full_avg["run"].apply(lambda r: pd.Series(extract_xy(r)))

    # Split data
    unfiltered = full_avg[~full_avg["is_filtered"]]
    filtered = full_avg[full_avg["is_filtered"]]

    # Compute shared color normalization
    vmin = full_avg[column].min()
    vmax = full_avg[column].max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    # Make the scatter plot
    plt.figure(figsize=(8, 6))

    # Plot unfiltered runs (semi-transparent)
    plt.scatter(
        unfiltered["x"],
        unfiltered["y"],
        c=unfiltered[column],
        cmap="viridis",
        norm=norm,
        alpha=0.5,
        edgecolor="none"
    )

    # Plot filtered runs (solid color, black border)
    scatter = plt.scatter(
        filtered["x"],
        filtered["y"],
        c=filtered[column],
        cmap="viridis",
        norm=norm,
        alpha=1.0,
        edgecolor="black",
        linewidth=0.5
    )

    # Shared colorbar
    plt.colorbar(scatter, label=column)
    plt.title(f"Run-wise {column} averages ({dataset})")
    plt.xlabel("forward truncation length")
    plt.ylabel("reverse truncation length")

    output_dir = os.path.join(result_dir, "analysis_mock/output")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{dataset}.{column}.png")

    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Plot saved to: {output_file}")

def fplot_results(full_table, filtered_table, column, dataset, result_dir):
    # Compute average per run
    full_avg = full_table.groupby("run")[column].mean().reset_index()
    filt_avg = filtered_table.groupby("run")[column].mean().reset_index()

    # Determine which runs are filtered
    filt_runs = set(filt_avg["run"])
    full_avg["is_filtered"] = full_avg["run"].isin(filt_runs)

    # Extract coordinates
    full_avg[["x", "y"]] = full_avg["run"].apply(lambda r: pd.Series(extract_xy(r)))

    # Make the scatter plot
    plt.figure(figsize=(8, 6))

    # Split data
    unfiltered = full_avg[~full_avg["is_filtered"]]
    filtered = full_avg[full_avg["is_filtered"]]

    # Plot unfiltered runs (semi-transparent)
    plt.scatter(
        unfiltered["x"],
        unfiltered["y"],
        c=unfiltered[column],
        cmap="viridis",
        alpha=0.5,
        edgecolor="none"
    )

    # Plot filtered runs (solid color, black border)
    scatter = plt.scatter(
        filtered["x"],
        filtered["y"],
        c=filtered[column],
        cmap="viridis",
        alpha=1.0,
        edgecolor="black",
        linewidth=0.5
    )

    plt.colorbar(scatter, label=column)
    plt.title(f"Run-wise {column} averages ({dataset})")
    plt.xlabel("forward truncation length")
    plt.ylabel("reverse truncation length")

    output_dir = os.path.join(result_dir, "analysis_mock/output")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{dataset}.{column}.png")

    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Plot saved to: {output_file}")


def main():
    args = parse_args()
    full_table, filtered_table, result_dir = read_files(args.dataset)
    plot_results(full_table, filtered_table, args.column, args.dataset, result_dir)


if __name__ == "__main__":
    main()

