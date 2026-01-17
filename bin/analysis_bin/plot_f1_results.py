#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Plot F1 results from per-sample CSV")
    parser.add_argument("--csv", required=True, help="CSV file produced by driver_compare_w_mock_allruns.py")
    parser.add_argument("--out", default="./", help="Output directory for plots")
    parser.add_argument("--rank", default="Genus", help="Taxonomic rank")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    df = pd.read_csv(args.csv)

    # Compute run-level mean F1 if desired
    df_mean = df.groupby("run_id")["f1"].mean().reset_index()
    print(df_mean)

    # --- Plotting ---
    plt.figure(figsize=(8, 6))

    # Normalize color by F1
    norm = plt.Normalize(vmin=df["f1"].min(), vmax=df["f1"].max())

    # Unfiltered points
    df_unfiltered = df[~df["filtered"]]
    plt.scatter(df_unfiltered["TP"], df_unfiltered["FP"],  # example axes
                c=df_unfiltered["f1"], cmap="viridis", norm=norm,
                alpha=0.5, edgecolor="none", label="Dropped")

    # Filtered points
    df_filtered = df[df["filtered"]]
    plt.scatter(df_filtered["TP"], df_filtered["FP"],
                c=df_filtered["f1"], cmap="viridis", norm=norm,
                alpha=1.0, edgecolor="black", linewidth=0.5, label="Kept")

    plt.colorbar(label=f"F1 score at {args.rank}")
    plt.xlabel("TP (example)")
    plt.ylabel("FP (example)")
    plt.title(f"Mock dataset F1 scores at {args.rank}")
    plt.legend()

    outfile = os.path.join(args.out, f"f1_scatter_plot_{args.rank}.png")
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Plot saved to {outfile}")

if __name__ == "__main__":
    main()

