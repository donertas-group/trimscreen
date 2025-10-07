#!/usr/bin/env python3

# The script expects the input files in a structured directory, as described.
# It uses read count as the x-axis and calculates detection stats and f1 summary score at the user-specified taxonomic rank.
# Make sure column names match exactly, especially for ranks (Genus, Family, etc.) — they are case-sensitive.
# example usage:
# ./compare_all_runs.py -M 13 -r run_180190 -R Genus --true true_composition.csv

import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt
import glob
import json

def parse_args():
    parser = argparse.ArgumentParser(description="Extract preprocessing statistics of a run, averaging all samples")
    parser.add_argument("-D", "--dataset", required=True, type=int, help="Dataset directory name")
    parser.add_argument("-C", "--column", required=True, type=str, help="Summary metric to compare (e.g., Genus_pasv)")
    return parser.parse_args()

def read_files(dataset):
    result_dir = "/scratch/shire/ssd/pipeline/16s_nf_pipeline"  # Modify if needed

    full_file = os.path.join(result_dir, dataset, "output/compare_runs/full_table.csv")
    full_table = pd.read_csv(full_file)

    filtered_file = os.path.join(result_dir, dataset, "output/compare_runs/filtered_table.csv")
    filtered_table = pd.read_csv(filtered_file)

    return full_table, filtered_table

def plot_results(summary_table, rank):
    output_file = os.path.join(result_dir, "analysis_mock/output", f"{dataset}.{column}.png")

    plt.figure(figsize=(10, 6))
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Plot saved to {output_file}")


def main():
    args = parse_args()
    full_table, filtered_table = read_files(args.dataset)
    
    plot_results(full_table, filtered_table, args.column)


if __name__ == "__main__":
    main()

