#!/usr/bin/env python
import pandas as pd
import sys
import math
import argparse
import json 
from itertools import combinations

#######################
# This script subsets the filtered table by keeping only one sample, either the first one of the replicates or the first sample when replicates are not present.
# Then it rank the runs by metrics provided, e.g. retained_reads_percent, shannon_genus, and rep_similarity.
# It prints out the optimal run based on the mean of all ranks, and prints out a report showing the rank scores under each metric.
#######################

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Find the best run")
    parser.add_argument("-i", "--input", required=True, help="filtered_table_csv")
    parser.add_argument("-m", "--metadata", required=False, help="Optional metadata tsv table with same format as required by nf-core/ampliseq")
    parser.add_argument("--metrics", required=True, nargs="+", help="Trimming-dependent metrics on which runs are evaluated and ranked")

    return parser.parse_args()

def main():
    args = parse_args()

    filtered_table_csv = args.input
    metadata_csv = args.metadata
    metrics = args.metrics

    df0 = pd.read_csv(filtered_table_csv)

    # get sampleIDs that are replicates
    rep_samples = None

    if metadata_csv is not None:
        metadata = pd.read_csv(metadata_csv)
    
        required_cols = {'condition', 'sampleID'}
        if not required_cols.issubset(metadata.columns):
            raise ValueError(f"Metadata must contain columns: {required_cols}")

        sample_ids = metadata[metadata['condition'] == 'sample']['sampleID'].unique()

        if len(sample_ids) == 0:
            raise ValueError("No samples found with condition == 'sample'")

        if 'is_replicate' in metadata.columns:
            rep_samples = metadata.loc[
                metadata['is_replicate'].astype(str).str.lower() == 'true', 'sampleID'
                ].unique()
            if len(rep_samples) > 0:
                df = df0[df0['sample'] == rep_samples[0]]
            else:
                df = df0[df0['sample'] == sample_ids[0]]
        else:
            df = df0[df0['sample'] == sample_ids[0]]

    else:
        sample_ids = df0['sample'].unique()
        if len(sample_ids) == 0:
            raise ValueError("No samples found in input table")

        df = df0[df0['sample'] == sample_ids[0]]


    # =========================
    # Rank runs by metrics
    # =========================

    # Ensure 'run' column exists
    if 'run' not in df.columns:
        raise ValueError("Input table must contain a 'run' column")

    ranks = pd.DataFrame({'run': df['run']})

    used_metrics = []

    for metric in metrics:
        if metric not in df.columns:
            print(f"Warning: metric '{metric}' not found in table, skipping", file=sys.stderr)
            continue

        # Skip metric if all values are None / NaN
        if df[metric].isna().all():
            print(f"Warning: metric '{metric}' contains only NaN, skipping", file=sys.stderr)
            continue

        # Higher value = better rank (rank 1 is best)
        ranks[metric + "_rank"] = df[metric].rank(
            ascending=False,
            method="average"
        )
        used_metrics.append(metric)

    if not used_metrics:
        raise ValueError("No valid metrics available for ranking")

    # =========================
    # Mean rank and optimal run
    # =========================

    rank_cols = [m + "_rank" for m in used_metrics]
    ranks["mean_rank"] = ranks[rank_cols].mean(axis=1)

    best_run_row = ranks.sort_values("mean_rank").iloc[0]
    best_run = best_run_row["run"]

    print(json.dumps([best_run]))

    # =========================
    # Write ranking report
    # =========================

    with open("report.txt", "w") as f:
        f.write("Run ranking report\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"Optimal run (lowest mean rank): {best_run}\n\n")

        for metric in used_metrics:
            f.write(f"Top 10 runs by {metric}:\n")
            f.write("-" * 40 + "\n")

            metric_df = (
                df[['run', metric]]
                .assign(rank=ranks[metric + "_rank"])
                .sort_values("rank")
                .head(10)
            )

            for _, row in metric_df.iterrows():
                f.write(
                    f"Run: {row['run']}, "
                    f"{metric}: {row[metric]}, "
                    f"rank: {row['rank']}\n"
                )

            f.write("\n")


if __name__ == "__main__":
    sys.exit(main())
