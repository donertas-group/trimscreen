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
    parser.add_argument("--metric_directions", nargs="+", help="Metric directions: + for higher-is-better, - for lower-is-better")
    parser.add_argument("--metric_weights", nargs="+", type=float, help="Optional metric weights (default = 1.0)")

    return parser.parse_args()

def main():
    args = parse_args()

    filtered_table_csv = args.input
    metadata_csv = args.metadata
    metrics = args.metrics
    directions = args.metric_directions
    weights = args.metric_weights

    df0 = pd.read_csv(filtered_table_csv)
    sample_ids = df0['sample'].unique()

    # get sampleIDs that are replicates
    rep_samples = None

    if metadata_csv is not None:
        metadata = pd.read_csv(metadata_csv)
    
        required_cols = {'condition', 'sampleID'}
        if not required_cols.issubset(metadata.columns):
            raise ValueError(f"Metadata must contain columns: {required_cols}")

        if 'is_replicate' in metadata.columns:
            rep_samples = metadata.loc[
                metadata['is_replicate'].astype(str).str.lower() == 'true', 'sampleID'
                ].unique()

            replicates = list(set(rep_samples).intersection(sample_ids))

            if len(replicates) > 0:
                df = df0[df0['sample'] == replicates[0]]
            else:
                df = df0[df0['sample'] == sample_ids[0]]
        else:
            df = df0[df0['sample'] == sample_ids[0]]

    else:
        df = df0[df0['sample'] == sample_ids[0]]

    # =========================
    # Score runs by metrics
    # =========================

    scores = pd.DataFrame({'run': df['run']})
    total_score = pd.Series(0.0, index=df.index)
    
    any_metric_used = False

    for i, metric in enumerate(metrics):
        if metric not in df.columns or df[metric].isna().all():
            continue

        values = df[metric]
        min_v, max_v = values.min(), values.max()

        if math.isclose(min_v, max_v):
            print(f"Warning: metric '{metric}' has no variance, skipping", file=sys.stderr)
            continue

        # Min–max scaling
        scaled = (values - min_v) / (max_v - min_v)
        any_metric_used = True

        # Direction
        direction = directions[i]
        if direction == "-":
            scaled = 1.0 - scaled

        # Weight
        weight = weights[i]
        total_score += weight * scaled

        scores[metric + "_scaled"] = scaled

    scores["total_score"] = total_score

    if not any_metric_used:
        print(json.dumps([]))
        return 0

    best_run = scores.sort_values("total_score", ascending=False).iloc[0]["run"]
    print(json.dumps([best_run]))

    # =========================
    # Write scoring report
    # =========================
    if any_metric_used:
        with open("report.txt", "w") as f:
            f.write("Run optimisation report\n")
            f.write("=" * 60 + "\n\n")

            f.write("Scoring method:\n")
            f.write("- Metrics are min–max scaled to [0, 1]\n")
            f.write("- Higher scaled value = better performance\n")
            f.write("- Negatively associated metrics are inverted\n")
            f.write("- Final score = weighted sum of scaled metrics\n\n")

            f.write(f"Optimal run (highest total score): {best_run}\n\n")

            # -------------------------
            # Metric summaries
            # -------------------------
            for i, metric in enumerate(metrics):
                if metric not in df.columns or metric + "_scaled" not in scores.columns:
                    continue

                direction = "higher is better" if directions[i] == "+" else "lower is better"
                weight = weights[i]

                f.write(f"Top 10 runs by {metric}:\n")
                f.write("-" * 40 + "\n")
                f.write(f"Direction: {direction}, Weight: {weight}\n\n")

                metric_df = (
                    pd.DataFrame({
                        "run": df["run"],
                        "raw_value": df[metric],
                        "scaled_value": scores[metric + "_scaled"],
                    })
                    .assign(rank=df[metric].rank(ascending=(directions[i] == "-"), method="average", na_option="bottom"))
                    .sort_values("scaled_value", ascending=False)
                    .head(10)
                )

                for _, row in metric_df.iterrows():
                    f.write(
                        f"Run: {row['run']}, "
                        f"{metric}: {row['raw_value']:.6g}, "
                        f"scaled: {row['scaled_value']:.4f}, "
                        f"rank: {row['rank']}\n"
                    )

                f.write("\n")

            # -------------------------
            # Overall score ranking
            # -------------------------
            f.write("Top 10 runs by total score:\n")
            f.write("-" * 40 + "\n")

            top_total = scores.sort_values("total_score", ascending=False).head(10)

            for _, row in top_total.iterrows():
                f.write(
                    f"Run: {row['run']}, "
                    f"total_score: {row['total_score']:.4f}\n"
                )

            f.write("\n")

            # -------------------------
            # Per-run breakdown (best run)
            # -------------------------
            f.write(f"Score breakdown for optimal run ({best_run}):\n")
            f.write("-" * 40 + "\n")

            best_idx = scores[scores["run"] == best_run].index[0]

            for i, metric in enumerate(metrics):
                col = metric + "_scaled"
                if col not in scores.columns:
                    continue

                f.write(
                    f"{metric}: "
                    f"raw={df.loc[best_idx, metric]:.6g}, "
                    f"scaled={scores.loc[best_idx, col]:.4f}, "
                    f"weight={weights[i]}\n"
                )

if __name__ == "__main__":
    sys.exit(main())
