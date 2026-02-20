#!/usr/bin/env python
import pandas as pd
import sys
import math
import argparse
import json

#######################
# This script scores runs per sample using provided metrics,
# then aggregates scores across samples using:
#   mean(total_score) - λ * std(total_score)
# where λ = 0.5.
#######################

LAMBDA = 0.5  # variance penalty weight

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Find the best run across samples")
    parser.add_argument("-i", "--input", required=True, help="sample summary table csv")
    parser.add_argument("-m", "--metadata", required=False,
                        help="Optional metadata tsv table (currently not used for filtering)")
    parser.add_argument("--metrics", required=True, nargs="+",
                        help="Trimming-dependent metrics on which runs are evaluated and ranked")
    parser.add_argument("--metric_directions", nargs="+",
                        help="Metric directions: + for higher-is-better, - for lower-is-better")
    parser.add_argument("--metric_weights", nargs="+", type=float,
                        help="Optional metric weights (default = 1.0)")

    return parser.parse_args()


def score_runs_within_sample(df_sample, metrics, directions, weights):
    """
    Score runs within a single sample.
    Scaling is performed only across runs of this sample.
    """
    scores = pd.DataFrame({
        "sample": df_sample["sample"],
        "run": df_sample["run"],
    })

    total_score = pd.Series(0.0, index=df_sample.index)
    any_metric_used = False

    for i, metric in enumerate(metrics):
        if metric not in df_sample.columns or df_sample[metric].isna().all():
            continue

        values = df_sample[metric]
        min_v, max_v = values.min(), values.max()

        if math.isclose(min_v, max_v):
            # no variance → metric carries no information for this sample
            continue

        scaled = (values - min_v) / (max_v - min_v)
        any_metric_used = True

        if directions[i] == "-":
            scaled = 1.0 - scaled

        weight = weights[i]
        total_score += weight * scaled

        scores[f"{metric}_scaled"] = scaled

    scores["total_score"] = total_score

    if not any_metric_used:
        return None

    return scores


def main():
    args = parse_args()

    df0 = pd.read_csv(args.input)

    metrics = args.metrics
    directions = args.metric_directions
    weights = args.metric_weights or [1.0] * len(metrics)

    if not (len(metrics) == len(directions) == len(weights)):
        raise ValueError("metrics, metric_directions, and metric_weights must have the same length")

    # -------------------------
    # Score runs per sample
    # -------------------------
    per_sample_scores = []

    for sample_id, df_sample in df0.groupby("sample"):
        scored = score_runs_within_sample(
            df_sample, metrics, directions, weights
        )
        if scored is not None:
            per_sample_scores.append(scored)

    if not per_sample_scores:
        print(json.dumps([]))
        return 0

    scores_all = pd.concat(per_sample_scores, ignore_index=True)

    # -------------------------
    # Aggregate across samples
    # -------------------------
    run_stats = (
        scores_all
        .groupby("run")["total_score"]
        .agg(["mean", "std", "count"])
        .rename(columns={
            "mean": "mean_score",
            "std": "std_score",
            "count": "n_samples"
        })
        .fillna({"std_score": 0.0})
    )

    run_stats["final_score"] = (
        run_stats["mean_score"] - LAMBDA * run_stats["std_score"]
    )

    run_stats = run_stats.sort_values("final_score", ascending=False)

    best_run = run_stats.index[0]
    print(json.dumps([best_run]))

    # =========================
    # Write scoring report
    # =========================
    with open("report.txt", "w") as f:
        f.write("Run optimisation report (multi-sample)\n")
        f.write("=" * 60 + "\n\n")

        f.write("Scoring method:\n")
        f.write("- Metrics min–max scaled per sample across runs\n")
        f.write("- Negatively associated metrics inverted\n")
        f.write("- Final per-sample score = weighted sum of scaled metrics\n")
        f.write(
            f"- Final run score = mean(total_score) - {LAMBDA} × sd(total_score)\n\n"
        )

        f.write("Evaluated metrics:\n")
        for i, metric in enumerate(metrics):
            if metric not in df0.columns:
                continue

            direction = "higher is better" if directions[i] == "+" else "lower is better"
            weight = weights[i]

            f.write(f"{metric} (direction: {direction}, weight: {weight})\n")

        f.write(f"\nSuggested run: {best_run}\n\n")

        # -------------------------
        # Overall ranking
        # -------------------------
        f.write("Top 10 runs by final score:\n")
        f.write("-" * 40 + "\n")

        for run, row in run_stats.head(10).iterrows():
            f.write(
                f"Run: {run}, "
                f"final_score: {row['final_score']:.4f}, "
                f"mean: {row['mean_score']:.4f}, "
                f"sd: {row['std_score']:.4f}, "
                f"samples: {int(row['n_samples'])}\n"
            )

        f.write("\n")

        # -------------------------
        # Per-sample breakdown for best run
        # -------------------------
        f.write(f"Per-sample scores for suggested run ({best_run}):\n")
        f.write("-" * 40 + "\n")

        best_run_samples = (
            scores_all[scores_all["run"] == best_run]
            .sort_values("total_score", ascending=False)
        )

        for _, row in best_run_samples.iterrows():
            f.write(
                f"Sample: {row['sample']}, "
                f"total_score: {row['total_score']:.4f}\n"
            )

        f.write("\n")


if __name__ == "__main__":
    sys.exit(main())

