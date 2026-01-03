#!/usr/bin/env python3
import os
import pandas as pd

# env
file_dir = "/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output" 

# ------------------------
# Load inputs
# ------------------------

selected = pd.read_csv(os.path.join(file_dir, "selected_runs.tsv"), sep="\t")
f1 = pd.read_csv(os.path.join(file_dir, "f1_scores.tsv"), sep="\t")

# Sanity checks
required_cols = {"dataset", "run", "f1"}
if not required_cols.issubset(f1.columns):
    raise ValueError(f"f1_scores.tsv must contain columns: {required_cols}")

# ------------------------
# Merge selection with F1
# ------------------------

df = f1.merge(
    selected,
    on="dataset",
    how="left",
    validate="many_to_one"
)

if df["selected_run"].isna().any():
    raise ValueError("Some datasets are missing selected runs")

# ------------------------
# Per-dataset evaluation
# ------------------------

results = []

for dataset, dfd in df.groupby("dataset"):
    dfd = dfd.copy()

    # Best possible F1
    idx_best = dfd["f1"].idxmax()
    best_f1 = dfd.loc[idx_best, "f1"]

    # Selected run
    sel_run = dfd["selected_run"].iloc[0]
    sel_row = dfd[dfd["run"] == sel_run]

    if sel_row.empty:
        raise ValueError(f"Selected run {sel_run} not found in F1 table for {dataset}")

    sel_f1 = sel_row["f1"].iloc[0]

    # Rank selected run by F1 (1 = best)
    dfd["f1_rank"] = dfd["f1"].rank(ascending=False, method="min")
    sel_rank = dfd.loc[dfd["run"] == sel_run, "f1_rank"].iloc[0]

    results.append({
        "dataset": dataset,
        "selected_run": sel_run,
        "selected_f1": sel_f1,
        "max_f1": best_f1,
        "f1_gap": best_f1 - sel_f1,
        "relative_efficiency": sel_f1 / best_f1 if best_f1 > 0 else float("nan"),
        "f1_rank": int(sel_rank),
        "n_runs": len(dfd)
    })

res = pd.DataFrame(results)

# ------------------------
# Output tables
# ------------------------

# Main manuscript / supplement table
res.sort_values("dataset").to_csv(
    os.path.join(file_dir, "validation_per_dataset.tsv"),
    sep="\t",
    index=False,
    float_format="%.4f"
)

# Summary statistics (useful for Results text)
summary = res[[
    "selected_f1",
    "max_f1",
    "f1_gap",
    "relative_efficiency",
    "f1_rank"
]].describe()

summary.to_csv(
    os.path.join(file_dir, "validation_summary_stats.tsv"),
    sep="\t",
    float_format="%.4f"
)

