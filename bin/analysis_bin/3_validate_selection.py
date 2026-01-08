#!/usr/bin/env python3
import os
import pandas as pd
import traceback

# ------------------------
# Paths
# ------------------------

file_dir = "/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output"
metadata_file = "/scratch/shire/data/nj/projects/trimscreen_manuscript/metadata/table_1_metadata.csv"

# ------------------------
# Load inputs
# ------------------------

selected = pd.read_csv(
    os.path.join(file_dir, "selected_runs.tsv"),
    sep="\t"
)

metadata = pd.read_csv(metadata_file)

# ------------------------
# Sanity checks (global)
# ------------------------

required_selected_cols = {"dataset", "selected_run"}
if not required_selected_cols.issubset(selected.columns):
    raise ValueError(f"selected_runs.tsv must contain columns: {required_selected_cols}")

required_meta_cols = {"dataset_label", "dataset_name_in_pipeline"}
if not required_meta_cols.issubset(metadata.columns):
    raise ValueError(f"Metadata must contain columns: {required_meta_cols}")

# ------------------------
# Dataset label → pipeline name map
# ------------------------

ds_to_pipeline = dict(
    zip(metadata["dataset_label"], metadata["dataset_name_in_pipeline"])
)

# ------------------------
# Per-dataset evaluation
# ------------------------

results = []

for _, row in selected.iterrows():
    dataset = row["dataset"]
    sel_run = row["selected_run"]

    try:
        # Map dataset name
        if dataset not in ds_to_pipeline:
            raise ValueError("dataset_name_in_pipeline not found in metadata")

        dataset_name_in_pipeline = ds_to_pipeline[dataset]

        merged_file = os.path.join(
            file_dir,
            f"merged_table.{dataset_name_in_pipeline}.csv"
        )

        if not os.path.exists(merged_file):
            raise FileNotFoundError(f"Missing merged file: {merged_file}")

        # Load ALL runs for this dataset
        merged = pd.read_csv(merged_file)

        required_cols = {"run_id", "f1_score"}
        if not required_cols.issubset(merged.columns):
            raise ValueError(
                f"{merged_file} must contain columns: {required_cols}"
            )

        # Best possible F1
        best_f1 = merged["f1_score"].max()

        # Selected run row
        sel_row = merged[merged["run_id"] == sel_run]
        if sel_row.empty:
            raise ValueError(
                f"Selected run {sel_run} not found in merged table"
            )

        sel_f1 = sel_row["f1_score"].iloc[0]

        # Rank selected run by F1 (1 = best)
        merged = merged.copy()
        merged["f1_rank"] = merged["f1_score"].rank(
            ascending=False,
            method="min"
        )

        sel_rank = int(
            merged.loc[merged["run_id"] == sel_run, "f1_rank"].iloc[0]
        )

        results.append({
            "dataset": dataset,
            "selected_run": sel_run,
            "selected_f1": sel_f1,
            "max_f1": best_f1,
            "f1_gap": best_f1 - sel_f1,
            "relative_efficiency": sel_f1 / best_f1 if best_f1 > 0 else float("nan"),
            "f1_rank": sel_rank,
            "n_runs": len(merged)
        })

    except Exception as e:
        print(f"[WARN] Skipping dataset '{dataset}': {e}")
        # Uncomment for full traceback if debugging:
        # traceback.print_exc()
        continue

# ------------------------
# Output tables
# ------------------------

res = pd.DataFrame(results)

if res.empty:
    print("[ERROR] No datasets processed successfully. No output written.")
else:
    res.sort_values("dataset").to_csv(
        os.path.join(file_dir, "validation_per_dataset.tsv"),
        sep="\t",
        index=False,
        float_format="%.4f"
    )

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

