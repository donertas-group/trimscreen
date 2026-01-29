#!/usr/bin/env python
# example usage:
# ./betadiv.py -D mock16 
# ./betadiv.py -D schirmer2015 --out_suffix .2

import os
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import argparse
import sys

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Calculate beta-diversity for mock dataset (per-sample, per-run)")
    parser.add_argument("-D", "--dataset_name", required=True, help="Dataset name")
    parser.add_argument("--out_suffix", default="", help="suffix string of pipeline output dir")
    parser.add_argument("--out", default="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output", help="Output dir")
    return parser.parse_args()

def clr_transform(df, pseudocount=1e-6):
    """Centered log-ratio (CLR) transformation per column (sample)."""
    df = df + pseudocount
    log_df = np.log(df)
    gm = log_df.mean(axis=0)
    clr = log_df.subtract(gm, axis=1)
    return clr

def main():
    args = parse_args()

    base_dir = f"/scratch/shire/ssd/pipeline/16s_nf_pipeline/{args.dataset_name}/output{args.out_suffix}"
    compare_runs_csv = os.path.join(base_dir, "compare_runs", "filtered_table.csv")

    # 1. Get run list
    df_runs = pd.read_csv(compare_runs_csv)
    run_list = df_runs['run'].dropna().unique().tolist()

    if not run_list:
        print("No runs found in filtered_table.csv")
        sys.exit(1)

    # 2. Load ASV tables
    asv_tables = []
    for run in run_list:
        path = os.path.join(base_dir, "runs", run, "dada2", "ASV_table.tsv.gz")
        if not os.path.exists(path):
            print(f"ASV table not found for run {run}: {path}")
            continue
        df = pd.read_csv(path, sep="\t")
        df = df.set_index("ASV_ID")
        df.columns = [f"{run}|{c}" for c in df.columns]  # mark both run and sample
        asv_tables.append(df)

    if not asv_tables:
        print("No ASV tables loaded")
        sys.exit(1)

    # 3. Merge ASV tables
    asv_table_full = asv_tables[0]
    for df in asv_tables[1:]:
        asv_table_full = asv_table_full.merge(df, left_index=True, right_index=True, how="outer")
    asv_table_full = asv_table_full.fillna(0)

    # 4. Extract sample/run names
    sample_names = sorted({c.split("|")[1] for c in asv_table_full.columns})
    run_names = sorted({c.split("|")[0] for c in asv_table_full.columns})

    results = []

    # 5. Compute per-run, per-sample mean Aitchison distance
    for run in run_names:
        run_cols = [c for c in asv_table_full.columns if c.startswith(f"{run}|")]
        run_table = asv_table_full[run_cols]

        if run_table.shape[1] < 2:
            continue  # skip runs with <2 samples

        # CLR transform
        clr_table = clr_transform(run_table)

        # Compute Euclidean distance on CLR = Aitchison distance
        dist_matrix = pd.DataFrame(
            squareform(pdist(clr_table.T, metric="euclidean")),
            index=clr_table.columns,
            columns=clr_table.columns
        )

        # Reduce to per-sample metric: mean distance to all other samples in the run
        for sample_col in clr_table.columns:
            sample_name = sample_col.split("|")[1]
            other_cols = [c for c in clr_table.columns if c != sample_col]
            mean_dist = dist_matrix.loc[sample_col, other_cols].mean()
            results.append({
                "run": run,
                "sample": sample_name,
                "mean_dist": mean_dist
            })

    if not results:
        print("No distances computed")
        sys.exit(1)

    df_results = pd.DataFrame(results)

    # 6. Save results
    out_path = os.path.join(args.out, f"{args.dataset_name}_per_sample_beta_diversity.csv")
    os.makedirs(args.out, exist_ok=True)
    df_results.to_csv(out_path, index=False)
    print(f"Per-sample beta-diversity metrics saved to {out_path}")


if __name__ == "__main__":
    sys.exit(main())

