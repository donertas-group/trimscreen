#!/usr/bin/env python
# example usage:
# ./betadiv.py -D mock16 
# ./betadiv.py -D schirmer2015 --out_suffix .2

import os
import re
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import argparse
import sys

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Calculate beta-diversity for mock dataset (per-sample, per-run)")
    parser.add_argument("-D", "--dataset_name", required=True, help="Dataset name")
    parser.add_argument("-R", "--rank", default="Genus", help="Taxonomic rank to calculate beta")
    parser.add_argument("--out_suffix", default="", help="suffix string of pipeline output dir")
    parser.add_argument("--out", default="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output", help="Output dir")

    return parser.parse_args()

def ruzicka(u, v):
    """Compute the Ruzicka (abundance-based Jaccard) dissimilarity between two vectors."""
    u, v = np.asarray(u), np.asarray(v)
    numerator = np.minimum(u, v).sum()
    denominator = np.maximum(u, v).sum()
    if denominator == 0:
        return 0.0
    return 1 - (numerator / denominator)

def main():
    args = parse_args()
    rank = args.rank

    base_dir = f"/scratch/shire/ssd/pipeline/16s_nf_pipeline/{args.dataset_name}/output{args.out_suffix}"
    compare_runs_csv = os.path.join(base_dir, "compare_runs", "full_table.csv")

    # 1. Get run list
    df_runs = pd.read_csv(compare_runs_csv)
    run_list = df_runs['run'].dropna().unique().tolist()

    # 2. Load ASV tables
    asv_tables = []
    for run in run_list:
        path = os.path.join(base_dir, "runs", run, "dada2", "ASV_table.tsv.gz")
        df = pd.read_csv(path, sep="\t")
        df = df.set_index("ASV_ID")
        df.columns = [f"{run}|{c}" for c in df.columns]  # mark both run and sample
        asv_tables.append(df)

    # 3. Merge ASV tables
    asv_table_full = asv_tables[0]
    for df in asv_tables[1:]:
        asv_table_full = asv_table_full.merge(df, left_index=True, right_index=True, how="outer")
    asv_table_full = asv_table_full.fillna(0)

    # 4. Merge taxonomy
    tax_tables = []
    for run in run_list:
        path = os.path.join(base_dir, "runs", run, "dada2", "ASV_tax_species.silva_138_2.tsv.gz")
        df_tax = pd.read_csv(path, sep="\t")
        tax_tables.append(df_tax)
    asv_tax_full = pd.concat(tax_tables).drop_duplicates(subset="ASV_ID")
    asv_tax_full = asv_tax_full.set_index("ASV_ID")

    # Align indices
    common = asv_table_full.index.intersection(asv_tax_full.index)
    asv_table_full = asv_table_full.loc[common]
    asv_tax_full = asv_tax_full.loc[common]

    # 5. Aggregate by taxonomic rank
    merged = asv_table_full.join(asv_tax_full[[rank]])
    agg_table = merged.groupby(rank).sum()

    ### --- CHANGED: Instead of summing across samples in each run ---
    ### We now extract the sample names and compute distances per sample

    # Extract sample names (without run prefix)
    sample_names = sorted({c.split("|")[1] for c in agg_table.columns})

    run_names = sorted({c.split("|")[0] for c in agg_table.columns})

    results = []  # store per-sample pairwise distances

    for sample in sample_names:
        # Get columns for this sample across all runs
        sample_cols = [c for c in agg_table.columns if c.endswith(f"|{sample}")]
        if len(sample_cols) < 2:
            continue  # skip samples not found in multiple runs

        sub = agg_table[sample_cols]
        sub.columns = [c.split("|")[0] for c in sub.columns]  # rename to just run name

        # compute ruzicka distances between runs for this sample
        dist_matrix = pd.DataFrame(
            squareform(pdist(sub.T, metric=ruzicka)),
            index=sub.columns,
            columns=sub.columns
        )

        # melt into long form
        dist_long = dist_matrix.stack().reset_index()
        dist_long.columns = ["run1", "run2", "ruzicka_dissimilarity"]
        dist_long = dist_long[dist_long["run1"] < dist_long["run2"]]  # avoid duplicates
        dist_long["sample"] = sample
        results.append(dist_long)

    if not results:
        print("No samples found across multiple runs.")
        sys.exit(1)

    all_dists = pd.concat(results, ignore_index=True)

    # Save all pairwise per-sample distances
    all_outfile = os.path.join(args.out, f"ruzicka_per_sample.{args.dataset_name}{args.out_suffix}.csv")
    all_dists.to_csv(all_outfile, index=False)
    print(f"Per-sample Ruzicka distances saved to: {all_outfile}")

    ### --- NEW: summarize per-run median dissimilarity across all samples ---
    run_scores = []
    for run in run_names:
        # collect distances where this run participated
        subset = all_dists[(all_dists.run1 == run) | (all_dists.run2 == run)]
        if subset.empty:
            continue
        median_val = subset["ruzicka_dissimilarity"].median()
        run_scores.append({"run_id": run, "median_dissimilarity": median_val})

    median_df = pd.DataFrame(run_scores)
    median_outfile = os.path.join(args.out, f"median_dissimilarity_per_sample.{args.dataset_name}{args.out_suffix}.csv")
    median_df.to_csv(median_outfile, index=False)
    print(f"Median across-run dissimilarity saved to: {median_outfile}")

if __name__ == "__main__":
    sys.exit(main())

