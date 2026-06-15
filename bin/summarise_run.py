#!/usr/bin/env python

import argparse
import sys
import pandas as pd
import numpy as np
import os
from skbio import diversity
from scipy.spatial.distance import pdist, squareform
import logging


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Compare runs and compute per-sample and per-run summaries including beta diversity.")

    parser.add_argument("-i", "--input", nargs=4, required=True,
                        help="summary_file asv_file tax_file run_name")

    parser.add_argument("-m", "--metadata", required=False, 
                        help="Metadata CSV with columns: sampleID, replicated, bio_sample")

    return parser.parse_args()


def clr_transform(df, pseudocount=1e-6):
    df = df + pseudocount
    log_df = np.log(df)
    gm = log_df.mean(axis=0)
    clr = log_df.subtract(gm, axis=1)
    return clr

def calculate_beta_diversity(asv_table, metadata):

    samples = asv_table.columns.tolist()

    # Keep only samples present in metadata
    metadata = metadata.set_index("sampleID")
    metadata = metadata.loc[metadata.index.intersection(samples)]

    if len(metadata) < 2:
        return None  # too few samples to compute any distances

    # Determine biological groups
    metadata["replicated"] = metadata["replicated"].astype(str).str.lower() == "true"

    metadata["bio_group"] = np.where(
        metadata["replicated"],
        metadata["bio_sample"],
        metadata.index
    )

    # Count groups
    group_counts = metadata["bio_group"].value_counts()

    n_samples = len(metadata)
    n_bio_groups = group_counts.shape[0]
    n_bio_samples = metadata["bio_sample"].nunique()

    # Conditions for each metric
    has_within = (group_counts >= 2).any()
    has_between = n_bio_groups >= 2

    # Subset ASV table
    asv_table = asv_table[metadata.index]

    # CLR transform
    clr_table = clr_transform(asv_table)

    # Distance matrix
    dist_df = pd.DataFrame(
        squareform(pdist(clr_table.T, metric="euclidean")),
        index=clr_table.columns,
        columns=clr_table.columns,
    )

    within = []
    between = []

    cols = dist_df.index.tolist()

    for i, s1 in enumerate(cols):
        g1 = metadata.loc[s1, "bio_group"]

        for j in range(i + 1, len(cols)):
            s2 = cols[j]
            g2 = metadata.loc[s2, "bio_group"]

            d = dist_df.loc[s1, s2]

            if g1 == g2:
                within.append(d)
            else:
                between.append(d)

    # Compute metrics independently
    mean_within = np.mean(within) if has_within and within else np.nan
    mean_between = np.mean(between) if has_between and between else np.nan

    # Ratio only if both exist
    if has_within and has_between and mean_within > 0:
        ratio = mean_between / mean_within
    else:
        ratio = np.nan

    return {
        "n_samples": n_samples,
        "n_bio_samples": int(n_bio_samples),
        #"n_biologically_different_samples": n_bio_groups,
        "mean_within_replicate_dist": mean_within,
        "mean_between_sample_dist": mean_between,
        "between_within_ratio": ratio,
    }



# -----------------------------
# Main processing
# -----------------------------

def process_run(summary_file, asv_file, tax_file, run, metadata):

    if not (os.path.exists(asv_file)
            and os.path.exists(tax_file)
            and os.path.exists(summary_file)):
        logging.warning(f"Missing files for run {run}")
        return None, None

    asv_tax = pd.read_csv(tax_file, sep="\t", index_col=0)
    asv_table = pd.read_csv(asv_file, sep="\t", index_col=0)
    summary_table = pd.read_csv(summary_file, sep="\t", index_col=0)

    merged_table = asv_table.merge(asv_tax, left_index=True, right_index=True)

    ranks = ["Phylum", "Family", "Genus", "Species"]
    sample_results = {}

    # -----------------------------
    # Alpha diversity per rank
    # -----------------------------

    for rank in ranks:
        if rank not in merged_table.columns:
            continue

        shannons = []

        for sample in asv_table.columns:
            valid = merged_table[
                (merged_table[sample] > 0) &
                (merged_table[rank].notna()) &
                (merged_table[rank].astype(str).str.strip() != "")
            ]

            rank_counts = valid.groupby(rank)[sample].sum()

            if rank_counts.sum() > 0 and len(rank_counts) > 1:
                sh = np.exp(diversity.alpha.shannon(rank_counts))
            else:
                sh = np.nan

            shannons.append(sh)

        sample_results[f"shannon_{rank}"] = shannons

    # -----------------------------
    # Per-sample metrics
    # -----------------------------

    sample_results["nasvs"] = (
        (asv_table > 0).sum(axis=0).tolist()
    )

    sample_results["nreads"] = (
        asv_table.sum(axis=0).tolist()
    )

    sam_df = pd.DataFrame(sample_results, index=asv_table.columns)
    sam_df["run"] = run

    # Add DADA2 columns (remain per-sample for now)
    sam_df["DADA2_input"] = summary_table["DADA2_input"].reindex(sam_df.index)

    if "lenfilter_output" in summary_table.columns:
        retained = summary_table["lenfilter_output"]
    else:
        retained = summary_table["nonchim"]

    sam_df["retained_reads_percent"] = (
        retained.reindex(sam_df.index)
        / sam_df["DADA2_input"].replace(0, np.nan)
    )

    # -----------------------------
    # Run-level summaries
    # -----------------------------

    run_summary = {
        "run": run,
        "nasvs": (asv_table > 0).any(axis=1).sum(),
        "DADA2_input": sam_df["DADA2_input"].sum(),
        "retained_reads_percent": retained.sum() / sam_df["DADA2_input"].replace(0, np.nan).sum()
    }

    # Beta diversity only if metadata provided
    if metadata is not None:
        beta_metrics = calculate_beta_diversity(asv_table, metadata)
        if beta_metrics is not None:
            run_summary.update(beta_metrics)

    # Mean Shannon across ranks
    for rank in ranks:
        if rank not in merged_table.columns:
            continue

        run_summary[f"mean_Shannon_{rank}"] = sam_df[f"shannon_{rank}"].mean()

    run_df = pd.DataFrame([run_summary])

    return sam_df, run_df


# -----------------------------
# Main
# -----------------------------

def main():
    args = parse_args()

    summary_file, asv_file, tax_file, run = args.input

    metadata = None
    if args.metadata is not None:
        metadata = pd.read_csv(args.metadata)

    sample_df, run_df = process_run(
        summary_file,
        asv_file,
        tax_file,
        run,
        metadata
    )

    if sample_df is None:
        sys.exit(1)

    sample_df.index.name = "sample"

    sample_df.to_csv(f"samplerun_summary.csv", index=True)
    run_df.to_csv(f"run_summary.csv", index=False)


if __name__ == "__main__":
    sys.exit(main())
