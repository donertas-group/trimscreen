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
        description="Compare runs. Create table from multiple files."
    )
    parser.add_argument(
        "-i", "--input", nargs="+", required=True, help="Input files"
    )
    parser.add_argument(
        "-m",
        "--metadata",
        required=False,
        help="Optional metadata tsv table with same format as required by nf-core/ampliseq",
    )
    return parser.parse_args()

def ruzicka_distance(x, y):
    num = np.minimum(x, y).sum()
    den = np.maximum(x, y).sum()
    if den == 0:
        return 0.0  # identical zero vectors
    return 1.0 - (num / den)

def calculate_mean_similarity(asv_table, rep_samples):
    """
    asv_table: ASVs x samples count table
    rep_samples: list of replicate sample IDs (columns)
    """

    if len(rep_samples) < 2:
        return None

    # Subset to replicates and drop ASVs with zero counts across all reps
    data = asv_table[rep_samples]
    data = data.loc[data.sum(axis=1) > 0]

    # calculate relative abundance
    data = data.div(data.sum(axis=0), axis=1)
    
    # Transpose: samples as rows, ASVs as columns
    data_T = data.T

    # Pairwise Ruzicka distances
    dist_matrix = pdist(data_T.values, metric=ruzicka_distance)

    # Convert to similarity
    similarities = 1.0 - dist_matrix

    return similarities.mean()


def process_run(
    summary_file,
    asv_file,
    tax_file,
    run,
    classifier_dir,
    ranks,
    rep_samples,
):
    """
    Process a single run directory to summarize ASV counts at specified taxonomic ranks.

    Returns:
    - A DataFrame summarizing the proportion of identified ASVs and the number of
      unique taxa at specified taxonomic ranks.
    """

    # Normalize rep_samples
    rep_samples = [] if rep_samples is None else list(rep_samples)

    if (
        os.path.exists(asv_file)
        and os.path.exists(tax_file)
        and os.path.exists(summary_file)
    ):
        asv_tax = pd.read_csv(tax_file, sep="\t", index_col=0)
        asv_table = pd.read_csv(asv_file, sep="\t", index_col=0)
        summary_table = pd.read_csv(summary_file, sep="\t", index_col=0)
    else:
        logging.warning(
            f"{classifier_dir}, ASV_table.tsv.gz, ASV_tax_*.tsv.gz "
            f"or overall_summary.tsv.gz not found in {run}"
        )
        return

    # Merge ASV and taxonomy tables
    merged_table = asv_table.merge(asv_tax, left_index=True, right_index=True)

    results = {}

    # Diversity metrics per rank
    for rank in ranks:
        if rank not in merged_table.columns:
            logging.warning(f"Rank '{rank}' not found in the table. Skipping.")
            continue

        ntaxa = []
        nasv = []
        shannons = []
        simpsons = []

        for sample in asv_table.columns:
            valid_asvs = merged_table[
                (merged_table[sample] > 0) & (merged_table[rank].notna()) & (merged_table[rank].astype(str).str.strip() != "")
            ]

            nasv.append(valid_asvs[sample].count())
            ntaxa.append(valid_asvs[rank].nunique())

            rank_nreads = valid_asvs.groupby(rank)[sample].sum()

            if rank_nreads.sum() > 0 and len(rank_nreads) > 1:
                shannons.append(np.exp(diversity.alpha.shannon(rank_nreads)))
                simpsons.append(np.exp(diversity.alpha.simpson(rank_nreads)))
            else:
                shannons.append(np.nan)
                simpsons.append(np.nan)

        results[rank] = ntaxa
        results[f"{rank}_nasv"] = nasv
        results[f"shannon_{rank}"] = shannons
        results[f"simpson_{rank}"] = simpsons

    # Per-sample metrics
    Nasvs = []
    Nreads = []

    samples = asv_table.columns.unique()

    for sample in samples:
        Nasvs.append(asv_table.loc[asv_table[sample] > 0, sample].count())
        Nreads.append(asv_table.loc[asv_table[sample] > 0, sample].sum())

    # Replicate similarity (Genus only, intentional)
    rep_sim = None

    replicates = list(set(rep_samples).intersection(samples))
  
    if len(replicates) >= 2:
        rep_sim = calculate_mean_similarity(asv_table, replicates)

    rep_similarity = [
        rep_sim if sample in replicates else None
        for sample in samples
    ]

    results["nasvs"] = Nasvs
    results["nreads"] = Nreads
    results["nasvs_in_run"] = (asv_table > 0).sum().tolist()
    results["rep_similarity"] = rep_similarity

    res_df = pd.DataFrame(results, index=asv_table.columns)
    res_df["run"] = run

    # Summary table metrics
    if "lenfilter_output" in summary_table.columns:
        res_df["DADA2_input"] = summary_table["DADA2_input"].reindex(res_df.index)

        res_df["retained_reads_percent"] = (
            summary_table["lenfilter_output"].reindex(res_df.index)
            / res_df["DADA2_input"].replace(0, np.nan)
        )
    else:
        res_df["DADA2_input"] = summary_table["DADA2_input"].reindex(res_df.index)
        res_df["retained_reads_percent"] = (
            summary_table["nonchim"].reindex(res_df.index)
            / res_df["DADA2_input"].replace(0, np.nan)
        )

    # Reorder columns (no 'sample' column; it's the index)
    res_df = res_df[["run"] + [c for c in res_df.columns if c != "run"]]

    return res_df


def main():
    args = parse_args()

    inputs = args.input
    metadata_csv = args.metadata

    rep_samples = []

    if metadata_csv is not None:
        metadata = pd.read_csv(metadata_csv)
        if "is_replicate" in metadata.columns:
            rep_samples = metadata.loc[
                metadata["is_replicate"]
                .astype(str)
                .str.lower()
                == "true",
                "sampleID",
            ].tolist()

    classifier_dir = "dada2"
    Ranks_to_analyse = ["Phylum", "Family", "Genus", "Species"]

    summary = process_run(
        inputs[0],
        inputs[1],
        inputs[2],
        inputs[3],
        classifier_dir,
        Ranks_to_analyse,
        rep_samples,
    )

    summary.index.name = "sample"
    summary.to_csv(f"{inputs[3]}_table.csv", sep=",", index=True)


if __name__ == "__main__":
    sys.exit(main())

