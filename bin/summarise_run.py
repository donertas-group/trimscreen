#!/usr/bin/env python

import argparse
import sys
import pandas as pd
import numpy as np
import os
from skbio import diversity
from scipy.spatial.distance import pdist, squareform


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


def ruzicka(u, v):
    """Compute the Ruzicka (abundance-based Jaccard) distance between two vectors."""
    u, v = np.asarray(u), np.asarray(v)
    denominator = np.maximum(u, v).sum()
    if denominator == 0:
        return 0.0
    return 1.0 - (np.minimum(u, v).sum() / denominator)


def calculate_mean_similarity(asv_table, asv_tax, run_id, rep_samples, rank):
    merged = asv_table.merge(asv_tax, left_index=True, right_index=True)

    # Drop NA or empty rank entries
    merged = merged[
        merged[rank].notna() & (merged[rank].astype(str).str.strip() != "")
    ]

    if merged.empty or len(rep_samples) < 2:
        return None

    # Aggregate counts at chosen rank
    agg = merged.groupby(rank)[rep_samples].sum()

    # Transpose for distance computation (samples as rows)
    agg_T = agg.T

    # Compute Ruzicka distance matrix
    dist_matrix = pdist(agg_T, metric=ruzicka)
    dist_df = pd.DataFrame(
        squareform(dist_matrix),
        index=agg_T.index,
        columns=agg_T.index,
    )

    # Mean pairwise distance → similarity
    mean_dist = dist_df.values[np.triu_indices_from(dist_df, k=1)].mean()
    similarity = 1.0 - mean_dist

    return similarity


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
        print(
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
            print(f"Rank '{rank}' not found in the table. Skipping.")
            continue

        ntaxa = []
        nasvs = []
        shannons = []
        simpsons = []

        for sample in asv_table.columns:
            valid_asvs = merged_table[
                (merged_table[sample] > 0) & (merged_table[rank].notna())
            ]

            nasvs.append(valid_asvs[sample].count())
            ntaxa.append(valid_asvs[rank].nunique())

            rank_nasvs = valid_asvs.groupby(rank)[sample].count()

            if rank_nasvs.sum() > 0 and len(rank_nasvs) > 1:
                shannons.append(diversity.alpha.shannon(rank_nasvs))
                simpsons.append(diversity.alpha.simpson(rank_nasvs))
            else:
                shannons.append(np.nan)
                simpsons.append(np.nan)

        results[rank] = ntaxa
        results[f"{rank}_nasv"] = nasvs
        results[f"shannon_{rank}"] = shannons
        results[f"simpson_{rank}"] = simpsons

    # Per-sample metrics
    Nasvs = []
    Nreads = []

    for sample in asv_table.columns:
        Nasvs.append(asv_table.loc[asv_table[sample] > 0, sample].count())
        Nreads.append(asv_table.loc[asv_table[sample] > 0, sample].sum())

    # Replicate similarity (Genus only, intentional)
    rep_sim = None
    if len(rep_samples) >= 2:
        rep_sim = calculate_mean_similarity(
            asv_table, asv_tax, run, rep_samples, rank="Genus"
        )

    rep_similarity = [
        rep_sim if sample in rep_samples else None
        for sample in asv_table.columns
    ]

    results["nasvs"] = Nasvs
    results["nreads"] = Nreads
    results["nasvs_in_run"] = (
        asv_table.count().reindex(asv_table.columns).tolist()
    )
    results["rep_similarity"] = rep_similarity

    res_df = pd.DataFrame(results, index=asv_table.columns)
    res_df["run"] = run

    # Summary table metrics
    if "lenfilter_output" in summary_table.columns:
        res_df["DADA2_input"] = summary_table["DADA2_input"]
        res_df["retained_reads_percent"] = (
            summary_table["lenfilter_output"]
            / summary_table["DADA2_input"].replace(0, np.nan)
        )
    else:
        res_df["retained_reads_percent"] = (
            summary_table["nonchim"]
            / summary_table["DADA2_input"].replace(0, np.nan)
        )

    for rank in ranks:
        res_df[f"{rank}_pasv"] = (
            res_df[f"{rank}_nasv"]
            / res_df["nasvs"].replace(0, np.nan)
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

