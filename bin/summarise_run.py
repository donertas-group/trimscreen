#!/usr/bin/env python
import argparse
import sys
import subprocess
import pandas as pd
import numpy as np
import os
from skbio import diversity

def parse_args(args=None):

    parser = argparse.ArgumentParser(description="Compare runs. Create table from multiple files.")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input files")

    return parser.parse_args()


def ruzicka(u, v):
    """Compute the Ruzicka (abundance-based Jaccard) similarity between two vectors."""
    u, v = np.asarray(u), np.asarray(v)
    numerator = np.minimum(u, v).sum()
    denominator = np.maximum(u, v).sum()
    if denominator == 0:  # avoid division by zero
        return 0.0
    return 1-(numerator / denominator)


def calculate_mean_similarity(asv_table, asv_tax, run_id, rep_samples, rank):

    merged = pd.merge(asv_table, asv_tax, on='ASV_ID')

    # Drop NA or empty rank entries
    merged = merged[merged[rank].notna() & (merged[rank].astype(str).str.strip() != '')]

    if merged.empty:
        return None

    meta_cols = ["ASV_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "confidence","sequence"]
    meta_cols = [c for c in meta_cols if c in merged.columns]

    # After finding rep_samples
    if len(rep_samples) < 2:
        return None  # or return 1.0 (perfect similarity), or 0.0

    # Aggregate counts at chosen rank
    agg = merged.groupby(rank)[rep_samples].sum()

    # Transpose for distance computation
    agg_T = agg.T  # samples as rows

    # Compute distance
    dist_matrix = pdist(agg_T, metric=ruzicka)
    dist_df = pd.DataFrame(squareform(dist_matrix), index=agg_T.index, columns=agg_T.index)

    # Compute mean dissimilarity & similarity
    mean_dist = dist_df.values[np.triu_indices_from(dist_df, k=1)].mean()
    similarity = 1 - mean_dist

    return similarity


def process_run(summary_file, asv_file, tax_file, run, classifier_dir, ranks, rep_samples):
    """
    Process a single run directory to summarize ASV counts at specified taxonomic ranks.
    Parameters:
    - asv_file: 'ASV_table.tsv.gz'
    - tax_file: 'ASV_tax_*.tsv.gz'
    - summary_file: 'overall_summary.tsv.gz'
    - ranks: a list of taxonomic ranks to analyse (e.g., 'Phylum', 'Class', etc.).
    Returns:
    - A DataFrame summarizing the proporatino of id'ed ASVs and the number of unique taxa at specified taxonomic ranks
    """
    # Load ASV and taxonomy tables

    if os.path.exists(asv_file) and os.path.exists(tax_file) and os.path.exists(summary_file):
        # Read the files 
        asv_tax = pd.read_csv(tax_file, sep='\t', index_col=0) # assuming there's exactly one match
        asv_table = pd.read_csv(asv_file, sep='\t', index_col=0)
        summary_table = pd.read_csv(summary_file, sep='\t', index_col=0)
    else:
        print(f"{classifier_dir}, ASV_table.tsv.gz, ASV_tax_*.tsv.gz or overall_summary.tsv.gz not found in {run}")
        return

    # Merge ASV and taxonomy tables on index (ASV)
    merged_table = asv_table.merge(asv_tax, left_index=True, right_index=True)

    # Create an empty dictionary to store results
    results = {}

    # Iterate through each rank
    for rank in ranks:
        if rank not in merged_table.columns:
            print(f"Rank '{rank}' not found in the table. Skipping.")
            continue

        # Initialize lists to store diversity metrics for each sample
        ntaxa = []
        nasvs = []
        shannons = []
        simpsons = []
 
        # Iterate through each sample column
        for sample in asv_table.columns.tolist():
            # Filter for ASVs with counts > 0 and corresponding taxonomic name not NaN
            valid_asvs = merged_table[(merged_table[sample] > 0) & (merged_table[rank].notna())]
            nasvs.append(valid_asvs[sample].count())

            # Calculate richness and alpha diversity at unique taxonomic identities
            unique_count = valid_asvs[rank].nunique()
            ntaxa.append(unique_count)

            # Group by phylum and sum the counts
            rank_nasvs = valid_asvs.groupby(rank)[sample].count()

            # Calculate Shannon diversity index using skbio
            shannon = diversity.alpha.shannon(rank_nasvs)
            shannons.append(shannon)

            simpson = diversity.alpha.simpson(rank_nasvs)
            simpsons.append(simpson)

        # Store the counts in the results dictionary
        results[rank] = ntaxa
        results[f'{rank}_nasv'] = nasvs
        results[f'shannon_{rank}'] = shannons
        results[f'simpson_{rank}'] = simpsons
        

    Nasvs = []
    Nreads = []
    rep_similarity = []

    for sample in asv_table.columns.tolist():
        Nasvs.append(asv_table[(asv_table[sample] > 0)][sample].count())
        Nreads.append(asv_table[(asv_table[sample] > 0)][sample].sum())

        if sample in rep_samples 
            similarity = calculate_mean_simiarlity(asv_table, asv_tax, run, rep_samples, rank="Genus")
            rep_similarity.append(similarity)
        else rep_similarity.append(None)

    results['nasvs'] = Nasvs
    results['nreads'] = Nreads
    results['nasvs_in_run'] = asv_table.count()
    results['rep_similarity'] = rep_similarity

    # Create a DataFrame from results
    res_df = pd.DataFrame(results, index=asv_table.columns.tolist())  # Use sample names as index
    res_df['run'] = run

    # Calculate and add more metrics from summary_table
    if {'lenfilter_output'}.issubset(summary_table.columns):
        res_df['DADA2_input'] = summary_table['DADA2_input']
        res_df['retained_reads_percent'] = summary_table['lenfilter_output']/ summary_table['DADA2_input'].replace(0, np.nan)
    else:
        res_df['retained_reads_percent'] = summary_table['nonchim'] / summary_table['DADA2_input'].replace(0, np.nan)

    for rank in ranks:
        #res_df[f'{rank}_pread'] = res_df[f'{rank}_nread'] / res_df['input_tax_filter'].replace(0, np.nan)
        res_df[f'{rank}_pasv'] = res_df[f'{rank}_nasv'] / res_df['nasvs'].replace(0, np.nan)
    
    # reorder columns
    #res_df = res_df[['sample', 'run'] + [c for c in res_df.columns if c not in ['sample', 'run']]]

    return res_df


def main():
    args = parse_args()
    
    inputs = args.input

    classifier_dir = "dada2"
    Ranks_to_analyse = ["Phylum","Family","Genus","Species"]
#### modify!
    rep_samples = # read from metadata #
####

    summary = process_run(inputs[0], inputs[1], inputs[2], inputs[3], classifier_dir, Ranks_to_analyse, rep_samples)
    summary.index.name = "sample"
    summary.to_csv(f"{inputs[3]}_table.csv", sep=",", index=True)


if __name__ == "__main__":
    sys.exit(main())

