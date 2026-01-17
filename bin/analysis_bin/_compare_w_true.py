#!/usr/bin/env python3

# The script expects the input files in a structured directory, as described.
# It uses read count as the x-axis and calculates detection stats and f1 summary score at the user-specified taxonomic rank.
# Make sure column names match exactly, especially for ranks (Genus, Family, etc.) — they are case-sensitive.

# example usage:
# ./compare_w_true.py -D mock13 -r run_180190 -R Genus (--true /path/to/true_composition.csv -X 20)

import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt
import glob
import json

def parse_args():
    parser = argparse.ArgumentParser(description="Compare 16S data with true mock composition.")
    parser.add_argument("-D", "--dataset_name", required=True, type=str, help="Mock community (e.g. mock02)")
    parser.add_argument("-r", "--run", required=True, type=str, help="Run ID")
    parser.add_argument("-R", "--rank", required=True, type=str, help="Taxonomic rank to compare (e.g., Genus)")
    parser.add_argument("--true", default="true_composition.csv", type=str, help="Path to true composition CSV file (full path or filename)")
    parser.add_argument("--out_suffix", default="", help="suffix string of pipeline output dir")
    parser.add_argument("-X", default=10, type=int, help="Target ASV read number filter value")
    return parser.parse_args()

def read_files(mock, run_id, true_comp_file, threshold, out_suffix=""):
    result_dir = "/scratch/shire/ssd/pipeline/16s_nf_pipeline"  # Modify if needed
    data_dir = "/scratch/shire/data/nj/raw_data/published"  # Modify if needed

    asv_tax_file = os.path.join(result_dir, f"{mock}", f"output{out_suffix}/runs", run_id, "dada2/ASV_tax.silva_138_2.tsv.gz")
    asv_table_file = os.path.join(result_dir, f"{mock}", f"output{out_suffix}/runs", run_id, "dada2/ASV_table.tsv.gz")
    summary_file = os.path.join(result_dir, f"{mock}", f"output{out_suffix}/runs", run_id, "overall_summary.tsv")

    # If user provided a full path, use it directly
    if os.path.isabs(true_comp_file):
        true_comp_path = true_comp_file
    else:
        # Otherwise use existing automatic directory structure
        if mock.startswith("mock"):
            true_comp_path = os.path.join(data_dir, "mockrobiota", mock, true_comp_file)
        else:
            true_comp_path = os.path.join(data_dir, mock, true_comp_file)

    asv_tax = pd.read_csv(asv_tax_file, sep='\t')
    asv_table = pd.read_csv(asv_table_file, sep='\t')

    # identify sample columns (all except ASV_ID)
    sample_cols = asv_table.columns.drop('ASV_ID')

    # compute max across samples and filter
    filtered_asv_table = asv_table[asv_table[sample_cols].max(axis=1) >= threshold]

    true_comp = pd.read_csv(true_comp_path)
    summary_table = pd.read_csv(summary_file, sep='\t')

    return asv_tax, filtered_asv_table, true_comp, summary_table

def prepare_taxa_sets(asv_tax, asv_table, true_comp, rank, sample_name):
    # Merge ASV table and tax info
    merged = pd.merge(asv_table[['ASV_ID', sample_name]], asv_tax, on='ASV_ID')
    merged = merged.rename(columns={sample_name: 'N_reads'})
    merged = merged[['ASV_ID', 'N_reads', rank]]

    # Drop NA or empty rank entries
    merged = merged[merged[rank].notna() & (merged[rank] != '')]

    # Group by rank and sum reads
    taxon_reads = merged.groupby(rank)[f'N_reads'].sum().reset_index()
    taxon_reads = taxon_reads.sort_values(by='N_reads', ascending=False)

    # Get true taxa at specified rank
    true_taxa = set(true_comp[rank].dropna().unique())

    return taxon_reads, true_taxa

def f1_score(y_correct, y_false, y_undetected):
    tp = y_correct
    fp = y_false
    fn = y_undetected
    
    if tp == 0 and (fp > 0 or fn > 0):
        return 0.0  # no true positives, but some errors
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    
    if precision + recall == 0:
        return 0.0
    return 2 * (precision * recall) / (precision + recall)


def evaluate_detection(taxon_reads, true_taxa):
    """
    Compute final detection statistics and F1 score.

    Parameters
    ----------
    taxon_reads : pd.DataFrame
        DataFrame with detected taxa at the chosen rank.
        Must have the taxon name in the first column.
    true_taxa : set
        Set of true taxa at the chosen rank.

    Returns
    -------
    tp : int Number of correctly detected taxa. 
    fp : int Number of falsely detected taxa.
    fn : int Number of undetected true taxa.
    f1 : float F1 score for the final detection state.
    """

    # All detected taxa
    detected_taxa = set(taxon_reads.iloc[:, 0]) if not taxon_reads.empty else set()

    tp = len(detected_taxa & true_taxa)
    fp = len(detected_taxa - true_taxa)
    fn = len(true_taxa - detected_taxa)

    f1 = f1_score(tp, fp, fn)

    return tp, fp, fn, f1

def main():
    args = parse_args()
    asv_tax, asv_table, true_comp, summary_table = read_files(
        args.dataset_name,
        args.run,
        args.true,
        args.X,
        args.out_suffix
    )

    sample_names = asv_table.columns[1:]

    for sample in sample_names:
        taxon_reads, true_taxa = prepare_taxa_sets(
            asv_tax,
            asv_table,
            true_comp,
            args.rank,
            sample
        )

        TP, FP, FN, f1 = evaluate_detection(taxon_reads, true_taxa)

        output = {
            "run": args.run,
            "sample": sample,
            "f1": f1,
            "TP": TP,
            "FP": FP,
            "FN": FN
        }

        print(json.dumps(output))

if __name__ == "__main__":
    main()

