#!/usr/bin/env python3

# The script expects the input files in a structured directory, as described.
# It uses read count as the x-axis and plots detection stats at the user-specified taxonomic rank.
# Make sure column names match exactly, especially for ranks (Genus, Family, etc.) — they are case-sensitive.
# example usage:
# ./compare_w_mock.py -M 13 -r run_180190 -R Genus --true_composition_file true_composition.csv

import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description="Compare 16S data with true mock composition.")
    parser.add_argument("-M", "--mock", required=True, type=int, help="Mock community number")
    parser.add_argument("-r", "--run", required=True, type=str, help="Run ID")
    parser.add_argument("-R", "--rank", required=True, type=str, help="Taxonomic rank to compare (e.g., Genus)")
    parser.add_argument("--true_composition_file", required=True, type=str, help="Path to true composition CSV file")
    return parser.parse_args()

def read_files(mock, run_id, true_comp_file):
    result_dir = "/scratch/shire/ssd/pipeline/16s_nf_pipeline"  # Modify if needed
    data_dir = "/scratch/shire/data/nj/raw_data/published/mockrobiota"  # Modify if needed

    asv_tax_file = os.path.join(result_dir, f"mock{mock}", "output/runs", run_id, "dada2/ASV_tax_species.silva_138.tsv")
    asv_table_file = os.path.join(result_dir, f"mock{mock}", "output/runs", run_id, "dada2/ASV_table.tsv")


    true_comp_path = os.path.join(data_dir, f"mock-{mock}", true_comp_file)

    asv_tax = pd.read_csv(asv_tax_file, sep='\t')
    asv_table = pd.read_csv(asv_table_file, sep='\t')
    # Rename the second column (which is mock-specific) to 'N_reads'
    read_col_name = asv_table.columns[1]
    asv_table = asv_table.rename(columns={read_col_name: 'N_reads'})

    true_comp = pd.read_csv(true_comp_path)

    return asv_tax, asv_table, true_comp

def prepare_taxa_sets(asv_tax, asv_table, true_comp, rank):
    # Merge ASV table and tax info
    merged = pd.merge(asv_table, asv_tax, on='ASV_ID')
    merged = merged[['ASV_ID', 'N_reads', rank]]

    # Drop NA or empty rank entries
    merged = merged[merged[rank].notna() & (merged[rank] != '')]

    # Group by rank and sum reads
    taxon_reads = merged.groupby(rank)[f'N_reads'].sum().reset_index()
    taxon_reads = taxon_reads.sort_values(by='N_reads', ascending=False)

    # Get true taxa at specified rank
    true_taxa = set(true_comp[rank].dropna().unique())

    return taxon_reads, true_taxa

def evaluate_detection(taxon_reads, true_taxa):
    x_reads = []
    y_correct = []
    y_false = []
    y_undetected = []

    seen = set()

    for _, row in taxon_reads.iterrows():
        taxon = row.iloc[0]
        reads = row['N_reads']
        seen.add(taxon)

        x_reads.append(reads)
        correct = len(seen & true_taxa)
        false = len(seen - true_taxa)
        undetected = len(true_taxa - seen)

        y_correct.append(correct)
        y_false.append(false)
        y_undetected.append(undetected)

    return x_reads, y_correct, y_false, y_undetected

def plot_results(x, correct, false, undetected, mock, run, rank):
    import os

    output_dir = "./output"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"detection_curve_mock{mock}.{run}.{rank}.png")

    plt.figure(figsize=(10, 6))
    plt.plot(x, correct, label=f"Correctly detected {rank}", color='green')
    plt.plot(x, false, label=f"Falsely detected {rank}", color='red')
    plt.plot(x, undetected, label=f"Undetected true {rank}", color='blue')
    plt.xscale('log')
    plt.xlabel("Read number (log scale)")
    plt.ylabel(f"Number of taxa at rank '{rank}'")
    plt.title(f"Taxon detection at {rank} level")
    plt.legend()
    plt.tight_layout()
    plt.grid(True)
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Plot saved to {output_file}")


def main():
    args = parse_args()
    asv_tax, asv_table, true_comp = read_files(args.mock, args.run, args.true_composition_file)
    taxon_reads, true_taxa = prepare_taxa_sets(asv_tax, asv_table, true_comp, args.rank)
    x, correct, false, undetected = evaluate_detection(taxon_reads, true_taxa)
    plot_results(x, correct, false, undetected, args.mock, args.run, args.rank)

if __name__ == "__main__":
    main()

