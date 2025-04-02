#!/usr/bin/env python
import pandas as pd
import sys
import math
import argparse

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Find the best run")
    parser.add_argument("-i", "--input", required=True, help="Input file")
    parser.add_argument("-m", "--metadata", required=False, help="Optional metadata tsv table with same format as required by nf-core/ampliseq")
    parser.add_argument("-t", "--taxlevels", required=True, nargs="+", help="Taxonomic levels on which runs are evaluated")

    return parser.parse_args()
 

def main():
    args = parse_args()

    filtered_table_csv = args.input
    metadata_tsv = args.metadata
    taxlevels = args.taxlevels

    # Read the CSV file and metadata table into pandas DataFrames
    df0 = pd.read_csv(filtered_table_csv)
    metadata = pd.read_csv(metadata_tsv, sep='\t')  # Assuming tab-separated values

    # Filter out the samples that are not marked as 'sample' in the metadata table
    sample_ids = metadata[metadata['condition'] == 'sample']['ID'].values

    # Filter df to keep only the rows with the matching sample IDs
    df = df0[df0['sample'].isin(sample_ids)]

    # Initialize dictionaries to track best runs for Genus and Phylum
    best1_runs = {}
    best2_runs = {}

    # Iterate over each sample to find the best runs based on Genus and Phylum
    for sample in df['sample'].unique():
        sample_data = df[df['sample'] == sample]

        # Find the best runs for Genus (all runs with the highest Genus value)
        max2_value = sample_data[taxlevels[1]].max()
        best2_rows = sample_data[sample_data[taxlevels[1]] == max2_value]
        for _, row in best2_rows.iterrows():
            best2_run_id = row['run']
            best2_runs[best2_run_id] = best2_runs.get(best2_run_id, 0) + 1

        # Find the best runs for Phylum (all runs with the highest Phylum value)
        max1_value = sample_data[taxlevels[0]].max()
        best1_rows = sample_data[sample_data[taxlevels[0]] == max1_value]
        for _, row in best1_rows.iterrows():
            best1_run_id = row['run']
            best1_runs[best1_run_id] = best1_runs.get(best1_run_id, 0) + 1

    # Convert the dictionaries into DataFrames for easier ranking
    df2 = pd.DataFrame(list(best2_runs.items()), columns=['run', 'value'])
    df1 = pd.DataFrame(list(best1_runs.items()), columns=['run', 'value'])

    # Rank the runs for Genus (with ties having the same rank)
    df2['rank'] = df2['value'].rank(method='min', ascending=False).astype(int)  # 'min' ties the rank

    # Rank the runs for Phylum (with ties having the same rank)
    df1['rank'] = df1['value'].rank(method='min', ascending=False).astype(int)  # 'min' ties the rank

    # Convert back to dictionaries
    ranking2 = dict(zip(df2['run'], df2['rank']))
    ranking1 = dict(zip(df1['run'], df1['rank']))

    # Combine ranks
    combined_ranks = {}
    for run in set(ranking2.keys()).union(ranking1.keys()):
        taxlevel2_rank = ranking2.get(run, len(ranking2) + 1)  # Assign lowest rank for missing runs
        taxlevel1_rank = ranking1.get(run, len(ranking1) + 1)  # Assign lowest rank for missing runs
        combined_ranks[run] = math.sqrt(taxlevel2_rank**2 + taxlevel1_rank**2)

    # Find the minimum total rank value
    min_rank = min(combined_ranks.values())

    # Identify all runs with the minimum rank
    bests = [(run, rank) for run, rank in combined_ranks.items() if rank == min_rank]
    best_runs = [run for run, rank in bests]
    print(best_runs)

    # Open the file for writing
    with open('report.txt', 'w') as f:

        # Output the results for the first taxonomic level
        f.write(f"Runs with highest {taxlevels[0]}:\n")
        for run, rank in sorted(ranking2.items(), key=lambda item: item[1])[:15]:
            f.write(f"{run}, rank: {rank}, best for {best2_runs.get(run, 0)} samples\n")

        # Output the results for the second taxonomic level
        f.write(f"\nRuns with highest {taxlevels[1]}:\n")
        for run, rank in sorted(ranking1.items(), key=lambda item: item[1])[:15]:
            f.write(f"{run}, rank: {rank}, best for {best1_runs.get(run, 0)} samples\n")

        # Print all the best runs
        f.write("\nBest Runs:\n")
        for run, rank in bests:
            f.write(f"Run: {run}, Combined Rank: {rank}\n")


if __name__ == "__main__":
    sys.exit(main())
