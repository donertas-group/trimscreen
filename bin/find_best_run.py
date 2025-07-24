#!/usr/bin/env python
import pandas as pd
import sys
import math
import argparse
import json 
from itertools import combinations

def get_step(df):
    # Extract the first three digits after 'run_' and convert to int
    lenf = df['run'].str.extract(r'run_(\d{3})')[0].astype(int).unique().tolist()

    # Generate all unique pairs and compute absolute differences
    diffs = [abs(a - b) for a, b in combinations(lenf, 2)]

    # Return the minimum difference
    return min(diffs) if diffs else None


def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Find the best run")
    parser.add_argument("-i", "--input", required=True, help="Input file")
    parser.add_argument("-m", "--metadata", required=False, help="Optional metadata tsv table with same format as required by nf-core/ampliseq")
    parser.add_argument("-t", "--taxlevels", required=True, nargs="+", help="Taxonomic levels on which runs are evaluated")

    return parser.parse_args()


def main():
    args = parse_args()

    filtered_table_csv = args.input
    metadata_csv = args.metadata
    taxlevels = args.taxlevels

    df0 = pd.read_csv(filtered_table_csv)
    if metadata_csv is not None:
        metadata = pd.read_csv(metadata_csv)
        sample_ids = metadata[metadata['condition'] == 'sample']['ID'].values
        df = df0[df0['sample'].isin(sample_ids)]
    else:
        df = df0
    
    # extract screening step
    step = get_step(df)
    N = math.ceil(5 / step) # min nruns for highest richness to be considered 

    # Initialize dictionaries to track best runs for Genus and Phylum
    best_runs = {}

    # Iterate over each sample to find its best runs based on two chosen ranks, e.g. Genus and Family
    for sample in df['sample'].unique():
        sample_data = df[df['sample'] == sample]

        summary_table = (
            sample_data
            .groupby([taxlevels[1], taxlevels[0]])
            .agg(nruns=('run', 'count'))
            .sort_values(by=[taxlevels[1], taxlevels[0]], ascending=False)
            .reset_index()
            .query('nruns >= @N')
        )

        # Get the top Family + Genus values
        if not summary_table.empty:
            top_family = summary_table.iloc[0][taxlevels[0]]
            top_genus = summary_table.iloc[0][taxlevels[1]]
            best_rows = sample_data[(sample_data[taxlevels[1]] == top_genus) & (sample_data[taxlevels[0]] == top_family)]
        else:
            print(f"Best runs not found for sample {sample}.")

    for _, row in best_rows.iterrows():
        best_run_id = row['run']
        best_runs[best_run_id] = best_runs.get(best_run_id, 0) + 1

    # Convert the dictionaries into DataFrames for easier ranking
    summ = pd.DataFrame(list(best_runs.items()), columns=['run', 'counts']).sort_values(by='counts', ascending=False)

    # output json formatted list as stdout
    print(json.dumps([summ.iloc[0]['run']]))

    # Open the file for writing the results
    with open('report.txt', 'w') as f:

        f.write(f"Runs giving highest richness at {taxlevels[0]} and {taxlevels[1]} levels:\n")
        for _, row in summ.head(15).iterrows():
            f.write(f"{row['run']}, best for {row['counts']} samples\n")


if __name__ == "__main__":
    sys.exit(main())
