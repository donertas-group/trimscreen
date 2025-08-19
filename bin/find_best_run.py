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
    return min(diffs) if diffs else 5


def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Find the best run")
    parser.add_argument("-i", "--input", required=True, help="filtered_table_csv")
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
        sample_ids = metadata[metadata['condition'] == 'sample']['sampleID'].values
        df = df0[df0['sample'].isin(sample_ids)]
    else:
        df = df0

    step = get_step(df)
    N = math.ceil(5 / step)

    def get_best_runs_by_sort_order(df, primary_tax, secondary_tax, N):
        best_runs = {}
        for sample in df['sample'].unique():
            sample_data = df[df['sample'] == sample]
            summary_table = (
                sample_data
                .groupby([primary_tax, secondary_tax])
                .agg(nruns=('run', 'count'))
                .sort_values(by=[primary_tax, secondary_tax], ascending=False)
                .reset_index()
                .query('nruns >= @N')
            )
            if not summary_table.empty:
                top_primary = summary_table.iloc[0][primary_tax]
                top_secondary = summary_table.iloc[0][secondary_tax]
                best_rows = sample_data[
                    (sample_data[primary_tax] == top_primary) &
                    (sample_data[secondary_tax] == top_secondary)
                ]
                for _, row in best_rows.iterrows():
                    best_run_id = row['run']
                    best_runs[best_run_id] = best_runs.get(best_run_id, 0) + 1
            else:
                print(f"Best runs not found for sample {sample}.")
        return pd.DataFrame(list(best_runs.items()), columns=['run', 'counts']).sort_values(by='counts', ascending=False)

    # First: maximize taxlevels[1] (Genus before Family, for example)
    summ_1 = get_best_runs_by_sort_order(df, taxlevels[1], taxlevels[0], N)
    # Second: maximize taxlevels[0] (Family before Genus, for example)
    summ_2 = get_best_runs_by_sort_order(df, taxlevels[0], taxlevels[1], N)

    # Print top run from first strategy as JSON
    print(json.dumps([summ_1.iloc[0]['run'], summ_2.iloc[0]['run']]))

    with open('report.txt', 'w') as f:
        f.write(f"Runs giving highest richness by prioritizing {taxlevels[1]} > {taxlevels[0]}:\n")
        for _, row in summ_1.head(15).iterrows():
            f.write(f"{row['run']}, best for {row['counts']} samples\n")

        f.write("\n" + "-" * 60 + "\n\n")

        f.write(f"Runs giving highest richness by prioritizing {taxlevels[0]} > {taxlevels[1]}:\n")
        for _, row in summ_2.head(15).iterrows():
            f.write(f"{row['run']}, best for {row['counts']} samples\n")

if __name__ == "__main__":
    sys.exit(main())
