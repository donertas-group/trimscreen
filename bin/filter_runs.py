#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
import argparse
import json

def parse_args(args=None):

    parser = argparse.ArgumentParser(description="Filter runs and output table.")
    parser.add_argument("-i", "--input", required=True, help="Input table")
    parser.add_argument("-c", "--columns", required=True, nargs="+", help="Columns to filter")

    return parser.parse_args()

def filter_runs(df, columns):
    import numpy as np

    # Convert specified columns to numeric
    df[columns] = df[columns].apply(pd.to_numeric, errors='coerce')

    # Track sets of "good" runs for each column
    results = {}

    for column in columns:
        # Count how many samples each run has
        run_counts = df['run'].value_counts()

        # Split into two DataFrames:
        # - one with multiple samples per run
        # - one with single sample per run
        multi_sample_runs = run_counts[run_counts > 1].index
        single_sample_runs = run_counts[run_counts == 1].index

        df_multi = df[df['run'].isin(multi_sample_runs)]
        df_single = df[df['run'].isin(single_sample_runs)]

        # -------------------------
        # Case 1: Multi-sample runs
        # -------------------------
        if not df_multi.empty:
            grouped = df_multi.groupby('run')[column].agg(['median', 'std']).reset_index()

            # Compute global stats for medians and stds
            median_mean = grouped['median'].mean()
            median_std = grouped['median'].std()
            median_lower, median_upper = median_mean - 1 * median_std, median_mean + 1 * median_std

            std_mean = grouped['std'].mean()
            std_std = grouped['std'].std()
            std_lower, std_upper = std_mean - 1 * std_std, std_mean + 1 * std_std

            # Keep runs within 1-sigma bounds
            filtered_multi = grouped[
                (grouped['median'] >= median_lower) & (grouped['median'] <= median_upper) &
                (grouped['std'] >= std_lower) & (grouped['std'] <= std_upper)
            ]

            good_multi_runs = set(filtered_multi['run'])
        else:
            good_multi_runs = set()

        # -------------------------
        # Case 2: Single-sample runs
        # -------------------------
        if not df_single.empty:
            values = df_single[column]
            val_mean = values.mean()
            val_std = values.std()
            lower, upper = val_mean - 1 * val_std, val_mean + 1 * val_std

            filtered_single = df_single[
                (df_single[column] >= lower) & (df_single[column] <= upper)
            ]

            good_single_runs = set(filtered_single['run'])
        else:
            good_single_runs = set()

        # Merge good runs for this column
        results[column] = good_multi_runs.union(good_single_runs)

    # Find runs that passed *all* filters
    common_runs = set.intersection(*results.values())

    # Return filtered DataFrame
    filtered_df = df[df['run'].isin(common_runs)]

    return filtered_df


def main():
    args = parse_args()

    file = args.input
    columns_to_filter = args.columns

    full_table = pd.read_csv(file)

    # filter runs by evaluating the median and sd of choosen columns
    filtered_table = filter_runs(full_table, columns_to_filter)


    # get rarefaction depth
    D = filtered_table['nreads'].min()

    # output as [runID, depth] as stdout for nf process
    out = [[run, int(D)] for run in filtered_table['run'].tolist()]
    print(json.dumps(out))
    
    filtered_table.to_csv("filtered_table.csv", index=False)


if __name__ == "__main__":
    sys.exit(main())

