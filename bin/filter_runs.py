#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
import argparse

def parse_args(args=None):

    parser = argparse.ArgumentParser(description="Filter runs and output table.")
    parser.add_argument("-i", "--input", required=True, help="Input table")
    parser.add_argument("-c", "--columns", required=True, nargs="+", help="Columns to filter")

    return parser.parse_args()


def filter_runs(df, columns):
    # Initialize an empty dictionary to store the results
    results = {}

    # Loop through each column provided
    for column in columns:

        # covert to numeric
        df[columns] = df[columns].apply(pd.to_numeric, errors='coerce')
        # Group by 'run' and calculate median and standard deviation
        grouped = df.groupby('run')[column].agg(['median', 'std']).reset_index()

        # Compute global 2-sigma range for median
        median_mean = grouped['median'].mean()
        median_std = grouped['median'].std()
        median_lower, median_upper = median_mean - 2 * median_std, median_mean + 2 * median_std

        # Compute global 2-sigma range for standard deviation
        std_mean = grouped['std'].mean()
        std_std = grouped['std'].std()
        std_lower, std_upper = std_mean - 2 * std_std, std_mean + 2 * std_std

        # Filter runs whose median and std fall within their respective 3-sigma range
        filtered_grouped = grouped[
            (grouped['median'] >= median_lower) & (grouped['median'] <= median_upper) &
            (grouped['std'] >= std_lower) & (grouped['std'] <= std_upper)
        ]

        # Store the filtered runs in the dictionary
        results[column] = set(filtered_grouped['run'].values)

    # Find the intersection of all top runs across columns
    common_runs = set.intersection(*results.values())

    # Filter the dataframe to keep only the rows corresponding to the common 'run's
    filtered_df = df[df['run'].isin(common_runs)]

    return filtered_df



def main():
    args = parse_args()

    file = args.input
    columns_to_filter = args.columns

    full_table = pd.read_csv(file)

    # filter runs by evaluating the median and sd of choosen columns
    filtered_table = filter_runs(full_table, columns_to_filter)

    filtered_table.to_csv("filtered_table.csv", index=False)


if __name__ == "__main__":
    sys.exit(main())





