#!/usr/bin/env python

import argparse
import pandas as pd
import sys

def parse_args(args=None):

    parser = argparse.ArgumentParser(description="Merge all run files into one.")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input files")

    return parser.parse_args()


def main():
    args = parse_args()
    files = args.input

    dfs = []
    for i, file in enumerate(files):
        df = pd.read_csv(file)
        
        if i == 0:
            dfs.append(df)  # Keep headers in the first file
        else:
            dfs.append(df[df.columns])  # Ensure same column order, drop headers

    # Concatenate all DataFrames
    merged_df = pd.concat(dfs, ignore_index=True)

    # Save to a new CSV file
    merged_df.to_csv("full_table.csv", index=False)

if __name__ == "__main__":
    sys.exit(main())

