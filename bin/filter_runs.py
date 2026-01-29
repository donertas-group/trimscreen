#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
import argparse
import json
from typing import List, Tuple, Dict

def parse_args(args=None):

    parser = argparse.ArgumentParser(description="Filter runs and output table.")
    parser.add_argument("-i", "--input", required=True, help="Input table")
    parser.add_argument("-r", "--lowest_relative_retention", type=float, default=0.8)
    parser.add_argument("-n", "--min_reads", type=int, default=2000)

    return parser.parse_args()

def find_good_runs(
    df: pd.DataFrame,
    min_reads: int = 2000,
    lowest_relative_retention: float = 0.8,
) -> List[str]:
    """
    Return runs where:
      1) Every row has nreads > min_reads
      2) For each sample, retained_reads_percent is at least
         lowest_relative_retention * max(retained_reads_percent) across runs
      3) No samples are filtered out within a run (strict run-level filter)

    Args:
      df: DataFrame with at least ['run', 'sample', 'nreads', 'retained_reads_percent']
      min_reads: minimum reads per sample
      lowest_relative_retention: fraction of max retained_reads_percent per sample

    Returns:
      List of run names that pass all criteria
    """

    # Ensure numeric columns
    df = df.copy()
    df['nreads'] = pd.to_numeric(df['nreads'], errors='coerce')
    df['retained_reads_percent'] = pd.to_numeric(
        df['retained_reads_percent'], errors='coerce'
    )

    # -----------------------------
    # Step 1: filter by min_reads
    # -----------------------------
    filtered = df[df['nreads'] > min_reads]

    # -----------------------------
    # Step 2: relative retention per sample
    # -----------------------------
    max_retention_per_sample = (
        filtered
        .groupby('sample')['retained_reads_percent']
        .max()
        .rename('max_retention')
    )

    filtered = filtered.merge(
        max_retention_per_sample,
        on='sample',
        how='left',
        validate='m:1'
    )

    filtered = filtered[
        filtered['retained_reads_percent']
        >= lowest_relative_retention * filtered['max_retention']
    ]

    # -----------------------------
    # Step 3: strict run-level filter
    # -----------------------------
    # Original sample counts per run
    original_counts = (
        df
        .groupby('run')['sample']
        .nunique()
    )

    # Remaining sample counts per run after filtering
    remaining_counts = (
        filtered
        .groupby('run')['sample']
        .nunique()
    )

    # Keep runs where no samples were lost
    good_runs = [
        run
        for run, n_samples in original_counts.items()
        if remaining_counts.get(run, 0) == n_samples
    ]

    return good_runs

def main():
    args = parse_args()

    file = args.input
    columns_to_filter = args.columns

    full_table = pd.read_csv(file)

    # find good runs by evaluating the per-sample median and sd of choosen columns
    good_runs = find_good_runs(
        full_table,
        min_reads=int(args.min_reads),
        lowest_relative_retention=float(args.lowest_relative_retention),
    )

    filtered_table = full_table[full_table['run'].isin(good_runs)]

    # get rarefaction depth
    D = filtered_table['nreads'].min()

    # output as [runID, depth] as stdout for nf process
    out = [[run, int(D)] for run in good_runs]
    print(json.dumps(out))
    
    filtered_table.to_csv("filtered_table.csv", index=False)


if __name__ == "__main__":
    sys.exit(main())

