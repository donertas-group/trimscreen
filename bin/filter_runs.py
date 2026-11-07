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
    parser.add_argument("-o", "--output", required=True, help="Output table")
    parser.add_argument("-r", "--lowest_relative_retention", type=float, default=0.5,
                        help="Minimum retained_reads_percent of a sample, as a fraction of "
                             "that sample's maximum across runs (default: 0.5)")
    parser.add_argument("-n", "--min_reads", type=int, default=2000,
                        help="Samples must have more than this many reads (nreads > min_reads) "
                             "(default: 2000)")
    return parser.parse_args()

def find_good_runs(
    df: pd.DataFrame,
    min_reads: int,
    lowest_relative_retention: float,
) -> List[str]:
    """
    Return the runs in which every sample passes both filters below.

    Filters (applied in this order):
      1) Read depth: rows with nreads <= min_reads (or non-numeric/missing
         nreads) are discarded.
      2) Relative retention: among the rows that survived step 1, the maximum
         retained_reads_percent is computed per sample across all runs.
         A row is kept only if its retained_reads_percent is >=
         lowest_relative_retention * that maximum. Because step 1 comes first,
         runs that fail the read-depth filter do not set the reference
         maximum. Missing/non-numeric retention values fail this filter.

    Run-level rule:
      A run is returned only if the number of unique samples it has after
      filtering equals the number it had originally, i.e. every sample in the
      run keeps at least one row that passes both filters. Counting is by
      unique sample name, so if a sample has several rows in a run, one
      passing row is enough for that sample to count as retained.

    Args:
      df: DataFrame with at least the columns
          ['run', 'sample', 'nreads', 'retained_reads_percent'].
      min_reads: exclusive lower bound on nreads (rows need nreads > min_reads).
      lowest_relative_retention: fraction (0-1) of the per-sample maximum
          retained_reads_percent that a row must reach to be kept.

    Returns:
      List of run names that pass all criteria.
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
    full_table = pd.read_csv(file)

    # find runs in which every sample has enough reads and sufficient
    # retention relative to the best run for that sample

    good_runs = find_good_runs(
        full_table,
        min_reads=int(args.min_reads),
        lowest_relative_retention=float(args.lowest_relative_retention),
    )

    filtered_table = full_table[full_table['run'].isin(good_runs)]

    # rarefaction depth: the minimum nreads across all rows of all good runs
    # (a single depth shared by every run)
    D = filtered_table['nreads'].min()

    # output as [runID, depth] as stdout for nf process
    out = [[run, int(D)] for run in good_runs]
    print(json.dumps(out))
    
    filtered_table.to_csv(args.output, index=False)


if __name__ == "__main__":
    sys.exit(main())

