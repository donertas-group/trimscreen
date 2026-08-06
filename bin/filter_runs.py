#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
import argparse
import json
from typing import List

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Filter runs and output table.")
    parser.add_argument("-i", "--input", required=True, help="Input table")
    parser.add_argument("-o", "--output", required=True, help="Output table")
    parser.add_argument("-r", "--lowest_relative_retention", type=float, default=0.5)
    parser.add_argument("-n", "--min_reads", type=int, default=2000)

    return parser.parse_args()

def find_good_runs(
    df: pd.DataFrame,
    min_reads: int,
    lowest_relative_retention: float,
) -> List[str]:
    """
    Return runs where:
      1) Every row (sample) within the run has nreads > min_reads
      2) The median retained_reads_percent of the run is at least
         lowest_relative_retention * max(median retained_reads_percent across all runs)
    """
    # Ensure numeric columns
    df = df.copy()
    df['nreads'] = pd.to_numeric(df['nreads'], errors='coerce')
    df['retained_reads_percent'] = pd.to_numeric(
        df['retained_reads_percent'], errors='coerce'
    )

    # -------------------------------------------------------------
    # Step 1: Find runs where ANY sample fails the absolute min_reads
    # -------------------------------------------------------------
    # If a single sample in a run is below min_reads, invalidate the whole run
    runs_failing_reads = df[df['nreads'] <= min_reads]['run'].unique()
    
    # -------------------------------------------------------------
    # Step 2: Evaluate run-level median retention
    # -------------------------------------------------------------
    # Calculate the median retention *per run*
    run_medians = df.groupby('run')['retained_reads_percent'].median()
    
    # Find the maximum run median across the entire dataset
    max_run_median = run_medians.max()
    
    # Determine the passing threshold
    threshold = lowest_relative_retention * max_run_median
    
    # Get runs that meet or exceed this median threshold
    runs_passing_retention = run_medians[run_medians >= threshold].index.tolist()

    # -------------------------------------------------------------
    # Step 3: Combine criteria
    # -------------------------------------------------------------
    # Keep runs that pass retention AND do not contain any low-read samples
    good_runs = [
        run for run in runs_passing_retention 
        if run not in runs_failing_reads
    ]

    return good_runs

def main():
    args = parse_args()

    full_table = pd.read_csv(args.input)

    good_runs = find_good_runs(
        full_table,
        min_reads=int(args.min_reads),
        lowest_relative_retention=float(args.lowest_relative_retention),
    )

    filtered_table = full_table[full_table['run'].isin(good_runs)]

    if filtered_table.empty:
        print(json.dumps([]))
        # Optional: write an empty dataframe or handle gracefully
        filtered_table.to_csv(args.output, index=False)
        return 0

    # Get rarefaction depth based ONLY on the selected good runs
    D = filtered_table['nreads'].min()

    # Output as [runID, depth] as stdout for nextflow/snakemake process
    out = [[run, int(D)] for run in good_runs]
    print(json.dumps(out))
    
    filtered_table.to_csv(args.output, index=False)

if __name__ == "__main__":
    sys.exit(main())
