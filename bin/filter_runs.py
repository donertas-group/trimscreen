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
    parser.add_argument("-c", "--columns", required=True, nargs="+", help="Columns to filter")

    return parser.parse_args()

def find_good_runs(df: pd.DataFrame, columns: List[str]) -> List[str]:
    """
    Return run(s) that keep the largest total number of samples after applying,
    for each column, a sample-level mean±std filter across runs.

    Args:
      df: DataFrame containing at least ['run','sample'] + columns
      columns: list of column names to apply the procedure to (e.g. ['Genus','Family'])
      min_run_samples: only consider runs that have >= this many samples (default 2)

    Returns:
      list of run names (ties are returned as multiple runs)
    """
    # defensive copy
    df = df.copy()

    # ensure numeric columns
    df[columns] = df[columns].apply(pd.to_numeric, errors='coerce')
    runs = df['run'].unique().tolist()
    # consider only runs with >= min_run_samples samples (this defines df_multi)
    #run_counts = df['run'].value_counts()
    #multi_sample_runs = run_counts[run_counts > 1].index.tolist()
    #df_multi = df[df['run'].isin(multi_sample_runs)].copy()

    #single_sample_runs = run_counts[run_counts == 1].index.tolist()
    #df_single = df[df['run'].isin(single_sample_runs)].copy()

    #if df_multi.empty:
    #    return []

    # initialize total retained-sample counts per run
    total_retained_per_run: Dict[str,int] = {run: 0 for run in runs}

    for col in columns:
        # compute per-sample mean and std across runs
        sample_stats = (
            df
            .groupby('sample')[col]
            .agg(mean='mean', std='std')
            .reset_index()
        )

        # treat NaN std (single observation for that sample) as 0
        sample_stats['std'] = sample_stats['std'].fillna(0.0)

        # compute lower/upper bounds (mean ± std)
        sample_stats['low'] = sample_stats['mean'] - sample_stats['std']
        sample_stats['high'] = sample_stats['mean'] + sample_stats['std']

        # join bounds back to df
        merged = df.merge(
            sample_stats[['sample','mean','std','low','high']],
            on='sample',
            how='left',
            validate='m:1'  # many rows in df to one row in sample_stats
        )

        # keep run-rows whose value is within sample_mean ± sample_std
        # (NaN values in the column will be treated as not within bounds)
        within_mask = (
            merged[col].notna() &
            (merged[col] >= merged['low']) &
            (merged[col] <= merged['high'])
        )
        kept = merged.loc[within_mask, ['run','sample']].drop_duplicates()

        # count how many distinct samples each run kept for THIS column
        if not kept.empty:
            per_run_counts = kept.groupby('run')['sample'].nunique()
            # add to totals (runs with 0 kept samples naturally add 0)
            for run, cnt in per_run_counts.items():
                total_retained_per_run[run] += int(cnt)
        # otherwise no run gains any counts for this column

    # find run(s) with the maximum total retained samples
    max_kept = max(total_retained_per_run.values())
    good_runs = [run for run, cnt in total_retained_per_run.items() if cnt == max_kept]

    return good_runs


def main():
    args = parse_args()

    file = args.input
    columns_to_filter = args.columns

    full_table = pd.read_csv(file)

    # find good runs by evaluating the per-sample median and sd of choosen columns
    good_runs = find_good_runs(full_table, columns_to_filter)
    filtered_table = full_table[full_table['run'].isin(good_runs)]

    # get rarefaction depth
    D = filtered_table['nreads'].min()

    # output as [runID, depth] as stdout for nf process
    out = [[run, int(D)] for run in good_runs]
    print(json.dumps(out))
    
    filtered_table.to_csv("filtered_table.csv", index=False)


if __name__ == "__main__":
    sys.exit(main())

