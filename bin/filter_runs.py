#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
import argparse
import json
from typing import List

# Column that identifies an individual biological sample across runs.
# NOTE: adjust this if your samplerun summary table uses a different
# column name for the sample identifier (e.g. "sample_id").
SAMPLE_COL = "sample"


def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Filter runs and output table.")
    parser.add_argument("-i", "--input", required=True, help="Input table")
    parser.add_argument("-o", "--output", required=True, help="Output table")
    parser.add_argument("-r", "--lowest_relative_retention", type=float, default=0.5)
    parser.add_argument("-n", "--min_reads", type=int, default=2000)
    parser.add_argument(
        "-s", "--expected_samples", default=None,
        help="Comma-separated list of all expected sample names (e.g. "
             "derived from the samplesheet). If provided, samples that "
             "never appear anywhere in the input table are reported as a "
             "warning (separate from the min_reads-driven errors)."
    )

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
      3) The run contains the FULL set of samples (i.e. no sample is
         missing from it), so that runs are actually comparable to one
         another. A run that is only missing part of the sample set would
         otherwise pass criteria 1 and 2 undetected.
    """
    if SAMPLE_COL not in df.columns:
        raise ValueError(
            f"Column '{SAMPLE_COL}' not found in input table; cannot check "
            "whether each run contains the full set of samples. Set "
            "SAMPLE_COL to the correct column name."
        )

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
    # Step 3: Find runs that don't contain the full set of samples
    # -------------------------------------------------------------
    # "Full set" is defined as every distinct sample seen anywhere in the
    # table. Note: if a sample is missing from ALL runs (e.g. it dropped
    # out of the table entirely), that can't be detected from this table
    # alone - only run-vs-run comparisons of samples present in the table.
    all_samples = set(df[SAMPLE_COL].dropna().unique())
    run_sample_sets = df.groupby('run')[SAMPLE_COL].apply(set)
    runs_missing_samples = [
        run for run in run_sample_sets.index
        if run_sample_sets[run] != all_samples
    ]

    # -------------------------------------------------------------
    # Step 4: Combine criteria
    # -------------------------------------------------------------
    # Keep runs that pass retention, do not contain any low-read samples,
    # and contain the full sample set
    good_runs = [
        run for run in runs_passing_retention
        if run not in runs_failing_reads
        and run not in runs_missing_samples
    ]

    return good_runs


def find_chronic_low_read_samples(df: pd.DataFrame, min_reads: int) -> List[str]:
    """
    Identify samples that have nreads <= min_reads in EVERY run in which
    they were sequenced.

    A sample that is low in one run but fine in another is just a
    run-specific issue. A sample that is *always* low, no matter which
    run it was sequenced in, can never let its run(s) pass the min_reads
    filter. If such a sample happens to be present in all runs, it will
    silently take every run down with it, leaving 0 good runs with no
    obvious explanation downstream. This function surfaces those samples
    so a specific, actionable error can be raised.
    """
    if SAMPLE_COL not in df.columns:
        # Can't attribute failures to a specific sample without a sample
        # identifier column; return no diagnosis rather than guessing.
        return []

    df = df.copy()
    df['nreads'] = pd.to_numeric(df['nreads'], errors='coerce')

    # True if the sample clears min_reads in at least one of its runs
    passes_somewhere = df.groupby(SAMPLE_COL)['nreads'].apply(
        lambda reads: (reads > min_reads).any()
    )

    chronic_low_samples = passes_somewhere[~passes_somewhere].index.tolist()

    return chronic_low_samples


def find_missing_from_table_samples(
    df: pd.DataFrame,
    expected_samples_str: str,
) -> List[str]:
    """
    Compare the full set of expected sample names (e.g. from the
    samplesheet, passed in as a comma-separated string) against the
    samples that actually appear anywhere in the input table (i.e. in at
    least one run, regardless of min_reads/retention filtering).

    Samples missing here were already dropped before reaching this module
    entirely - e.g. by upstream ASV length filtering - which is a
    different problem than a sample failing the min_reads/retention
    criteria applied by this script. This is reported as a warning only;
    it does not affect which runs are selected.
    """
    if not expected_samples_str:
        return []

    expected_samples = {
        s.strip() for s in expected_samples_str.split(',') if s.strip()
    }
    observed_samples = (
        set(df[SAMPLE_COL].dropna().unique()) if SAMPLE_COL in df.columns else set()
    )

    return sorted(expected_samples - observed_samples)


def main():
    args = parse_args()

    full_table = pd.read_csv(args.input)

    good_runs = find_good_runs(
        full_table,
        min_reads=int(args.min_reads),
        lowest_relative_retention=float(args.lowest_relative_retention),
    )

    chronic_low_samples = find_chronic_low_read_samples(
        full_table,
        min_reads=int(args.min_reads),
    )

    missing_from_table_samples = find_missing_from_table_samples(
        full_table,
        args.expected_samples,
    )

    filtered_table = full_table[full_table['run'].isin(good_runs)]

    if filtered_table.empty:
        # Line 1: good runs (empty here)
        # Line 2: chronic low-read samples (error-level, see caller)
        # Line 3: samples missing from the table entirely (warning-level)
        print(json.dumps([]))
        print(json.dumps(chronic_low_samples))
        print(json.dumps(missing_from_table_samples))
        # Optional: write an empty dataframe or handle gracefully
        filtered_table.to_csv(args.output, index=False)
        return 0

    # Get rarefaction depth based ONLY on the selected good runs
    D = filtered_table['nreads'].min()

    # Output as [runID, depth] as stdout for nextflow/snakemake process,
    # followed by chronic low-read samples (error-level) and samples
    # missing from the table entirely (warning-level), so the caller can
    # diagnose both kinds of failure/anomaly.
    out = [[run, int(D)] for run in good_runs]
    print(json.dumps(out))
    print(json.dumps(chronic_low_samples))
    print(json.dumps(missing_from_table_samples))

    filtered_table.to_csv(args.output, index=False)

if __name__ == "__main__":
    sys.exit(main())

