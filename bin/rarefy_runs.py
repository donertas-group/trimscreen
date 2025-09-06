#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import argparse

def parse_args(args=None):

    parser = argparse.ArgumentParser(description="Rarefy abundance table of filtered runs and output table.")
    parser.add_argument("-i", "--input", required=True, help="ASV table to be rarefied")
    parser.add_argument("-d", "--depth", type=int, required=True, help="Number of reads per sample to rarefy to.")

    return parser.parse_args()


def rarefy_runs(asv_table: pd.DataFrame, rarefy_depth: int, random_state=None) -> pd.DataFrame:
    """
    Rarefy an ASV abundance table to a given depth (reads per sample).
    
    Parameters
    ----------
    asv_table : pd.DataFrame
        ASV abundance table with ASVs as rows and samples as columns.
    rarefy_depth : int
        Target depth (number of reads per sample).
    random_state : int, optional
        Seed for reproducibility.
    
    Returns
    -------
    pd.DataFrame
        Rarefied abundance table (same shape as input).
        Samples with fewer reads than rarefy_depth are filled with NaNs.
    """
    rng = np.random.default_rng(random_state)
    rarefied = pd.DataFrame(0, index=asv_table.index, columns=asv_table.columns)
    
    for sample in asv_table.columns:
        counts = asv_table[sample].values
        total_reads = counts.sum()
        
        if total_reads < rarefy_depth:
            # Cannot rarefy this sample, not enough reads
            rarefied[sample] = np.nan
            continue
        
        # Expand into vector of read-level ASV identities
        asv_ids = np.repeat(np.arange(len(counts)), counts)
        
        # Subsample without replacement
        chosen = rng.choice(asv_ids, size=rarefy_depth, replace=False)
        
        # Count back into ASV space
        subsampled_counts = np.bincount(chosen, minlength=len(counts))
        
        rarefied[sample] = subsampled_counts
    
    return rarefied


def main():
    args = parse_args()

    file = args.input
    depth = args.depth

    asv_table = pd.read_csv(file, sep="\t", index_col=0)
    asv_table = asv_table.fillna(0).astype(int)
    # rarefy runs 
    rarefied_table = rarefy_runs(asv_table, depth, random_state=42)

    rarefied_table.to_csv("ASV_table_rarefied.tsv.gz", sep="\t", index=True, compression='gzip')


if __name__ == "__main__":
    sys.exit(main())

