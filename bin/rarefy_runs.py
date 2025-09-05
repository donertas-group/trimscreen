#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
import argparse

def parse_args(args=None):

    parser = argparse.ArgumentParser(description="Rarefy abundance table of filtered runs and output table.")
    parser.add_argument("-i", "--input", required=True, help="ASV table to be rarefied")
    parser.add_argument("-d", "--depth", required=True, help="Rarefaction depth")

    return parser.parse_args()

def rarefy_runs(abundance_table, rarefy_depth):

    return abundance_table


def main():
    args = parse_args()

    file = args.input
    depth = args.depth

    asv_table = pd.read_csv(file, sep="\t")

    # rarefy runs 
    rarefied_table = rarefy_runs(asv_table, depth)

    rarefied_table.to_csv("ASV_table_rarefied.tsv.gz", sep="\t", index=False, compression='gzip')


if __name__ == "__main__":
    sys.exit(main())

