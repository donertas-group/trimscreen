#!/usr/bin/env python
import argparse

import subprocess
import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from skbio import diversity

from _process_run import process_run
from _process_run import filter_runs
from _process_run import add_group_statistics

#import config

def parse_args(args=None):

    parser = argparse.ArgumentParser(description="Compare runsi. Create table from multiple files.")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input files")

    return parser.parse_args()




def main():
    args = parse_args()
    
    files = args.input

classifier_dir = "dada2"

Ranks_to_analyse = ["Phylum","Family","Genus","Species"]











if __name__ == "__main__":
    sys.exit(main())

# Use the variables defined in config.py
#workdir = config.workdir
#run_dirs = config.run_dirs
#outdir = config.outdir

#params_file = config.params_file


all_runs = pd.concat([process_run(run_dir, classifier_dir, Ranks_to_analyse) for run_dir in run_dirs])


# reset index ('sample') to column
all_runs['run_ID'] = all_runs['run'].str[-8:]
all_runs = all_runs.reset_index().rename(columns={'index': 'sample'}) # reset index ('sample') to column

# Standardise numbers of taxa by sample
#all_runs_std = standardise_samples(all_runs, Ranks_to_analyse)


if os.path.exists(params_file):
    params_table = pd.read_csv(params_file, sep=',', index_col=0)
else:
    print(f"{params_files} not found")

full_table = pd.merge(all_runs, params_table, left_on='run_ID', right_on='Unique_ID')

# Drop the extra run_ID columns as it's not needed anymore
full_table.drop(columns=['run_ID'], inplace=True)


column_to_summarise=["Species","Genus","Family","Phylum"]
#full_table = add_group_statistics(full_table, column_to_group="run", columns=column_to_summarise)

# Reorder columns to place 'sample' and 'run' first
full_table = full_table[['sample', 'run','trunclenf','trunclenr'] + [col for col in full_table.columns if col not in ['sample', 'run','trunclenf','trunclenr']]]

# Save the DataFrame as a CSV file with the new filename
full_table.to_csv(os.path.join(outdir, "full_table.csv"), index=False)
