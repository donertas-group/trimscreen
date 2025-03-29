import pandas as pd
import numpy as np
import os
import glob
from _process_run import filter_runs
import config

# Use the variables defined in config.py
workdir = config.workdir
run_dirs = config.run_dirs
outdir = config.outdir
classifier_dir = config.classifier_dir

full_table = pd.read_csv(os.path.join(outdir, "full_table.csv"))

# filter runs by evaluating the median and sd of choosen columns
filtered_table = filter_runs(full_table, ['retained_percent','Phylum_pasv','Genus_pasv'])
filtered_table.to_csv(os.path.join(outdir, "filtered_table.csv"), index=False)


