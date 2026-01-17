#!/usr/bin/env python3

# example usage
# ./driver_reps_similarity.allruns.py -D hc227_v3v4 --f1_file f1_scores_hc227_v3v4_Genus.txt --median_distance_file median_distances_per_sample.hc227_v3v4.csv
# ./driver_reps_similarity.allruns.py -D schirmer2015 --out_suffix .2 --f1_file f1_scores_schirmer2015.2_Genus.txt --median_distance_file median_distances_per_sample.schirmer2015.2.csv


import subprocess
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd


def main():
    import pandas as pd
    parser = argparse.ArgumentParser(description="Aggregate compare_w_mock.py results across runs")
    parser.add_argument("-D", required=True, help="Dataset name")#choices=["hc227_v3v4","mock13-15","mock03_05","mock20_22","mock21_23","tourlousse2022","schirmer2015"],
    parser.add_argument("-R", default="Genus", help="Taxnomic rank")
    parser.add_argument("--out", default="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output", help="Output plot filename")
    parser.add_argument("--out_suffix", default="", help="suffix string of pipeline output dir")
    parser.add_argument("--f1_file", type=str, required=True, help="txt file listing f1 scores and runIDs")
    args = parser.parse_args()

    fu_tb_path = f"/scratch/shire/ssd/pipeline/16s_nf_pipeline/{args.D}/output{args.out_suffix}/compare_runs/full_table.csv"
    full_tab = pd.read_csv(fu_tb_path)
    
    fi_tb_path = f"/scratch/shire/ssd/pipeline/16s_nf_pipeline/{args.D}/output{args.out_suffix}/compare_runs/filtered_table.csv"
    filt_tab = pd.read_csv(fi_tb_path)

    f1_df = pd.read_csv(os.path.join(args.out, args.f1_file))

    runs = full_tab['run'].astype(str).unique().tolist()
    runs_filtered = filt_tab['run'].astype(str).unique().tolist()
    

    # ensure column names match (case-insensitive check)
    expected_cols = {"run", "f1"}
    if not expected_cols.issubset(f1_df.columns):
        raise ValueError(f"f1_file must have columns {expected_cols1}, but found {f1_df.columns.tolist()}")

    merged = pd.merge(full_tab, f1_df, on=["run","sample"], how="inner")
    merged.to_csv(os.path.join(args.out,f"merged_table.{args.D}{args.out_suffix}.csv"), sep=",",index=False)


if __name__ == "__main__":
    main()
