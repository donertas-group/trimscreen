#!/usr/bin/env python3
import subprocess
import json
import argparse
import os
import pandas as pd

def run_compare(D, run_id, R, out_suffix, true_file):
    """Call _compare_w_true.py and capture per-sample JSON output."""
    try:
        # _compare_w_true.py now prints one JSON per sample
        result = subprocess.run(
            ["./_compare_w_true.py", "-D", D, "-r", run_id, "-R", R, "--out_suffix", out_suffix, "--true", true_file],
            capture_output=True,
            text=True,
            check=True
        )
        # Each line is JSON per sample
        sample_results = [json.loads(line) for line in result.stdout.strip().splitlines()]
        return sample_results

    except (subprocess.CalledProcessError, json.JSONDecodeError) as e:
        print(f"[ERROR] run {run_id} failed: {e}")
        return []

def main():
    parser = argparse.ArgumentParser(description="Aggregate _compare_w_true.py results across runs")
    parser.add_argument("-D", required=True, help="Mock dataset")
    parser.add_argument("-R", default="Genus", help="Taxonomic rank")
    parser.add_argument("--true", default="true_composition.csv", help="True composition CSV file")
    parser.add_argument("--out", default="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output", help="Output directory")
    parser.add_argument("--out_suffix", default="", help="Pipeline output suffix")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # Load run lists from full table and filtered table
    fu_tb_path = f"/scratch/shire/ssd/pipeline/16s_nf_pipeline/{args.D}/output{args.out_suffix}/compare_runs/full_table.csv"
    fi_tb_path = f"/scratch/shire/ssd/pipeline/16s_nf_pipeline/{args.D}/output{args.out_suffix}/compare_runs/filtered_table.csv"

    fu_tb = pd.read_csv(fu_tb_path)
    fi_tb = pd.read_csv(fi_tb_path)

    runs_all = fu_tb['run'].astype(str).unique()
    runs_filtered = set(fi_tb['run'].astype(str).unique())

    # Output CSV file
    output_csv = os.path.join(args.out, f"f1_scores_{args.D}{args.out_suffix}_{args.R}.csv")

    all_records = []

    for run_id in runs_filtered: # alt: runs_all
        sample_data = run_compare(args.D, run_id, args.R, args.out_suffix, args.true)
        if not sample_data:
            continue

        for rec in sample_data:
            rec["run_id"] = run_id
#            rec["filtered"] = run_id in runs_filtered
            all_records.append(rec)

    if not all_records:
        print("No data collected. Exiting.")
        return

    # Save all per-sample results to CSV
    df_all = pd.DataFrame(all_records)
    df_all.to_csv(output_csv, index=False)
    print(f"Per-sample results saved to {output_csv}")

if __name__ == "__main__":
    main()

