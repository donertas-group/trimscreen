#!/usr/bin/env python3
# example usage
# ./driver_compare_w_mock.allruns.py -D schirmer2015 --out_suffix .2 
# ./driver_compare_w_mock.allruns.py -D mock13-15 -R Genus --true /scratch/shire/data/nj/raw_data/published/mockrobiota/mock-13/true_composition_mock_13_14_15.csv
import subprocess
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

def run_compare(D, run_id, R, out_suffix, true_file):
    """Call _compare_w_true.py and capture (x, f1) output."""
    try:
        result = subprocess.run(
            ["./_compare_w_true.py", "-D", D, "-r", run_id, "-R", R, "--out_suffix", out_suffix, "--true", true_file],
            capture_output=True,
            text=True,
            check=True
        )
        # Try to parse JSON
        data = json.loads(result.stdout.strip())
        f1 = data.get("f1", [])
        f1_mean = data.get("f1_mean", [])
        return f1, f1_mean

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] _compare_w_true.py failed for run {run_id}: {e.stderr.strip()}")
        return None, None
    except json.JSONDecodeError as e:
        print(f"[ERROR] Invalid JSON output for run {run_id}: {e}")
        print("Raw output was:", result.stdout.strip())
        return None, None


def main():
    import pandas as pd
    parser = argparse.ArgumentParser(description="Aggregate _compare_w_true.py results across runs")
    parser.add_argument("-D", required=True, help="Mock dataset number")
    parser.add_argument("-R", default="Genus", help="Taxnomic rank")
    parser.add_argument("--true", default="true_composition.csv", help="true_composition.csv file")
    parser.add_argument("--out", default="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output", help="Output plot filename")
    parser.add_argument("--out_suffix", default="", help="suffix string of pipeline output dir")
    args = parser.parse_args()

    fu_tb_path = f"/scratch/shire/ssd/pipeline/16s_nf_pipeline/{args.D}/output{args.out_suffix}/compare_runs/full_table.csv"
    fu_tb = pd.read_csv(fu_tb_path)
    
    fi_tb_path = f"/scratch/shire/ssd/pipeline/16s_nf_pipeline/{args.D}/output{args.out_suffix}/compare_runs/filtered_table.csv"
    fi_tb = pd.read_csv(fi_tb_path)

    runs = fu_tb['run'].astype(str).unique().tolist()
    runs_filtered = fi_tb['run'].astype(str).unique().tolist()

    # Store data for all runs and filtered runs separately
    data_all = []
    data_filtered = []

    f1_filename = f"f1_scores_{args.D}{args.out_suffix}_{args.R}.txt"

    with open(os.path.join(args.out, f1_filename), "w") as f:
        # add header
        f.write("run_id\ttrunclenf\ttrunclenr\tf1_score\tf1_score_mean\n")

        for run_id in runs:
            f1, f1_mean = run_compare(args.D, run_id, args.R, args.out_suffix, args.true)
            if not f1:
                continue  # Skip this run safely
            print(f1)
            # Extract lenf, lenr from run id (last 6 digits)
            suffix = run_id[-6:]
            lenf, lenr = int(suffix[:3]), int(suffix[3:])

            # Save results
            f.write(f"{run_id}\t{lenf}\t{lenr}\t{f1}\t{f1_mean}\n")

            # Add to appropriate dataset for plotting
            record = {"run_id": run_id, "lenf": lenf, "lenr": lenr, "f1": f1, "f1_mean": f1_mean}
            data_all.append(record)

            if run_id in runs_filtered:
                data_filtered.append(record)

    print(f"All results saved to {os.path.join(args.out, f1_filename)}")

    # --- Prepare data for plotting ---
    import pandas as pd
    df_all = pd.DataFrame(data_all)
    df_filtered = pd.DataFrame(data_filtered)

    # --- Define color normalization ---
    norm = plt.Normalize(vmin=df_all["f1"].min(), vmax=df_all["f1"].max())

    # --- Plotting following script 2 style ---
    plt.figure(figsize=(8, 6))

    # Plot unfiltered runs (semi-transparent)
    plt.scatter(
        df_all["lenf"],
        df_all["lenr"],
        c=df_all["f1"],
        cmap="viridis",
        norm=norm,
        alpha=0.5,
        edgecolor="none",
        label="Dropped"
    )

    # Plot filtered runs (solid color, black border)
    plt.scatter(
        df_filtered["lenf"],
        df_filtered["lenr"],
        c=df_filtered["f1"],
        cmap="viridis",
        norm=norm,
        alpha=1.0,
        edgecolor="black",
        linewidth=0.5,
        label="Kept"
    )

    # Colorbar + labels
    plt.colorbar(label=f"F1 score at {args.R}")
    plt.xlabel("Forward read length")
    plt.ylabel("Reverse read length")
    plt.title(f"Mock {args.D}")
    plt.legend()

    # Save plot
    outfile = os.path.join(args.out, f"f1_scatter_{args.D}{args.out_suffix}_{args.R}.png")
    plt.savefig(outfile, dpi=300)
    plt.close()

    print(f"Plot saved to {outfile}")

if __name__ == "__main__":
    main()

