#!/usr/bin/env python3
# example usage
# ./driver_compare_w_mock.allruns.py -M mock13-15 -R Genus --true /scratch/shire/data/nj/raw_data/published/mockrobiota/mock-13/true_composition_mock_13_14_15.csv
import subprocess
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

def run_compare(M, run_id, R, true_file):
    """Call compare_w_true.py and capture (x, f1) output."""
    try:
        result = subprocess.run(
            ["./compare_w_true.py", "-M", M, "-r", run_id, "-R", R, "--true", true_file],
            capture_output=True,
            text=True,
            check=True
        )
        # Try to parse JSON
        data = json.loads(result.stdout.strip())
#        x = data.get("x", [])
        f1 = data.get("f1", [])
        # Validate the content
   #     if not isinstance(f1, float):
   #         print(f"[WARN] Empty or invalid output for run {run_id}. Skipping.")
   #         return None
        return f1

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] compare_w_true.py failed for run {run_id}: {e.stderr.strip()}")
        return None
    except json.JSONDecodeError as e:
        print(f"[ERROR] Invalid JSON output for run {run_id}: {e}")
        print("Raw output was:", result.stdout.strip())
        return None


def main():
    import pandas as pd
    parser = argparse.ArgumentParser(description="Aggregate compare_w_true.py results across runs")
    parser.add_argument("-M", required=True, help="Mock dataset number")
    parser.add_argument("-R", required=True, help="Taxnomic rank")
    #parser.add_argument("-X", type=int, required=True, help="Target ASV read number filter value")
    parser.add_argument("--true", default="true_composition.csv", help="true_composition.csv file")
    #parser.add_argument("--runs", nargs="+", required=True, help="List of run IDs (e.g. run_123456)")
    parser.add_argument("--out", default="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output", help="Output plot filename")
    args = parser.parse_args()

    fu_tb_path = f"/scratch/shire/ssd/pipeline/16s_nf_pipeline/{args.M}/output/compare_runs/full_table.csv"
    fu_tb = pd.read_csv(fu_tb_path)
    
    fi_tb_path = f"/scratch/shire/ssd/pipeline/16s_nf_pipeline/{args.M}/output/compare_runs/filtered_table.csv"
    fi_tb = pd.read_csv(fi_tb_path)

    runs = fu_tb['run'].astype(str).unique().tolist()
    runs_filtered = fi_tb['run'].astype(str).unique().tolist()

    # Store data for all runs and filtered runs separately
    data_all = []
    data_filtered = []

    with open(os.path.join(args.out, f"f1_scores_{args.M}_{args.R}.txt"), "w") as f:
        # add header
        f.write("run_id\ttrunclenf\ttrunclenr\tasv_abund_threshold\tf1_score\n")

        for run_id in runs:
            f1 = run_compare(args.M, run_id, args.R, args.true)
            if not f1:
                continue  # Skip this run safely
            print(f1)
            # Extract lenf, lenr from run id (last 6 digits)
            suffix = run_id[-6:]
            lenf, lenr = int(suffix[:3]), int(suffix[3:])

            # Save results
            f.write(f"{run_id}\t{lenf}\t{lenr}\t{f1}\n")

            # Add to appropriate dataset
            record = {"run_id": run_id, "lenf": lenf, "lenr": lenr, "f1": f1}
            data_all.append(record)
            if run_id in runs_filtered:
                data_filtered.append(record)

    print(f"All results saved to {os.path.join(args.out, 'f1_scores.txt')}")

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
    plt.title(f"Mock {args.M}")
    plt.legend()

    # Save plot
    outfile = os.path.join(args.out, f"f1_scatter_{args.M}_{args.R}.png")
    plt.savefig(outfile, dpi=300)
    plt.close()

    print(f"Plot saved to {outfile}")

if __name__ == "__main__":
    main()

