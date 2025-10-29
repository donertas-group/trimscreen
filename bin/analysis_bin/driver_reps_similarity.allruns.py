#!/usr/bin/env python3
# example usage
# ./driver_reps_similarity.allruns.py -D hc227_v3v4 -R Genus 
import subprocess
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

def run_compare(D, run_id, R):
    """Call betadiv_reps.py and capture (run, similarity) output."""
    try:
        result = subprocess.run(
            ["./betadiv_reps.py", "-D", D, "-r", run_id, "-R", R],
            capture_output=True,
            text=True,
            check=True
        )
        # Try to parse JSON
        data = json.loads(result.stdout.strip())
        similarity = data.get("similarity", [])
        
        # Validate the content
        if not isinstance(similarity, float):
            print(f"[WARN] Empty or invalid output for run {run_id}. Skipping.")
            return None
        return similarity

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] betadiv_reps.py failed for run {run_id}: {e.stderr.strip()}")
        return None
    except json.JSONDecodeError as e:
        print(f"[ERROR] Invalid JSON output for run {run_id}: {e}")
        print("Raw output was:", result.stdout.strip())
        return None


def main():
    import pandas as pd
    parser = argparse.ArgumentParser(description="Aggregate compare_w_mock.py results across runs")
    parser.add_argument("-D", required=True, choices=["hc227_v3v4","mock13-15","mock20_22","mock21_23"], help="Dataset name")
    parser.add_argument("-R", required=True, help="Taxnomic rank")
    parser.add_argument("--out", default="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output", help="Output plot filename")
    args = parser.parse_args()

    fu_tb_path = f"/scratch/shire/ssd/pipeline/16s_nf_pipeline/{args.D}/output/compare_runs/full_table.csv"
    fu_tb = pd.read_csv(fu_tb_path)
    
    fi_tb_path = f"/scratch/shire/ssd/pipeline/16s_nf_pipeline/{args.D}/output/compare_runs/filtered_table.csv"
    fi_tb = pd.read_csv(fi_tb_path)

    runs = fu_tb['run'].astype(str).unique().tolist()
    runs_filtered = fi_tb['run'].astype(str).unique().tolist()

    # Store data for all runs and filtered runs separately
    data_all = []
    data_filtered = []

    # --- Prepare data for plotting ---

    for run_id in runs:
        similarity = run_compare(args.D, run_id, args.R)
        if not similarity:
            continue  # Skip this run safely

        # Extract lenf, lenr from run id (last 6 digits)
        suffix = run_id[-6:]
        lenf, lenr = int(suffix[:3]), int(suffix[3:])

        # Add to appropriate dataset
        record = {"run_id": run_id, "lenf": lenf, "lenr": lenr, "similarity": similarity}
        data_all.append(record)
        print(record)
        if run_id in runs_filtered:
            data_filtered.append(record)

    import pandas as pd
    df_all = pd.DataFrame(data_all)
    df_filtered = pd.DataFrame(data_filtered)

    # --- Define color normalization ---
    norm = plt.Normalize(vmin = df_all["similarity"].min(), vmax = df_all["similarity"].max())

    # --- Plotting following script 2 style ---
    plt.figure(figsize=(8, 6))

    # Plot unfiltered runs (semi-transparent)
    plt.scatter(
        df_all["lenf"],
        df_all["lenr"],
        c=df_all["similarity"],
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
        c=df_filtered["similarity"],
        cmap="viridis",
        norm=norm,
        alpha=1.0,
        edgecolor="black",
        linewidth=0.5,
        label="Kept"
    )

    # Colorbar + labels
    plt.colorbar(label=f"Similarity at {args.R}")
    plt.xlabel("Forward read length")
    plt.ylabel("Reverse read length")
    plt.title(f"{args.D}")
    plt.legend()

    # Save plot
    outfile = os.path.join(args.out, f"reps_similarity_{args.D}_{args.R}.png")
    plt.savefig(outfile, dpi=300)
    plt.close()

    print(f"Plot saved to {outfile}")

if __name__ == "__main__":
    main()

