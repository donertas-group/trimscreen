#!/usr/bin/env python3

# example usage
# ./driver_reps_similarity.allruns.py -D hc227_v3v4 --f1_file f1_scores_hc227_v3v4_Genus.txt --median_distance_file median_distances_per_sample.hc227_v3v4.csv
# ./driver_reps_similarity.allruns.py -D schirmer2015 --out_suffix .2 --f1_file f1_scores_schirmer2015.2_Genus.txt --median_distance_file median_distances_per_sample.schirmer2015.2.csv

# choices=["hc227_v3v4","mock13-15","mock03_05","mock20_22","mock21_23","tourlousse2022","schirmer2015"]

import subprocess
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

def run_compare(D, run_id, R, out_suffix):
    """Call _betadiv_reps.py and capture (run, similarity) output."""
    try:
        result = subprocess.run(
            ["./_betadiv_reps.py", "-D", D, "-r", run_id, "-R", R, "--out_suffix", out_suffix],
            capture_output=True,
            text=True,
            check=True
        )
        # Try to parse JSON
        data = json.loads(result.stdout.strip())
        similarity = data.get("similarity", [])
        
        # Validate the content
        #if not isinstance(similarity, float):
        #    print(f"[WARN] Empty or invalid output for run {run_id}. Skipping.")
        #    return None
        return similarity

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] _betadiv_reps.py failed for run {run_id}: {e.stderr.strip()}")
        return None
    except json.JSONDecodeError as e:
        print(f"[ERROR] Invalid JSON output for run {run_id}: {e}")
        print("Raw output was:", result.stdout.strip())
        return None


def main():
    import pandas as pd
    parser = argparse.ArgumentParser(description="Aggregate compare_w_mock.py results across runs")
    parser.add_argument("-D", required=True, help="Dataset name")#choices=["hc227_v3v4","mock13-15","mock03_05","mock20_22","mock21_23","tourlousse2022","schirmer2015"],
    parser.add_argument("-R", default="Genus", help="Taxnomic rank")
    parser.add_argument("--out", default="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output", help="Output plot filename")
    parser.add_argument("--out_suffix", default="", help="suffix string of pipeline output dir")
    parser.add_argument("--f1_file", type=str, required=False, help="txt file listing f1 scores and runIDs")
    parser.add_argument("--median_distance_file", type=str, required=False, help="csv file listing median distances and runIDs")
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

    # --- Prepare data for plotting ---

    for run_id in runs:
        similarity = run_compare(args.D, run_id, args.R, args.out_suffix)

        # take metrics from the first sample among the replicates
        shannon = fu_tb.loc[fu_tb['run'] == run_id, 'shannon_Genus'].iloc[0]
        preads = fu_tb.loc[fu_tb['run'] == run_id, 'retained_reads_percent'].iloc[0]
        nasvs = fu_tb.loc[fu_tb['run'] == run_id, 'nasvs'].iloc[0]

        # Extract lenf, lenr from run id (last 6 digits)
        suffix = run_id[-6:]
        lenf, lenr = int(suffix[:3]), int(suffix[3:])

        # Add to appropriate dataset
        record = {"run_id": run_id, "lenf": lenf, "lenr": lenr, "retained_reads_percent": preads, "nasvs": nasvs, "similarity": similarity, "shannon_Genus": shannon}
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
        alpha=0.2,
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
    outfile = os.path.join(args.out, f"reps_similarity_{args.D}{args.out_suffix}_{args.R}.png")
    plt.savefig(outfile, dpi=300)
    plt.close()

    print(f"Plot saved to {outfile}")


    # ---- Merge with f1 file and plot ----
    if args.f1_file and args.median_distance_file:
        f1_df = pd.read_csv(os.path.join(args.out, args.f1_file), sep='\t')
        med_dist = pd.read_csv(os.path.join(args.out, args.median_distance_file))

        # ensure column names match (case-insensitive check)
        expected_cols1 = {"run_id", "f1_score"}
        if not expected_cols1.issubset(f1_df.columns):
            raise ValueError(f"f1_file must have columns {expected_cols1}, but found {f1_df.columns.tolist()}")

        expected_cols2 = {"run_id", "median_distance"}
        if not expected_cols2.issubset(med_dist.columns):
            raise ValueError(f"f1_file must have columns {expected_cols2}, but found {med_dist.columns.tolist()}")

        merged = pd.merge(df_all, f1_df, on="run_id", how="inner")
        merged = pd.merge(merged, med_dist, on="run_id", how="inner")

        merged.to_csv(os.path.join(args.out,f"merged_table.{args.D}{args.out_suffix}.csv"), sep=",",index=False)

        if merged.empty:
            print("[WARN] No overlapping run_ids found between df_all and f1_file.")
        else:
            xvars = ["shannon_Genus", "retained_reads_percent", "similarity", "median_distance"]

            for x_var in xvars:
                if x_var not in merged.columns:
                    print(f"[WARN] Column '{x_var}' not found in merged DataFrame. Skipping.")
                    continue
                merged[x_var] = pd.to_numeric(merged[x_var])
                plt.figure(figsize=(8, 6))
                plt.scatter(merged[x_var], merged["f1_score"], alpha=0.2)
                plt.xlabel(x_var.capitalize().replace("_", " "))
                plt.ylabel("F1 Score")
                plt.title(f"{x_var.capitalize()} vs F1 Score ({args.D}, Rank={args.R})")
                plt.grid(True)

                scatterfile = os.path.join(args.out, f"{x_var}_vs_f1_{args.D}{args.out_suffix}.png")
                os.makedirs(os.path.dirname(scatterfile), exist_ok=True)
                plt.savefig(scatterfile, dpi=300, bbox_inches="tight")
                plt.close()

                print(f"[INFO] Scatter plot saved to {scatterfile}")

if __name__ == "__main__":
    main()

