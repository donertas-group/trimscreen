#!/usr/bin/env python3
import subprocess
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os

# example usage
# ./driver_compare_w_mock.py -M 16 -R Genus -X 50 --true true_composition.csv --runs run_180132 run_176132 run_174132 run_180130 run_174130 run_178132 run_172132 run_176130 run_178131


def run_compare(M, run_id, R, true_file):
    """Call compare_w_mock.py and capture (x, f1) output."""
    # Run compare_w_mock.py as a subprocess
    result = subprocess.run(
        ["./compare_w_mock.py", "-M", M, "-r", run_id, "-R", R, "--true", true_file],
        capture_output=True,
        text=True,
        check=True
    )
    # Assume compare_w_mock.py prints JSON like: {"x": [...], "f1": [...]}
    data = json.loads(result.stdout.strip())
    return data["x"], data["f1"]

def main():
    parser = argparse.ArgumentParser(description="Aggregate compare_w_mock.py results across runs")
    parser.add_argument("-M", required=True, help="Mock dataset number")
    parser.add_argument("-R", required=True, help="Taxnomic rank")
    parser.add_argument("-X", type=int, required=True, help="Target ASV read number filter value")
    parser.add_argument("--true", default="true_composition.csv", help="true_composition.csv file")
    parser.add_argument("--runs", nargs="+", required=True, help="List of run IDs (e.g. run_123456)")
    parser.add_argument("--out", default="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output", help="Output plot filename")
    args = parser.parse_args()

    lenf_vals, lenr_vals, f1_vals = [], [], []

    with open(os.path.join(args.out, f"f1_scores_mock{args.M}_{args.R}_min{args.X}.txt"), "w") as f:
        # add header
        f.write("run_id\ttrunclenf\ttrunclenr\tasv_abund_threshold\tf1_score\n")

        for run_id in args.runs:
            x, f1 = run_compare(args.M, run_id, args.R, args.true)

            # Find index of closest and no smaller value to X
            idx = int(np.argmin([abs(val - args.X) for val in x]))
            x = np.array(x) 
            mask = x >= args.X
            if np.any(mask):
                idx = np.where(mask)[0][np.argmin(x[mask] - args.X)]
            else:
                idx = None  # no value >= X

            x_prime, f1_prime = x[idx], f1[idx]

            # Extract lenf, lenr from run id (last 6 digits)
            suffix = run_id[-6:]
            lenf, lenr = int(suffix[:3]), int(suffix[3:])

            lenf_vals.append(lenf)
            lenr_vals.append(lenr)
            f1_vals.append(f1_prime)

            f.write(f"{run_id}\t{lenf}\t{lenr}\t{x_prime}\t{f1_prime}\n")
           # print(f"{run_id}: closest X={x_prime}, f1={f1_prime}")

    print(f"All results saved to {os.path.join(args.out, 'f1_scores.txt')}")

    # Scatter plot
    plt.scatter(lenf_vals, lenr_vals, c=f1_vals, cmap="viridis", s=40)
    plt.colorbar(label=f"f1 score at {args.R}")
    plt.xlabel("Forward read length")
    plt.ylabel("Reverse read length")
    plt.title(f"Mock {args.M}, min {args.X} reads/ASV")

    outfile = os.path.join(args.out, f"f1_scatter_mock{args.M}_{args.R}_min{args.X}.rarefied.png")
    plt.savefig(outfile, dpi=300)

if __name__ == "__main__":
    main()

