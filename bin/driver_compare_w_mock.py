#!/usr/bin/env python3
import subprocess
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt

# example usage
# ./driver_compare_w_mock.py -M 16 -R Genus -X 50 --true true_composition_mock16.csv --runs run_180132 run_176132 run_174132 run_180130 run_174130 run_178132 run_172132 run_176130 run_178131


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
    parser.add_argument("-X", type=float, required=True, help="Target ASV read number filter value")
    parser.add_argument("--true", required=True, help="true_composition.csv file")
    parser.add_argument("--runs", nargs="+", required=True, help="List of run IDs (e.g. run_123456)")
    parser.add_argument("--out", default="./output/f1_scatter.png", help="Output plot filename")
    args = parser.parse_args()

    aaa_vals, bbb_vals, f1_vals = [], [], []

    for run_id in args.runs:
        x, f1 = run_compare(args.M, run_id, args.R, args.true)

        # Find index of closest value to X
        idx = int(np.argmin([abs(val - args.X) for val in x]))
        x_prime, f1_prime = x[idx], f1[idx]

        # Extract aaa, bbb from run id (last 6 digits)
        suffix = run_id[-6:]
        aaa, bbb = int(suffix[:3]), int(suffix[3:])

        aaa_vals.append(aaa)
        bbb_vals.append(bbb)
        f1_vals.append(f1_prime)

        print(f"{run_id}: closest X={x_prime}, f1={f1_prime}")

    # Scatter plot
    plt.scatter(aaa_vals, bbb_vals, c=f1_vals, cmap="viridis", s=80)
    plt.colorbar(label="f1 score")
    plt.xlabel("aaa (first 3 digits of run suffix)")
    plt.ylabel("bbb (last 3 digits of run suffix)")
    plt.title(f"Runs at X={args.X}")
    plt.savefig(args.out, dpi=300)

if __name__ == "__main__":
    main()

