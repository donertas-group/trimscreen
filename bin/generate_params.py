#!/usr/bin/env python
import argparse
import itertools
import pandas as pd
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate trimming parameter combinations with constraints."
    )
    parser.add_argument("--marker_size_min", type=int, required=True,
                        help="Minimum marker size.")
    parser.add_argument("--minimum_overlap", type=int, required=True,
                        help="Minimum required overlap.")
    parser.add_argument("--step_size", type=int, required=True,
                        help="Step size for generating f and r.")
    parser.add_argument("--read_length", type=int, required=True,
                        help="Maximum read length.")
    parser.add_argument("--outdir", "-o", type=str, required=True,
                        help="Output directory.")

    return parser.parse_args()

def main():
    args = parse_args()

    total_min = args.marker_size_min + args.minimum_overlap
    min_val = 50 #args.minimum_overlap
    max_val = args.read_length
    step = args.step_size

    # Generate valid values for f and r
    f_values = list(range(min_val, max_val + 1, step))
    r_values = list(range(min_val, max_val + 1, step))

    # Filter combinations based on constraints
    valid_combinations = [
        (f, r) for f, r in itertools.product(f_values, r_values)
        if f + r >= total_min
    ]

    run_data = [
        {"runID": f"run_{f:03d}{r:03d}", "trunclenf": f, "trunclenr": r}
        for f, r in valid_combinations
    ]

    os.makedirs(args.outdir, exist_ok=True)
    output_path = os.path.join(args.outdir, "summary_params_settings.csv")

    df = pd.DataFrame(run_data)
    df.to_csv(output_path, index=False)

if __name__ == "__main__":
    sys.exit(main())

