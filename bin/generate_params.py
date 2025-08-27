#!/usr/bin/env python3
import argparse
import itertools
import pandas as pd
import os
import sys
import numpy as np

########## Introduction ###########
# The mandatory inputs (marker_size_min, minimum_overlap, step, read_length, FW_primer_len, RV_primer_len) are required.
# The optional inputs (trunclenf_range, trunclenr_range) will override the mandatory logic if both are provided.
# The output is the same in both cases: a summary_params_settings.csv file.

#Mandatory algorithm (default):

#python script.py --marker_size_min 250 --minimum_overlap 20 --step_size 10 \
#--read_length 300 --FW_primer_len 17 --RV_primer_len 21 -o results/

#Optional algorithm (when both -f and -r are provided):

#python script.py --marker_size_min 250 --minimum_overlap 20 --step_size 10 \
#--read_length 300 --FW_primer_len 17 --RV_primer_len 21 \
#-f 100:10:200 -r 100,150,200 -o results/

###################################

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate trimming parameter combinations based on user input."
    )

    # Required arguments (for mandatory algorithm)
    parser.add_argument("--marker_size_min", type=int, required=True,
                        help="Minimum marker size.")
    parser.add_argument("--minimum_overlap", type=int, 
                        help="Minimum required overlap.")
    parser.add_argument("--step_size", type=int, required=True,
                        help="Step size for generating values.")
    parser.add_argument("--read_length", type=int, required=True,
                        help="Maximum read length.")
    parser.add_argument("--FW_primer_len", type=int, required=True,
                        help="Forward primer length.")
    parser.add_argument("--RV_primer_len", type=int, required=True,
                        help="Reverse primer length.")
    parser.add_argument("--outdir", "-o", type=str, required=True,
                        help="Output directory.")

    # Optional arguments 
    parser.add_argument("--trunclenf_range", "-f", type=str,
                        help="Format: 'min:step:max', 'value1,value2', or 'value'.")
    parser.add_argument("--trunclenr_range", "-r", type=str,
                        help="Format: 'min:step:max', 'value1,value2', or 'value'.")

    return parser.parse_args()

def parse_range(value):
    """Parse range input in the format min:step:max or comma-separated values."""
    if ':' in value:  # Range format
        try:
            start, step, end = map(int, value.split(':'))
            return [start + i * step for i in range((end - start) // step + 1)]
        except ValueError:
            raise ValueError(f"Invalid range format: {value}")
    elif ',' in value:  # Comma-separated
        return [int(v.strip()) for v in value.split(',')]
    else:  # Single value
        return [int(value)]

def run_optional_algorithm(trunclenf_range, trunclenr_range, outdir):
    f_values = parse_range(trunclenf_range)
    r_values = parse_range(trunclenr_range)
    combinations = list(itertools.product(f_values, r_values))
    run_data = [{"runID": f"run_{f}{r}", "trunclenf": f, "trunclenr": r} for f, r in combinations]

    os.makedirs(outdir, exist_ok=True)
    output_path = os.path.join(outdir, "summary_params_settings.csv")
    pd.DataFrame(run_data).to_csv(output_path, index=False)

def run_mandatory_algorithm(args):
    total_min = args.marker_size_min + args.minimum_overlap - 10 # allow 10bp buffer for min_marker_size
    min_val = 50 # consistent with dada2
    max_val = args.read_length
    step = args.step_size
    f_len = args.FW_primer_len
    r_len = args.RV_primer_len

    f_values = [x for x in (np.arange(min_val, max_val + 1, step) - f_len) if x > 0]
    r_values = [x for x in (np.arange(min_val, max_val + 1, step) - r_len) if x > 0]

    valid_combinations = [
        (f, r) for f, r in itertools.product(f_values, r_values)
        if f + r >= total_min
    ]

    run_data = [{"runID": f"run_{f:03d}{r:03d}", "trunclenf": f, "trunclenr": r} for f, r in valid_combinations]

    os.makedirs(args.outdir, exist_ok=True)
    output_path = os.path.join(args.outdir, "summary_params_settings.csv")
    pd.DataFrame(run_data).to_csv(output_path, index=False)

def main():
    args = parse_args()

    if args.trunclenf_range and args.trunclenr_range:
        run_optional_algorithm(args.trunclenf_range, args.trunclenr_range, args.outdir)
    else:
        run_mandatory_algorithm(args)

if __name__ == "__main__":
    sys.exit(main())

