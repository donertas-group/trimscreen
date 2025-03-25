#!/usr/bin/env python
import argparse
import itertools
import uuid
import pandas as pd
import os
import sys
def parse_args():
    parser = argparse.ArgumentParser(
        description="Creating a CSV file listing all pairs of trimming parameters to be screened.")
    parser.add_argument(
        "--trunclenf_range", "-f", 
        help="Format: 'min:step:max', 'value1,value2', or 'value'.", 
        type=str, required=True)
    parser.add_argument(
        "--trunclenr_range", "-r", 
        help="Format: 'min:step:max', 'value1,value2', or 'value'.", 
        type=str, required=True)
    parser.add_argument("--outdir", "-o", type=str, required=True,
                        help="Output directory specified by the workflow.")

    return parser.parse_args()

def parse_range(value):
    """Parse range input in the format min:step:max or comma-separated values."""
    if ':' in value:  # Range input
        try:
            start, step, end = map(int, value.split(':'))
            return [start + i * step for i in range((end - start) // step + 1)]
        except ValueError:
            raise ValueError(f"Invalid range format: {value}")
    elif ',' in value:  # Comma-separated values
        return [int(v.strip()) for v in value.split(',')]
    else:  # Single value
        return [int(value)]

def main():
    args = parse_args()
    
    f_values = parse_range(args.trunclenf_range)
    r_values = parse_range(args.trunclenr_range)
    
    combinations = list(itertools.product(f_values, r_values))
    
    # Assign unique 8-digit runIDs
    #run_data = [{"runID": str(uuid.uuid4())[:8], "trunclenf": f, "trunclenr": r} for f, r in combinations]
    run_data = [{"runID": f"run_{f}{r}", "trunclenf": f, "trunclenr": r} for f, r in combinations]

    #os.makedirs(os.path.join(args.outdir, "generate_params"), exist_ok=True)  # Ensure output directory exists
    output_path = os.path.join(args.outdir, "summary_params_settings.csv")

    df = pd.DataFrame(run_data)
    df.to_csv(output_path, index=False)

if __name__ == "__main__":
    sys.exit(main())
