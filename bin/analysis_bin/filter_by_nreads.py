#!/usr/bin/env python3

# usage
# ./filter_by_nreads.py -i merged_table.tourlousse2022.1.csv -o merged_table.filtered.tourlousse2022.1.csv


import csv
import argparse
from collections import defaultdict
from pathlib import Path


MIN_NREADS = 2000
WDIR = Path("/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output")



def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter runs where any sample has nreads < 5000"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input CSV file (e.g. merged_table.sam.csv)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV file (e.g. merged_table.filtered.sam.csv)"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = WDIR / args.input
    output_path = WDIR / args.output

    with input_path.open(newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames

    # Collect nreads per run
    run_to_nreads = defaultdict(list)
    for row in rows:
        run = row["run"]
        nreads = int(row["nreads"])
        run_to_nreads[run].append(nreads)

    # Keep runs where the minimum nreads >= threshold
    valid_runs = {
        run for run, values in run_to_nreads.items()
        if min(values) >= MIN_NREADS
    }

    # Filter rows
    filtered_rows = [
        row for row in rows
        if row["run"] in valid_runs
    ]

    # Write output
    with output_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(filtered_rows)


if __name__ == "__main__":
    main()

