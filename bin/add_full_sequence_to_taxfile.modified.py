#!/usr/bin/env python3
import pandas as pd
import sys, os, gzip

if len(sys.argv) != 4:
    exit("Usage: add_full_sequence_to_taxfile.py <ASV_tax.tsv(.gz)> <ASV_seqs.fasta(.gz)> <outfile.tsv(.gz)>")

# Read taxonomy table
taxfile = sys.argv[1]
tax = pd.read_csv(taxfile, sep="\t", header=0, compression="infer")
tax.drop(columns="sequence", inplace=True)

# Read FASTA file
seqs = pd.DataFrame(columns=["id", "sequence"])
seq = ""
name = ""
fasta_file = sys.argv[2]
open_func = gzip.open if fasta_file.endswith(".gz") else open
with open_func(fasta_file, "rt") as reader:
    for line in reader:
        if line.startswith(">"):
            if seq and name:
                seqs = seqs.append({"id": name, "sequence": seq}, ignore_index=True)
                seq = ""
            name = line.lstrip(">").rstrip()
        else:
            seq += line.rstrip("\n")
if seq and name:
    seqs = seqs.append({"id": name, "sequence": seq}, ignore_index=True)

# Merge and write output
tax = tax.set_index("ASV_ID").join(seqs.set_index("id"), how="outer")
outfile = sys.argv[3]
compression = "gzip" if outfile.endswith(".gz") else None
tax.to_csv(outfile, sep="\t", na_rep="", index_label="ASV_ID", compression=compression)

