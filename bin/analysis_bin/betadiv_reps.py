#!/usr/bin/env python3

# example usage:
# ./betadiv_reps.py -D hc227_v3v4 -r run_225225 -R Genus

import argparse
import pandas as pd
import os
import numpy as np
import glob
import json
from scipy.spatial.distance import pdist, squareform

def parse_args():
    parser = argparse.ArgumentParser(description="Compare 16S data with true mock composition.")
    parser.add_argument("-D", "--dataset", required=True, type=str, help="Dataset name")
    parser.add_argument("-r", "--run", required=True, type=str, help="Run ID")
    parser.add_argument("-R", "--rank", required=True, type=str, help="Taxonomic rank to compare (e.g., Genus)")
    return parser.parse_args()

def ruzicka(u, v):
    """Compute the Ruzicka (abundance-based Jaccard) similarity between two vectors."""
    u, v = np.asarray(u), np.asarray(v)
    numerator = np.minimum(u, v).sum()
    denominator = np.maximum(u, v).sum()
    if denominator == 0:  # avoid division by zero
        return 0.0
    return 1-(numerator / denominator)

def calculate_mean_similarity(dataset, run_id, rank):
    result_dir = "/scratch/shire/ssd/pipeline/16s_nf_pipeline"  # Modify if needed

    asv_tax_file = os.path.join(result_dir, f"{dataset}", "output/runs", run_id, "dada2/ASV_tax.silva_138_2.tsv.gz")
    asv_table_file = os.path.join(result_dir, f"{dataset}", "output/runs", run_id, "dada2/ASV_table.tsv.gz")

    asv_tax = pd.read_csv(asv_tax_file, sep='\t')
    asv_table = pd.read_csv(asv_table_file, sep='\t')

    merged = pd.merge(asv_table, asv_tax, on='ASV_ID')

    # Drop NA or empty rank entries
    merged = merged[merged[rank].notna() & (merged[rank].astype(str).str.strip() != '')]

    if merged.empty:
        return None

    meta_cols = ["ASV_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "confidence","sequence"]
    meta_cols = [c for c in meta_cols if c in merged.columns]

    # Identify sample columns
    sample_cols = [c for c in merged.columns if c not in meta_cols]

    # Aggregate counts at chosen rank
    agg = merged.groupby(rank)[sample_cols].sum()

    # Transpose for distance computation
    agg_T = agg.T  # samples as rows

    # Compute distance
    dist_matrix = pdist(agg_T, metric=ruzicka)
    dist_df = pd.DataFrame(squareform(dist_matrix), index=agg_T.index, columns=agg_T.index)

    # Compute mean dissimilarity & similarity
    mean_dist = dist_df.values[np.triu_indices_from(dist_df, k=1)].mean()
    similarity = 1 - mean_dist
        
    return similarity


def main():
    args = parse_args()
    similarity = calculate_mean_similarity(args.dataset, args.run, args.rank)
    
    output = {"run": args.run, "similarity": similarity}
    print(json.dumps(output))

if __name__ == "__main__":
    main()

