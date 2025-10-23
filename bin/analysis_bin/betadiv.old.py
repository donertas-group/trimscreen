#!/usr/bin/env python

# example usage: ./betadiv.py -D mock16 -R Genus
import os
import re
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from skbio.stats import subsample_counts  # skbio is efficient for rarefaction
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch
import argparse
import sys

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Calculate beta-diversity for mock dataset")
    parser.add_argument("-D", "--dataset_name", required=True, help="Dataset name")
    parser.add_argument("-R", "--rank", required=True, help="Taxonomic rank to calculate beta")
    parser.add_argument("--out", default="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output", help="Output dir")

    return parser.parse_args()

def ruzicka(u, v):
    """Compute the Ruzicka (abundance-based Jaccard) similarity between two vectors."""
    u, v = np.asarray(u), np.asarray(v)
    numerator = np.minimum(u, v).sum()
    denominator = np.maximum(u, v).sum()
    if denominator == 0:  # avoid division by zero
        return 0.0
    return 1-(numerator / denominator)


def main():
    args = parse_args()
    rank = args.rank

    # 1. Set sample ID
    sampleID = re.sub(r'(\D+)(\d+)', r'\1_\2', args.dataset_name)
    base_dir = "/scratch/shire/ssd/pipeline/16s_nf_pipeline/" + args.dataset_name + "/output"
    compare_runs_csv = os.path.join(base_dir, "compare_runs", "filtered_table.csv")

    # 2. Get run list
    df_runs = pd.read_csv(compare_runs_csv)
    run_list = df_runs['run'].dropna().unique().tolist()

    # 3. Read ASV tables and merge
    asv_tables = []
    for run in run_list:
        path = os.path.join(base_dir, "runs", run, "dada2", "ASV_table.tsv.gz")
        df = pd.read_csv(path, sep="\t", usecols=["ASV_ID", sampleID])
        df = df.rename(columns={sampleID: run})
        asv_tables.append(df)

    # 4. Full outer join on ASV_ID
    asv_table_full = asv_tables[0]
    for df in asv_tables[1:]:
        asv_table_full = pd.merge(asv_table_full, df, on="ASV_ID", how="outer")

    asv_table_full = asv_table_full.fillna(0)
    asv_table_full.set_index("ASV_ID", inplace=True)
    asv_table_full_filtered = asv_table_full[asv_table_full.max(axis=1) >= 20]

    # 5 Rarefy each run to the same sequencing depth (per column)
    # Convert to integers to ensure counts are valid
    asv_table_full_filtered = asv_table_full_filtered.astype(int)
    #asv_table_full_filtered.to_csv(os.path.join(args.out, 'asv_table_full_filtered.csv'), index=True)

    # Determine the minimum sequencing depth (excluding 0)
    min_depth = asv_table_full_filtered.sum(axis=0).min()
    print(f"rarefy to {min_depth} reads:")

    # Perform rarefaction per run (column)
    rarefied_asv_table = asv_table_full_filtered.apply(lambda col: pd.Series(subsample_counts(col.values, n=min_depth), index=col.index), axis=0)
    #rarefied_asv_table.to_csv(os.path.join(args.out, 'rarefied_asv_table.csv'), index=True)

    # 6. Merge ASV tax tables
    tax_tables = []
    for run in run_list:
        path = os.path.join(base_dir, "runs", run, "dada2", "ASV_tax_species.silva_138_2.tsv.gz")
        df_tax = pd.read_csv(path, sep="\t")
        tax_tables.append(df_tax)

    asv_tax_full = pd.concat(tax_tables).drop_duplicates(subset="ASV_ID")
    asv_tax_full.set_index("ASV_ID", inplace=True)

    # Ensure tax and count table match
    #asv_tax_full.to_csv(os.path.join(args.out, 'asv_tax_full.csv'), index=True)
    common_index = asv_table_full_filtered.index.intersection(asv_tax_full.index)
    asv_table_full_filtered = asv_table_full_filtered.loc[common_index]
    asv_tax_full = asv_tax_full.loc[common_index]

    # 7. Aggregate counts at a given taxonomic rank (e.g., "Genus")
    merged = rarefied_asv_table.join(asv_tax_full[[rank]])
    agg_table = merged.groupby(rank).sum()

    #agg_table.to_csv('agg_table.csv', index=False)


    # 8. Calculate Bray-Curtis and Ruzicka beta diversity
    bray_curtis_matrix = pd.DataFrame(
        squareform(pdist(agg_table.T, metric='braycurtis')),
        index=agg_table.columns,
        columns=agg_table.columns
    )

    ruzicka_matrix = pd.DataFrame(
        squareform(pdist(agg_table.T, metric=ruzicka)),
        index=agg_table.columns,
        columns=agg_table.columns
    )

    betadiv_matrix = ruzicka_matrix

    # 9. Plot heatmap
    # Assume betadiv_matrix is a square DataFrame with ASV_IDs as index and columns
    # Step 1: Convert to condensed form for clustering
    condensed_dist = ssd.squareform(betadiv_matrix.values)

    # Step 2: Compute linkage
    linkage = sch.linkage(condensed_dist, method='average')  # or other method

    # Step 3: Get the order of rows/columns from dendrogram
    dendro = sch.dendrogram(linkage, no_plot=True)
    ordered_indices = dendro['leaves']

    # Step 4: Get ASV_IDs in the new order
    ordered_labels = betadiv_matrix.index[ordered_indices]

    # Step 5: Reorder the matrix using label-based indexing
    ordered_matrix = betadiv_matrix.loc[ordered_labels, ordered_labels]

    # Step 6: Plot the heatmap and save to PDF
    plt.figure(figsize=(12, 12))
    sns.heatmap(ordered_matrix, cmap="viridis", xticklabels=False, yticklabels=False)
    plt.title(f"Beta Diversity Heatmap (Clustered) at {rank}")
    plt.tight_layout()
    plt.savefig(os.path.join(args.out, f"betadiv_heatmap_clustered_{args.dataset_name}_{rank}.jpg"), dpi=300, bbox_inches='tight')


    # 10. Find the most central cluster
    # Step 1: Find the "central" run — lowest average distance to all others
    avg_distances = betadiv_matrix.mean(axis=0)
    central_run = avg_distances.idxmin()
    print(f"Central run (most consensus-like): {central_run}")

    # Step 2: Cut the dendrogram into clusters (you can try other t values too)
    from scipy.cluster.hierarchy import fcluster
    # Try a distance threshold; tweak if necessary
    cluster_assignments = fcluster(linkage, t=0.02, criterion='distance')  

    # Step 3: Score each cluster
    cluster_scores = {}
    unique_clusters = np.unique(cluster_assignments)
    labels = betadiv_matrix.index.to_numpy()

    for c in unique_clusters:
        member_indices = np.where(cluster_assignments == c)[0]
        member_labels = labels[member_indices]
        
        # Intra-cluster distance
        intra_dist = betadiv_matrix.loc[member_labels, member_labels]
        mean_intra = intra_dist.values[np.triu_indices_from(intra_dist.values, k=1)].mean()
        
        # Mean distance to central run
        to_central = betadiv_matrix.loc[member_labels, central_run].mean()
        
        # Record score
        cluster_scores[c] = {
            'size': len(member_labels),
            'mean_to_central': to_central,
            'intra_cluster_distance': mean_intra,
            'members': list(member_labels)
        }

    # Step 4: Filter and rank clusters (e.g., only consider clusters with ≥3 members)
    min_cluster_size = 3
    filtered_clusters = {k: v for k, v in cluster_scores.items() if v['size'] >= min_cluster_size}

    # Rank by proximity to central run, then tightness
    best_cluster = min(filtered_clusters.items(), key=lambda x: (x[1]['mean_to_central'], x[1]['intra_cluster_distance']))
    best_cluster_id, best_cluster_info = best_cluster

    print(f"Best cluster ID: {best_cluster_id}")
    print(f"Cluster members: {best_cluster_info['members']}")
    print(f"Size: {best_cluster_info['size']}")
    print(f"Mean distance to central run: {best_cluster_info['mean_to_central']:.4f}")
    print(f"Mean intra-cluster distance: {best_cluster_info['intra_cluster_distance']:.4f}")

    # Step 5 (optional): Export best cluster members
    pd.Series(best_cluster_info['members']).to_csv(os.path.join(args.out, 'best_cluster_runs.csv'), index=False, header=False)

if __name__ == "__main__":
    sys.exit(main())


