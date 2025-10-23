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

    # 1. Set dataset directories
    base_dir = f"/scratch/shire/ssd/pipeline/16s_nf_pipeline/{args.dataset_name}/output"
    compare_runs_csv = os.path.join(base_dir, "compare_runs", "filtered_table.csv")

    # 2. Get run list
    df_runs = pd.read_csv(compare_runs_csv)
    run_list = df_runs['run'].dropna().unique().tolist()

    # 3. Read ASV tables for each run
    asv_tables = []
    for run in run_list:
        path = os.path.join(base_dir, "runs", run, "dada2", "ASV_table.tsv.gz")
        df = pd.read_csv(path, sep="\t")

        # Identify all sample columns (everything except ASV_ID)
        sample_cols = [col for col in df.columns if col != "ASV_ID"]

        # Rename each sample column to include the run prefix
        df = df.rename(columns={col: f"{run}_{col}" for col in sample_cols})

        asv_tables.append(df)

    # 4. Full outer join on ASV_ID
    asv_table_full = asv_tables[0]
    for df in asv_tables[1:]:
        asv_table_full = pd.merge(asv_table_full, df, on="ASV_ID", how="outer")

    asv_table_full = asv_table_full.fillna(0)
    asv_table_full.set_index("ASV_ID", inplace=True)
    asv_table_full.to_csv(os.path.join(args.out, f'asv_table_full.{args.dataset_name}.csv'), index=True)
    
    # 5. Filter low abundance ASVs
    asv_table_full_filtered = asv_table_full[asv_table_full.max(axis=1) >= 20]

    # Convert to integers
    asv_table_full_filtered = asv_table_full_filtered.astype(int)
    asv_table_full_filtered.to_csv(os.path.join(args.out, f'asv_table_full_filtered.{args.dataset_name}.csv'), index=True)

    # 6. Rarefy each sample to the same depth
    min_depth = asv_table_full_filtered.sum(axis=0).min()
    print(f"Rarefying to {min_depth} reads per sample...")

    rarefied_asv_table = asv_table_full_filtered.apply(
        lambda col: pd.Series(subsample_counts(col.values, n=min_depth), index=col.index),
        axis=0
    )

    # 7. Merge taxonomic information
    tax_tables = []
    for run in run_list:
        path = os.path.join(base_dir, "runs", run, "dada2", "ASV_tax_species.silva_138_2.tsv.gz")
        df_tax = pd.read_csv(path, sep="\t")
        tax_tables.append(df_tax)

    asv_tax_full = pd.concat(tax_tables).drop_duplicates(subset="ASV_ID")
    asv_tax_full.set_index("ASV_ID", inplace=True)

    # Align taxa and abundance tables
    common_index = rarefied_asv_table.index.intersection(asv_tax_full.index)
    rarefied_asv_table = rarefied_asv_table.loc[common_index]
    asv_tax_full = asv_tax_full.loc[common_index]

    # 8. Aggregate counts by rank
    merged = rarefied_asv_table.join(asv_tax_full[[rank]])
    agg_table = merged.groupby(rank).sum()

    # Collapse sample columns by run
    agg_table_by_run = agg_table.groupby(
        agg_table.columns.str.extract(r'(^[^_]+)')[0],  # extract run prefix
        axis=1
    ).sum()

    # 9. compute Bray-Curtis between runs
    bray_curtis_matrix = pd.DataFrame(
        squareform(pdist(agg_table_by_run.T, metric='braycurtis')),
        index=agg_table_by_run.columns,
        columns=agg_table_by_run.columns
    )

    # bray_curtis_matrix.to_csv(os.path.join(args.out, 'bray_curtis_matrix.csv'))
        
    ruzicka_matrix = pd.DataFrame(
        squareform(pdist(agg_table_by_run.T, metric=ruzicka)),
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


