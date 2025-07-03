#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from skbio.stats import subsample_counts  # skbio is efficient for rarefaction
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch


# 1. Set sample ID
sampleID = "mock_21"
base_dir = "../mock21/output"
compare_runs_csv = os.path.join(base_dir, "compare_runs", "filtered_table.csv")

# 2. Get run list
df_runs = pd.read_csv(compare_runs_csv)
run_list = df_runs['run'].dropna().unique().tolist()

# 3. Read ASV tables and merge
asv_tables = []
for run in run_list:
    path = os.path.join(base_dir, "runs", run, "dada2", "ASV_table.tsv")
    df = pd.read_csv(path, sep="\t", usecols=["ASV_ID", sampleID])
    df = df.rename(columns={sampleID: run})
    asv_tables.append(df)

# 4. Full outer join on ASV_ID
asv_table_full = asv_tables[0]
for df in asv_tables[1:]:
    asv_table_full = pd.merge(asv_table_full, df, on="ASV_ID", how="outer")

asv_table_full = asv_table_full.fillna(0)
asv_table_full.set_index("ASV_ID", inplace=True)

# 4.5 Rarefy each run to the same sequencing depth (per column)
# Convert to integers to ensure counts are valid
asv_table_full = asv_table_full.astype(int)
#asv_table_full.to_csv('./output/asv_table_full.csv', index=False)

# Determine the minimum sequencing depth (excluding 0)
min_depth = asv_table_full.sum(axis=0).min()
print("rarefy to N of reads:")
print(min_depth)
# Perform rarefaction per run (column)
rarefied_table = asv_table_full.apply(lambda col: pd.Series(subsample_counts(col.values, n=min_depth), index=col.index), axis=0)
#rarefied_table.to_csv('./output/rarefied_table.csv', index=False)


# 5. Merge ASV tax tables
tax_tables = []
for run in run_list:
    path = os.path.join(base_dir, "runs", run, "dada2", "ASV_tax_species.silva_138.tsv")
    df_tax = pd.read_csv(path, sep="\t")
    tax_tables.append(df_tax)

asv_tax_full = pd.concat(tax_tables).drop_duplicates(subset="ASV_ID")
asv_tax_full.set_index("ASV_ID", inplace=True)

# Ensure tax and count table match
asv_tax_full = asv_tax_full.loc[asv_table_full.index]
#asv_tax_full.to_csv('./output/asv_tax_full.csv', index=False)

# 6. Aggregate counts at a given taxonomic rank (e.g., "Genus")
rank = "Genus"
rank = "Family"
merged = rarefied_table.join(asv_tax_full[[rank]])
agg_table = merged.groupby(rank).sum()

#agg_table.to_csv('agg_table.csv', index=False)


# 7. Calculate Bray-Curtis beta diversity
# Transpose so rows = samples, columns = features
bray_curtis_matrix = pd.DataFrame(
    squareform(pdist(agg_table.T, metric='braycurtis')),
    index=agg_table.columns,
    columns=agg_table.columns
)

#bray_curtis_matrix.to_csv('bray_curtis_matrix.csv', index=True)

# 8. Plot heatmap
# Assume bray_curtis_matrix is a square DataFrame with ASV_IDs as index and columns
# Step 1: Convert to condensed form for clustering
condensed_dist = ssd.squareform(bray_curtis_matrix.values)

# Step 2: Compute linkage
linkage = sch.linkage(condensed_dist, method='average')  # or other method

# Step 3: Get the order of rows/columns from dendrogram
dendro = sch.dendrogram(linkage, no_plot=True)
ordered_indices = dendro['leaves']

# Step 4: Get ASV_IDs in the new order
ordered_labels = bray_curtis_matrix.index[ordered_indices]

# Step 5: Reorder the matrix using label-based indexing
ordered_matrix = bray_curtis_matrix.loc[ordered_labels, ordered_labels]

# Step 6: Plot the heatmap and save to PDF
plt.figure(figsize=(12, 12))
sns.heatmap(ordered_matrix, cmap="viridis", xticklabels=False, yticklabels=False)
plt.title(f"Bray-Curtis Beta Diversity Heatmap (Clustered) at {rank}")
plt.tight_layout()
plt.savefig(f"./output/bray_curtis_heatmap_clustered_{rank}.jpg", dpi=300, bbox_inches='tight')

