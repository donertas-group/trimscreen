import pandas as pd
import os
import shutil

import config

# Use the variables defined in config.py
workdir = config.workdir
outdir = config.outdir

# Paths to files. Note the the metadata tsv file has the same format as required by nf-core/ampliseq
filtered_table_csv = f"{outdir}/filtered_table.csv"
metadata_tsv = f"{workdir}/input/metadata.tsv"

# Read the CSV file and metadata table into pandas DataFrames
df0 = pd.read_csv(filtered_table_csv)
metadata = pd.read_csv(metadata_tsv, sep='\t')  # Assuming tab-separated values

# Filter out the samples that are not marked as 'sample' in the metadata table
sample_ids = metadata[metadata['condition'] == 'sample']['ID'].values

# Filter df to keep only the rows with the matching sample IDs
df = df0[df0['sample'].isin(sample_ids)]


# Initialize dictionaries to track best runs for Genus and Phylum
best_genus_runs = {}
best_phylum_runs = {}

# Iterate over each sample to find the best runs based on Genus and Phylum
for sample in df['sample'].unique():
    sample_data = df[df['sample'] == sample]

    # Find the best runs for Genus (all runs with the highest Genus value)
    max_genus_value = sample_data['Genus'].max()
    best_genus_rows = sample_data[sample_data['Genus'] == max_genus_value]
    for _, row in best_genus_rows.iterrows():
        best_genus_run_id = row['run']
        best_genus_runs[best_genus_run_id] = best_genus_runs.get(best_genus_run_id, 0) + 1

    # Find the best runs for Phylum (all runs with the highest Phylum value)
    max_phylum_value = sample_data['Phylum'].max()
    best_phylum_rows = sample_data[sample_data['Phylum'] == max_phylum_value]
    for _, row in best_phylum_rows.iterrows():
        best_phylum_run_id = row['run']
        best_phylum_runs[best_phylum_run_id] = best_phylum_runs.get(best_phylum_run_id, 0) + 1

# Convert the dictionaries into DataFrames for easier ranking
genus_df = pd.DataFrame(list(best_genus_runs.items()), columns=['run', 'value'])
phylum_df = pd.DataFrame(list(best_phylum_runs.items()), columns=['run', 'value'])

# Rank the runs for Genus (with ties having the same rank)
genus_df['rank'] = genus_df['value'].rank(method='min', ascending=False).astype(int)  # 'min' ties the rank

# Rank the runs for Phylum (with ties having the same rank)
phylum_df['rank'] = phylum_df['value'].rank(method='min', ascending=False).astype(int)  # 'min' ties the rank

# Convert back to dictionaries
genus_ranks = dict(zip(genus_df['run'], genus_df['rank']))
phylum_ranks = dict(zip(phylum_df['run'], phylum_df['rank']))

# Combine ranks
combined_ranks = {}
for run in set(genus_ranks.keys()).union(phylum_ranks.keys()):
    genus_rank = genus_ranks.get(run, len(genus_ranks) + 1)  # Assign lowest rank for missing runs
    phylum_rank = phylum_ranks.get(run, len(phylum_ranks) + 1)  # Assign lowest rank for missing runs
    combined_ranks[run] = genus_rank + phylum_rank

# Determine the best run based on the highest sum rank value
best_run = min(combined_ranks, key=combined_ranks.get)

# Output the results
print("Runs with highest Genus:")
for run, rank in sorted(genus_ranks.items(), key=lambda item: item[1])[:15]:
    print(f"{run}, rank: {rank}, best for {best_genus_runs.get(run, 0)} samples")

print("\nRuns with highest Phylum:")
for run, rank in sorted(phylum_ranks.items(), key=lambda item: item[1])[:15]:
    print(f"{run}, rank: {rank}, best for {best_phylum_runs.get(run, 0)} samples")

# Find and print the best run(s) based on the combined rank
# Find the minimum total rank value
min_rank = min(combined_ranks.values())

# Identify all runs with the minimum rank
best_runs = [(run, rank) for run, rank in combined_ranks.items() if rank == min_rank]

# Print all the best runs
print("\nBest Runs:")
for run, rank in best_runs:
    trunclenf = df[df['run'] == run]['trunclenf'].values[0]
    trunclenr = df[df['run'] == run]['trunclenr'].values[0]
    print(f"Run: {run}, Combined Rank: {rank}, trunclenf: {trunclenf}, trunclenr: {trunclenr}")

# Move the directory of the best run to the output directory
#best_run_directory = os.path.join(outdir,"runs", best_run)

#if os.path.exists(best_run_directory):
#    destination_dir = os.path.join(outdir, "best_run")
#    shutil.copytree(best_run_directory, destination_dir, dirs_exist_ok=True)
#    print(f"\n{best_run} is the best run, run result output to {destination_dir}.")
#else:
#    print(f"\nRun directory {best_run_directory} does not exist.")
