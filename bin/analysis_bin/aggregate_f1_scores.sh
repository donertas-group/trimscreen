#!/usr/bin/env bash

set -euo pipefail

# Input files
SELECTED_RUNS="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output/selected_runs.tsv"
METADATA="/scratch/shire/data/nj/projects/trimscreen_manuscript/metadata/table_1_metadata.csv"
MERGED_DIR="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output"

# Output file
OUT="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output/f1_scores.tsv"

# Write header
echo -e "dataset\trun\tf1" > "$OUT"

# Read selected_runs.tsv row by row (skip header)
tail -n +2 "$SELECTED_RUNS" | while IFS=$'\t' read -r dataset selected_run rest; do

    # Get dataset_name_in_pipeline from metadata
    dataset_name_in_pipeline=$(awk -F',' -v ds="$dataset" '
        NR==1 {
            for (i=1; i<=NF; i++) {
                if ($i=="dataset_label") dsl=i
                if ($i=="dataset_name_in_pipeline") dnip=i
            }
            next
        }
        $dsl==ds { print $dnip; exit }
    ' "$METADATA")

    # Skip if not found
    [[ -z "$dataset_name_in_pipeline" ]] && continue

    merged_file="$MERGED_DIR/merged_table.${dataset_name_in_pipeline}.csv"

    # Skip if merged table does not exist
    [[ ! -f "$merged_file" ]] && continue

    # Get f1_score for selected_run
    f1_score=$(awk -F',' -v run="$selected_run" '
        NR==1 {
            for (i=1; i<=NF; i++) {
                if ($i=="run_id") rid=i
                if ($i=="f1_score") f1=i
            }
            next
        }
        $rid==run { print $f1; exit }
    ' "$merged_file")

    # Write result if f1_score found
    if [[ -n "$f1_score" ]]; then
        echo -e "${dataset}\t${selected_run}\t${f1_score}" >> "$OUT"
    fi

done

