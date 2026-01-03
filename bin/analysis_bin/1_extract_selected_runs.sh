#!/usr/bin/env bash

# Output table
OUT="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output/selected_runs.tsv"
echo -e "dataset\tselected_run" > "$OUT"

# EDIT THIS LIST: dataset_id <TAB> path_to_report.txt
while read -r DATASET REPORT; do
    if [[ ! -f "$REPORT" ]]; then
        echo "Warning: $REPORT not found, skipping $DATASET" >&2
        continue
    fi

    RUN=$(grep -m1 "Optimal run (highest total score):" "$REPORT" \
          | sed 's/.*: //')

    echo -e "${DATASET}\t${RUN}" >> "$OUT"

done <<EOF
hc227	/scratch/shire/ssd/pipeline/16s_nf_pipeline/hc227_v3v4/output/compare_runs/report.txt
bokulich_b3_even	/scratch/shire/ssd/pipeline/16s_nf_pipeline/bokulich2015/output.1/compare_runs/report.txt
kozich	/scratch/shire/ssd/pipeline/16s_nf_pipeline/mock13-15/output/compare_runs/report.txt
EOF

