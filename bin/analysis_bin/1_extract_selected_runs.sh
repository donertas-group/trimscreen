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
callahan	/scratch/shire/ssd/pipeline/16s_nf_pipeline/mock12/output/compare_runs/report.txt
hc227	/scratch/shire/ssd/pipeline/16s_nf_pipeline/hc227_v3v4/output/compare_runs/report.txt
bokulich_b3_even	/scratch/shire/ssd/pipeline/16s_nf_pipeline/bokulich2015/output.1/compare_runs/report.txt
bokulich_b3_staggered	/scratch/shire/ssd/pipeline/16s_nf_pipeline/bokulich2015/output.2/compare_runs/report.txt
bokulich_b5_even	/scratch/shire/ssd/pipeline/16s_nf_pipeline/bokulich2015/output.5/compare_runs/report.txt
bokulich_b5_staggered	/scratch/shire/ssd/pipeline/16s_nf_pipeline/bokulich2015/output.6/compare_runs/report.txt
gohl_1	/scratch/shire/ssd/pipeline/16s_nf_pipeline/gohl2016/output.1/compare_runs/report.txt
gohl_2	/scratch/shire/ssd/pipeline/16s_nf_pipeline/gohl2016/output.2/compare_runs/report.txt
kozich	/scratch/shire/ssd/pipeline/16s_nf_pipeline/mock13-15/output/compare_runs/report.txt
schirmer_1	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.1/compare_runs/report.txt
schirmer_2	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.2/compare_runs/report.txt
schirmer_3	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.3/compare_runs/report.txt
schirmer_4	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.4/compare_runs/report.txt
schirmer_5	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.4/compare_runs/report.txt
schirmer_6	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.6/compare_runs/report.txt
schirmer_7	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.7/compare_runs/report.txt
schirmer_8	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.8/compare_runs/report.txt
schirmer_9	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.9/compare_runs/report.txt
schirmer_10	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.10/compare_runs/report.txt
schirmer_11	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.11/compare_runs/report.txt
schirmer_12	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.12/compare_runs/report.txt
schirmer_13	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.13/compare_runs/report.txt
schirmer_18	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.18/compare_runs/report.txt
schirmer_19	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.19/compare_runs/report.txt
schirmer_20	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.20/compare_runs/report.txt
schirmer_21	/scratch/shire/ssd/pipeline/16s_nf_pipeline/schirmer2015/output.21/compare_runs/report.txt
tourlousse2022.1	/scratch/shire/ssd/pipeline/16s_nf_pipeline/tourlousse2022/output.1/compare_runs/report.txt
EOF

