#!/usr/bin/env bash

numbers=(1 2 5 6)
dataset_name="bokulich2015"

for n in "${numbers[@]}"; do
    s=".${n}"
    ./driver_compare_w_true.allruns.py -D "${dataset_name}" --out_suffix "${s}" --true /scratch/shire/data/nj/raw_data/published/mockrobiota/mock03/true_composition.csv
    ./driver_merge_stats.py \
        -D "${dataset_name}" \
        --out_suffix "${s}" \
        --f1_file "f1_scores_${dataset_name}${s}_Genus.csv"

done



