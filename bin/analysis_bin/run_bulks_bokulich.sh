#!/usr/bin/env bash

numbers=(1 2 5 6)
dataset_name="bokulich2015"

for n in "${numbers[@]}"; do
    s=".${n}"
    ./driver_compare_w_mock.allruns.py -D "${dataset_name}" --out_suffix "${s}" --true /scratch/shire/data/nj/raw_data/published/mockrobiota/mock03/true_composition.csv
    ./betadiv.py -D "${dataset_name}" --out_suffix "${s}"
    ./driver_reps_similarity.allruns.py \
        -D "${dataset_name}" \
        --out_suffix "${s}" \
        --f1_file "f1_scores_${dataset_name}${s}_Genus.txt" \
        --median_distance_file "median_distances_per_sample.${dataset_name}${s}.csv"

done



