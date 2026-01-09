#!/usr/bin/env bash

numbers=(1 2 3 4 5 6 7 8 9 12 13)
dataset_name="schirmer2015"

for n in "${numbers[@]}"; do
    s=".${n}"
    ./driver_compare_w_mock.allruns.py -D "${dataset_name}" --out_suffix "${s}"
    ./betadiv.py -D "${dataset_name}" --out_suffix "${s}"
    ./driver_reps_similarity.allruns.py \
        -D "${dataset_name}" \
        --out_suffix "${s}" \
        --f1_file "f1_scores_${dataset_name}${s}_Genus.txt" \
        --median_distance_file "median_distances_per_sample.${dataset_name}${s}.csv"

done



