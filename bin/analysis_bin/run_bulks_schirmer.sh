#!/usr/bin/env bash

numbers=(2 4 5 6 7 8 9 10 11 12 17 18 19 20)
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



