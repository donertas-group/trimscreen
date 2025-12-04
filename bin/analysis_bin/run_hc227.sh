#!/usr/bin/env bash

dataset_name="hc227_v3v4"

    ./driver_compare_w_mock.allruns.py -D "${dataset_name}" 
    ./betadiv.py -D "${dataset_name}" 
    ./driver_reps_similarity.allruns.py \
        -D "${dataset_name}" \
        --f1_file "f1_scores_${dataset_name}_Genus.txt" \
        --median_distance_file "median_distances_per_sample.${dataset_name}.csv"



