#!/usr/bin/env bash

dataset_name="mock13-15"

    ./driver_compare_w_mock.allruns.py -D "${dataset_name}" --true "/scratch/shire/data/nj/raw_data/published/mockrobiota/mock13/true_composition_mock_13_14_15.csv"
    ./betadiv.py -D "${dataset_name}" 
    ./driver_reps_similarity.allruns.py \
        -D "${dataset_name}" \
        --f1_file "f1_scores_${dataset_name}_Genus.txt" \
        --median_distance_file "median_distances_per_sample.${dataset_name}.csv"



