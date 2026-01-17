#!/usr/bin/env bash

dataset_name="mock13-15"

    ./driver_compare_w_true.allruns.py -D "${dataset_name}" --true "/scratch/shire/data/nj/raw_data/published/mockrobiota/mock13/true_composition.csv"
    ./driver_merge_stats.py \
        -D "${dataset_name}" \
        --f1_file "f1_scores_${dataset_name}_Genus.csv" 



