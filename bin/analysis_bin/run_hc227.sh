#!/usr/bin/env bash

dataset_name="hc227_v3v4"

    ./driver_compare_w_true.allruns.py -D "${dataset_name}" 
    ./driver_merge_stats.py \
        -D "${dataset_name}" \
        --f1_file "f1_scores_${dataset_name}_Genus.csv" 



