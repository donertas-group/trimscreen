#!/usr/bin/env bash

numbers=(1 2 3 4 5 6 7 8 9 10 11 12 13 18 19 20 21)
dataset_name="schirmer2015"

for n in "${numbers[@]}"; do
    s=".${n}"
    ./driver_compare_w_true.allruns.py -D "${dataset_name}" --out_suffix "${s}"
    ./driver_merge_stats.py \
        -D "${dataset_name}" \
        --out_suffix "${s}" \
        --f1_file "f1_scores_${dataset_name}${s}_Genus.csv" 

done



