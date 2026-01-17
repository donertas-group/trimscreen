#!/usr/bin/env bash

numbers=(1)
#numbers=(2 3 4 5 6 7 8 9)
dataset_name="tourlousse2022"

for n in "${numbers[@]}"; do
    s=".${n}"
    ./driver_compare_w_true.allruns.py -D "${dataset_name}" --out_suffix "${s}"
    ./driver_merge_stats.py \
        -D "${dataset_name}" \
        --out_suffix "${s}" \
        --f1_file "f1_scores_${dataset_name}${s}_Genus.csv" 

done



