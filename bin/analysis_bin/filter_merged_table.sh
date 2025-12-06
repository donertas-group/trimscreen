#!/usr/bin/env bash
set -euo pipefail

# Base paths
outdir="/scratch/shire/ssd/pipeline/16s_nf_pipeline/analysis_mock/output"

#d="callahan";D="mock12"
#d="schirmer2015"; D="schirmer2015"
#d="hc227";D="hc227_v3v4"
d="kozich";D="mock13-15"


# Array of numeric suffixes — may be empty
numbers=() #1 2 4 5 6 7 8 9 11 12 17 18 19 20)

# If the array is empty, we run once with an empty suffix
if [ ${#numbers[@]} -eq 0 ]; then
    numbers=("")
fi

for n in "${numbers[@]}"; do
    # Build suffix
    if [[ -z "$n" ]]; then
        s=""
    else
        s=".$n"
    fi

    input_file="${outdir}/merged_table.${D}${s}.csv"
    ref_file="/scratch/shire/ssd/pipeline/16s_nf_pipeline/${D}/output${s}/compare_runs/filtered_table.csv"
    output_file="${outdir}/merged_table_filtered.${D}${s}.csv"

    echo "Processing suffix '$s' ..."
    echo "  Input:  $input_file"
    echo "  Ref:    $ref_file"
    echo "  Output: $output_file"

    # Extract the set of allowed run values from the reference file
    # Then filter the input file based on run_id
    awk -F',' '
        NR==FNR {
            if (FNR==1) {
                # Find the column index for the "run" field in reference file
                for (i=1; i<=NF; i++) if ($i=="run") run_col=i
                next
            }
            allowed[$run_col] = 1
            next
        }
        NR!=FNR {
            if (FNR==1) {
                # Find the column index for "run_id" in the input file
                for (i=1; i<=NF; i++) if ($i=="run_id") id_col=i
                print
                next
            }
            if ($id_col in allowed) print
        }
    ' "$ref_file" "$input_file" > "$output_file"

    echo "Finished: $d"
    echo
done

