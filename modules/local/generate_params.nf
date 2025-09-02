process GENERATE_PARAMS {
    label 'process_single'

    input:
      val marker_size_min
      val minimum_overlap
      val screen_interval
      val read_length
      val FW_primer_len
      val RV_primer_len
      val trunclenf_range
      val trunclenr_range
  
    output:
      path "summary_params_settings.csv", emit: params_csv

    script:
    def trunclen_arg = (trunclenf_range != '' && 
                       trunclenr_range != '') ?
                       "--trunclenf_range ${trunclenf_range} --trunclenr_range ${trunclenr_range}" : ""
    """
    set -euo pipefail
    set -C  # prevent overwriting files

    generate_params.py --marker_size_min $marker_size_min \
        --FW_primer_len $FW_primer_len --RV_primer_len $RV_primer_len \
        --minimum_overlap $minimum_overlap \
        --screen_interval $screen_interval \
        ${trunclen_arg} \
        --read_length $read_length -o .

    # Check if file has meaningful content (more than just a newline)
    file_size=\$(stat -c%s summary_params_settings.csv)

    if [ "\$file_size" -le 1 ]; then
        echo "ERROR: CSV file is empty or contains only a newline"
        exit 1
    fi
    """
}

