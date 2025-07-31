process GENERATE_PARAMS {

    input:
      val marker_size_min
      val minimum_overlap
      val step_size
      val read_length
      val FW_primer_len
      val RV_primer_len

    output:
      path "summary_params_settings.csv", emit: params_csv

    script:
    """
    generate_params.py --marker_size_min $marker_size_min \
        --FW_primer_len $FW_primer_len --RV_primer_len $RV_primer_len \
        --minimum_overlap $minimum_overlap \
        --step_size $step_size \
        --read_length $read_length -o .
    """
}

