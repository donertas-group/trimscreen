process GENERATE_PARAMS {

    input:
      val marker_size_min
      val minimum_overlap
      val step_size
      val read_length
      val FW_primer_len
      val RV_primer_len
      val trunclenf_range
      val trunclenr_range
  
    output:
      path "summary_params_settings.csv", emit: params_csv

    script:
    def trunclen_arg = (trunclenf_range != "" && 
                       trunclenr_range != "") ?
                       "--trunclenf_range ${trunclenf_range} --trunclenr_range ${trunclenr_range}" : ""
    """
    generate_params.py --marker_size_min $marker_size_min \
        --FW_primer_len $FW_primer_len --RV_primer_len $RV_primer_len \
        --minimum_overlap $minimum_overlap \
        --step_size $step_size \
        ${trunclen_arg} \
        --read_length $read_length -o .
    """
}

