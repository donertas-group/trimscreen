process WRITE_AMPLISEQ_COMMAND {
    label 'process_single'

    publishDir "${params.outdir}",
        mode: "$params.publish_dir_mode",
        pattern: "*.sh",
        enabled: !params.skip_run_comparison

    input:
    val(best_run)

    output:
    path("run_ampliseq.sh")

    script:
    // Fail early if the mandatory ampliseq inputs are missing
    if( !params.input || !params.metadata )
        error "WRITE_AMPLISEQ_COMMAND: 'input' and 'metadata' must be set in params.yaml"
    if( !params.skip_cutadapt && !(params.FW_primer && params.RV_primer) )
        error "WRITE_AMPLISEQ_COMMAND: set 'FW_primer' and 'RV_primer', or 'skip_cutadapt: true', in params.yaml"

    // Primers are only written if cutadapt is not skipped
    def primer_args = params.skip_cutadapt
        ? "--skip_cutadapt"
        : "--FW_primer '${params.FW_primer}' --RV_primer '${params.RV_primer}'"

    """
    # Extract run name like run_123456 from file
    run_id="${best_run}"

    # Split into X and Y
    X=\${run_id:4:3}
    Y=\${run_id:7:3}

    # Write script
    cat <<- EOF > run_ampliseq.sh
    #!/usr/bin/bash
    nextflow run nf-core/ampliseq \\\\
        --input '${file(params.input).toUriString()}' \\\\
        --metadata '${file(params.metadata).toUriString()}' \\\\
        ${primer_args} \\\\
        --trunclenf \$X \\\\
        --trunclenr \$Y \\\\
        --outdir results_ampliseq
    EOF

    chmod +x run_ampliseq.sh
    """
}
