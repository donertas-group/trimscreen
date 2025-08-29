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
    """
    # Extract run name like run_123456 from file
    run_id="${best_run}"

    # Split into X and Y
    X=\${run_id:4:3}
    Y=\${run_id:7:3}

    # Write script
    cat <<- EOF > run_ampliseq.sh
    #!/usr/bin/bash
    nextflow run nf-core/ampliseq --FW_primer \$X --RV_primer \$Y
    EOF

    chmod +x run_ampliseq.sh
    """
}
