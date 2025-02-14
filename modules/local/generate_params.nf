process GENERATE_PARAMS {
    input:
    str truncrangef
    str truncranger
    val outdir
 
    output:
    path "${outdir}/generate_params/summary_params_settings.csv"

    script:
    """
    generate_params.py --trunclenf_range $trunclenf_range --trunclenr_range $trunclenr_range --outdir $outdir
    """
}

