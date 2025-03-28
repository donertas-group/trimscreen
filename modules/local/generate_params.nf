process GENERATE_PARAMS {

    input:
    val trunclenf_range
    val trunclenr_range

    output:
    path "summary_params_settings.csv", emit: params_csv

    script:
    """ 
    generate_params.py -f $trunclenf_range -r $trunclenr_range -o .
    """
}

