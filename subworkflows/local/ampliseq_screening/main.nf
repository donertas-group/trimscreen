include { AMPLISEQ_SIMPLIFIED                                       } from '../ampliseq_simplified/main'


workflow AMPLISEQ_SCREENING {
    take:
    ch_samplesheet
    
    main:
    // generate sets of parameters based on input ranges
    GENERATE_PARAMS(params.trunclenf_range, params.trunclenr_range)

    // create a channel with parameters as input to ampliseq (simplified from nf-core)
    ch_params = GENERATE_PARAMS.out.params_csv
    .splitCsv(header: true, sep: ',')
    .map { row -> tuple(row.runID, row.trunclenf, row.trunclenr) }

    AMPLISEQ_SIMPLIFIED(ch_samplesheet, ch_params)
}

/* Subworkflow
workflow AMPLISEQ_SIMPLE {
    take:
    ch_params

    main:
    PROCESS_1(ch_params) | view()
    // PROCESS_2(PROCESS_1.out)

    //emit:
    // Define your outputs here
}

process PROCESS_1 {
    input:
    tuple val(runID), val(trunclenf), val(trunclenr) // ... more parameters

    output:
    stdout

    script:
    """
    echo "$runID $trunclenf $trunclenr"
    """
}
*/
process GENERATE_PARAMS {

    def out_path = file(params.outdir).toString() + '/generate_params/'

    publishDir "$out_path", mode: 'copy'

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

