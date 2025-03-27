include { AMPLISEQ_SIMPLIFIED                                       } from '../ampliseq_simplified/main'
include { COMPARE_RUNS                                              } from '../compare_runs/main'

include { GENERATE_PARAMS                                           } from '../../../modules/local/generate_params'

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


    ch_samplesheet.combine(ch_params)
    .map { meta, read1, read2, _, runID, trunclenf, trunclenr -> 
    def new_meta = meta + [ sample: meta.id, id: "${meta.id}.${runID}", runID: runID, run: runID, trunclenf: trunclenf, trunclenr: trunclenr]  
    tuple(new_meta, read1, read2, _)}
    .set{ ch_samplesheet_w_params }
    

    AMPLISEQ_SIMPLIFIED(ch_samplesheet_w_params)
    ch_stats = AMPLISEQ_SIMPLIFIED.out.runs_summary
    ch_asv = AMPLISEQ_SIMPLIFIED.out.runs_asv_table
    ch_tax = AMPLISEQ_SIMPLIFIED.out.runs_asv_tax



    if (true){//!params.skip_run_comparison) {
        COMPARE_RUNS (
        ch_stats,
        ch_asv,
        ch_tax
        )  
    }



}
/*
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
}*/

