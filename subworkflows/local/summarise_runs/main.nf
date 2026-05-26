include { SUMMARISE_RUN } from '../../../modules/local/summarise_run'
include { MERGE_SUMMARIES as MERGE_RUN_SUMMARIES } from '../../../modules/local/merge_summaries'
include { MERGE_SUMMARIES as MERGE_SAMPLERUN_SUMMARIES } from '../../../modules/local/merge_summaries'

workflow SUMMARISE_RUNS {
    take:
    ch_run_data // [ meta, summary, asv_table, asv_tax ]
    ch_metadata

    main:
    SUMMARISE_RUN (ch_run_data, ch_metadata)
    
    ch_run_table = SUMMARISE_RUN.out.csv
    
    MERGE_SAMPLERUN_SUMMARIES (
        ch_run_table.map { meta, samplerun_csv, run_csv -> samplerun_csv }.collect(), 
        'samplerun_summaries' 
    ) 

    MERGE_RUN_SUMMARIES (
        ch_run_table.map { meta, samplerun_csv, run_csv -> run_csv }.collect(), 
        'run_summaries' 
    ) 

    emit:
    full_table = MERGE_SAMPLERUN_SUMMARIES.out.csv
}
