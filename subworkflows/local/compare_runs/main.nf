include { SUMMARISE_RUN                    } from '../../../modules/local/summarise_run'
include { MERGE_RUN_SUMMARIES              } from '../../../modules/local/merge_run_summaries'
include { FILTER_RUNS                      } from '../../../modules/local/filter_runs'
include { FIND_BEST_RUN                    } from '../../../modules/local/find_best_run'

workflow COMPARE_RUNS {
    take:
    runs_summary 
    runs_asv_table
    runs_asv_tax

    main:
    runs_summary
        .combine(runs_asv_table, by:0)
        .combine(runs_asv_tax, by:0)
        .set{ch_run_data}

    SUMMARISE_RUN (ch_run_data)
    ch_run_table = SUMMARISE_RUN.out.csv
    
    MERGE_RUN_SUMMARIES (
        ch_run_table.map { it[1] }.collect()
    ) 

    full_table = MERGE_RUN_SUMMARIES.out.csv
    
    FILTER_RUNS( full_table )
    filtered_table = FILTER_RUNS.out.csv 

    if (params.metadata) {
    ch_metadata = Channel.fromPath("${params.metadata}", checkIfExists: true)
    }

    FIND_BEST_RUN (filtered_table, ch_metadata)
    

    emit:
    best_runs = FIND_BEST_RUN.out

}
