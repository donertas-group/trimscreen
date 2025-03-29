include { SUMMARISE_RUN        } from '../../../modules/local/summarise_run'
include { MERGE_RUN_SUMMARIES        } from '../../../modules/local/merge_run_summaries'
//include { COMPARE_RUNS_FILTER       } from '../../../modules/local/compare_runs_filter'
//include { COMPARE_RUNS_DECIDE       } from '../../../modules/local/compare_runs_decide'

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

  //  COMPARE_RUNS_FILTER ()
    //COMPARE_RUNS_DECIDE ()

   // emit:

}
