include { SUMMARISE_RUN                    } from '../../../modules/local/summarise_run'
include { MERGE_RUN_SUMMARIES              } from '../../../modules/local/merge_run_summaries'
include { FILTER_RUNS                      } from '../../../modules/local/filter_runs'
include { RAREFY_RUNS                      } from '../../../modules/local/rarefy_runs'
include { RANK_BASED_OPTIMISATION          } from '../../../modules/local/rank_based_optimisation'

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

    if (params.metadata) {
        ch_metadata = Channel.fromPath("${params.metadata}", checkIfExists: true)
    } else {
        ch_metadata = Channel.value(file('NO_FILE'))  // Creates a single-item channel
    }

    ch_run_data
        .combine(ch_metadata)
        .set { ch_run_data_with_metadata }

    SUMMARISE_RUN (ch_run_data, ch_metadata)
    
    ch_run_table = SUMMARISE_RUN.out.csv
    
    MERGE_RUN_SUMMARIES (
        ch_run_table.map { it[1] }.collect()
    ) 

    full_table = MERGE_RUN_SUMMARIES.out.csv

    FILTER_RUNS( full_table )

    ch_runs_filtered = FILTER_RUNS.out.filtered
        .map {stdout, csv -> stdout}
        .splitJson()
        .map { [it[0], it[1]] }

    logged = false
    ch_runs_filtered.subscribe { value ->
        // only log for the first element
        if (!logged) {
            def reads = value[1]
            log.info "Rarefy filtered runs to ${reads} reads per sample."
            logged = true
        }
    }

    ch_run_data
        .map { meta, stats, asv, tax -> [meta.runID, meta, asv] }
        .join ( ch_runs_filtered )
        .map { run, meta, asv, depth -> [meta, asv, depth]}
        .set{ ch_runs_to_rarefy }

   // ch_runs_to_rarefy.view()


   
    RAREFY_RUNS(ch_runs_to_rarefy)
 
    ASV_table_rarefied = RAREFY_RUNS.out.tsv

    filtered_table = FILTER_RUNS.out.filtered
    .map {stdout, csv -> csv} 
    RANK_BASED_OPTIMISATION (filtered_table, ch_metadata)

    emit:
    best_runs = RANK_BASED_OPTIMISATION.out

}
