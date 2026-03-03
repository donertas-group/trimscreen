include { SUMMARISE_RUN                                              } from '../../../modules/local/summarise_run'
include { SUMMARISE_RUN as SUMMARISE_RAREFIED_RUN                    } from '../../../modules/local/summarise_run'
include { MERGE_SUMMARIES as MERGE_RUN_SUMMARIES                     } from '../../../modules/local/merge_summaries'
include { MERGE_SUMMARIES as MERGE_SAMPLERUN_SUMMARIES               } from '../../../modules/local/merge_summaries'
include { MERGE_SUMMARIES as MERGE_RAREFIED_RUN_SUMMARIES            } from '../../../modules/local/merge_summaries'
include { MERGE_SUMMARIES as MERGE_RAREFIED_SAMPLERUN_SUMMARIES      } from '../../../modules/local/merge_summaries'
include { FILTER_RUNS                                                } from '../../../modules/local/filter_runs'
include { RAREFY_RUNS                                                } from '../../../modules/local/rarefy_runs'
include { SCORE_BASED_OPTIMISATION                                   } from '../../../modules/local/score_based_optimisation'

workflow COMPARE_RUNS {
    take:
    runs_summary 
    runs_asv_table
    runs_asv_tax

    main:
    ch_run_data = runs_summary
        .combine(runs_asv_table, by:0)
        .combine(runs_asv_tax, by:0)

    if (params.metadata) {
        ch_metadata = Channel.fromPath("${params.metadata}", checkIfExists: true)
    } else {
        ch_metadata = Channel.value(file('NO_FILE'))  // Creates a single-item channel
    }

    SUMMARISE_RUN (ch_run_data, ch_metadata)
    
    ch_run_table = SUMMARISE_RUN.out.csv
    
    MERGE_SAMPLERUN_SUMMARIES (
        ch_run_table.map { meta, samplerun_csv, run_csv -> tuple(meta, samplerun_csv) }.collect(), 
        'samplerun_summaries' 
    ) 

    MERGE_RUN_SUMMARIES (
        ch_run_table.map { meta, samplerun_csv, run_csv -> tuple(meta, run_csv) }.collect(), 
        'run_summaries' 
    ) 

    full_table = MERGE_SAMPLERUN_SUMMARIES.out.csv


    //  
    // If run comparision is to be based on rarefied runs, filter by reads retention first and then rarefy
    //
    if (params.compare_rarefied_runs) {

        FILTER_RUNS( full_table, '.filtered' )

        ch_runs_filtered = FILTER_RUNS.out.filtered
            .map {stdout, csv -> stdout}
            .splitJson()
            .map { [it[0], it[1]] }


        // Safety check: count items and validate
        ch_runs_filtered
            .toList()
            .map { runs_list ->
                def count = runs_list.size()
                def run_ids = runs_list.collect { id, depth -> id }
                
                if (count == 0) {
                    error "No run is left after filtering, use `--trunclenf_range`, `--trunclenr_range` or decrease `--screen_interval` to screen more runs"
                } else if (count == 1) {
                    error "Only one run (${run_ids[0]}) is left after filtering, use `--trunclenf_range`, `--trunclenr_range` or decrease `--screen_interval` to screen more runs"
                } else if (count == 2) {
                    log.warn "Only two runs (${run_ids.join(', ')}) are left after filtering, consider using `--trunclenf_range`, `--trunclenr_range` or decrease `--screen_interval` to screen more runs"
                } else {
                    log.info "${count} runs left after filtering"
                }
                
                return runs_list
            }
            .flatMap { v -> v }
            .set { ch_runs_filtered_valid }


        logged = false
        ch_runs_filtered_valid.subscribe { value ->
            // only log for the first element
            if (!logged) {
                def reads = value[1]
                log.info "Rarefy filtered runs to ${reads} reads per sample."
                logged = true
            }
        }

        ch_runs_to_rarefy = ch_run_data
            .map { meta, stats, asv, tax -> [meta.runID, meta, asv] }
            .join ( ch_runs_filtered_valid )
            .map { run, meta, asv, depth -> [meta, asv, depth]}

       
        RAREFY_RUNS(ch_runs_to_rarefy)
     
        asv_table_rarefied = RAREFY_RUNS.out.tsv

        ch_rarefied_run_data = runs_summary
            .combine(asv_table_rarefied, by:0)
            .combine(runs_asv_tax, by:0)

        SUMMARISE_RAREFIED_RUN( ch_rarefied_run_data, ch_metadata ) 
        ch_rarefied_run_table = SUMMARISE_RAREFIED_RUN.out.csv


        MERGE_RAREFIED_SAMPLERUN_SUMMARIES(
            ch_rarefied_run_table.map { meta, samplerun_csv, run_csv -> tuple(meta, samplerun_csv) }.collect(), 
            'samplerun_summaries.filtered.rarefied')

        MERGE_RAREFIED_RUN_SUMMARIES(
            ch_rarefied_run_table.map { meta, samplerun_csv, run_csv -> tuple(meta, run_csv) }.collect(), 
            'run_summaries.filtered.rarefied')

        rarefied_table = MERGE_RAREFIED_SAMPLERUN_SUMMARIES.out.csv

        SCORE_BASED_OPTIMISATION (rarefied_table, ch_metadata)
   
    } else {

        SCORE_BASED_OPTIMISATION (full_table, ch_metadata)

    }

    emit:
    best_runs = SCORE_BASED_OPTIMISATION.out

}
