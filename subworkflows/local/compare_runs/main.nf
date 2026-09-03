include { SUMMARISE_RUN as SUMMARISE_RAREFIED_RUN } from '../../../modules/local/summarise_run'
include { MERGE_SUMMARIES as MERGE_RAREFIED_RUN_SUMMARIES } from '../../../modules/local/merge_summaries'
include { MERGE_SUMMARIES as MERGE_RAREFIED_SAMPLERUN_SUMMARIES } from '../../../modules/local/merge_summaries'
include { FILTER_RUNS } from '../../../modules/local/filter_runs'
include { RAREFY_RUNS } from '../../../modules/local/rarefy_runs'
include { SCORE_BASED_OPTIMISATION } from '../../../modules/local/score_based_optimisation'

workflow COMPARE_RUNS {
    take:
    full_table
    ch_run_data
    runs_summary
    runs_asv_tax
    run_qtrim_string
    ch_metadata
    ch_samplesheet
 
    main:
    // If run comparison is to be based on rarefied runs, filter by reads retention first and then rarefy

    ch_all_samples = ch_samplesheet
        .map { meta, read1, read2 -> meta.id }
        .collect()
        .map { ids -> ids.join(',') }

    if (params.compare_rarefied_runs) {

        FILTER_RUNS( full_table, '.filtered', params.min_reads, ch_all_samples )

        // FILTER_RUNS stdout now contains three JSON lines:
        //   1) good runs
        //   2) "chronic" low-read samples (never clear min_reads in ANY
        //      run) - error-level, stops the pipeline if it wipes out all
        //      runs
        //   3) samples missing from the table entirely, vs. the
        //      samplesheet - warning-level only, pipeline continues
        // See filter_runs.py for details.
        ch_filter_result = FILTER_RUNS.out.filtered
            .map { stdout, csv ->
                def lines = stdout.trim().readLines()
                def good_runs = new groovy.json.JsonSlurper().parseText(lines[0])
                def chronic_low_samples = lines.size() > 1
                    ? new groovy.json.JsonSlurper().parseText(lines[1])
                    : []
                def missing_from_table_samples = lines.size() > 2
                    ? new groovy.json.JsonSlurper().parseText(lines[2])
                    : []
                [good_runs, chronic_low_samples, missing_from_table_samples]
            }

        // Safety check: count items and validate
        ch_filter_result
            .map { good_runs, problem_samples, missing_from_table_samples ->

                if (missing_from_table_samples.size() > 0) {
                    def sample_str = missing_from_table_samples.size() == 1
                        ? "sample ${missing_from_table_samples[0]} is"
                        : "samples ${missing_from_table_samples[0..-2].join(', ')} and ${missing_from_table_samples[-1]} are"
                    log.warn "${sample_str} missing from all the runs. Please check ASV length filtering parameters."
                }

                def count = good_runs.size()
                def run_ids = good_runs.collect { id, depth -> id }

                if (count == 0 && problem_samples.size() > 0) {
                    // More specific diagnosis: chronic low-read samples
                    // present in every run wiped out all runs, which would
                    // otherwise just look like an unexplained "0 runs left".
                    def sample_str = problem_samples.size() == 1
                        ? "sample ${problem_samples[0]}"
                        : "samples ${problem_samples[0..-2].join(', ')} and ${problem_samples[-1]}"
                    def verb = problem_samples.size() == 1 ? 'has' : 'have'
                    error "${sample_str} ${verb} nreads < ${params.min_reads} among all runs. No run is left after filtering because of this. Consider removing poorly sequenced samples from samplesheet, or lower `--min_reads` parameter."
                } else if (count == 0) {
                    error "No run is left after filtering, use `--trunclenf_range`, `--trunclenr_range` or decrease `--screen_interval` to screen more runs"
                } else if (count == 1) {
                    error "Only one run (${run_ids[0]}) is left after filtering, use `--trunclenf_range`, `--trunclenr_range` or decrease `--screen_interval` to screen more runs"
                } else if (count == 2) {
                    log.warn "Only two runs (${run_ids.join(', ')}) are left after filtering. Continue. However, consider using `--trunclenf_range`, `--trunclenr_range` or decrease `--screen_interval` to screen more runs"
                } else {
                    log.warn "${count} runs left after filtering. Continue."
                }
                return good_runs
            }
            .flatMap { v -> v }
            .set { ch_runs_filtered_valid }

        logged = false
        ch_runs_filtered_valid.subscribe { value ->
            if (!logged) {
                def reads = value[1]
                log.info "Rarefy filtered runs to ${reads} reads per sample."
                logged = true
            }
        }

        ch_runs_to_rarefy = ch_run_data
            .map { meta, stats, asv, tax -> [meta.run, meta, asv] }
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
            ch_rarefied_run_table.map { meta, samplerun_csv, run_csv -> samplerun_csv }.collect(), 
            'samplerun_summaries.filtered.rarefied'
        )

        MERGE_RAREFIED_RUN_SUMMARIES(
            ch_rarefied_run_table.map { meta, samplerun_csv, run_csv -> run_csv }.collect(), 
            'run_summaries.filtered.rarefied'
        )

        rarefied_table = MERGE_RAREFIED_SAMPLERUN_SUMMARIES.out.csv

        SCORE_BASED_OPTIMISATION (rarefied_table, ch_metadata, run_qtrim_string)
   
    } else {
        SCORE_BASED_OPTIMISATION (full_table, ch_metadata, run_qtrim_string)
    }

    emit:
    best_runs = SCORE_BASED_OPTIMISATION.out
}

