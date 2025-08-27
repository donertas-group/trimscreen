include { AMPLISEQ_SIMPLIFIED                                       } from '../ampliseq_simplified/main'
include { AMPLISEQ_SIMPLIFIED as AMPLISEQ_SIMPLIFIED_RERUN } from '../ampliseq_simplified/main'
include { GENERATE_PARAMS                                           } from '../../../modules/local/generate_params'
include { COMPARE_RUNS                  } from '../../../subworkflows/local/compare_runs/main'
include { CREATE_LINK                   } from '../../../modules/local/create_link'

workflow AMPLISEQ_SCREENING {
    take:
    ch_samplesheet
    
    main:
    // generate sets of parameters based on input ranges
    FW_primer_len = params.FW_primer ? params.FW_primer.size() : 0
    RV_primer_len = params.RV_primer ? params.RV_primer.size() : 0

    ch_read_length = ch_samplesheet
            .first()
            .map { meta, readfw, readrv -> 
                // Return the first FASTQ file (readfw)
                return readfw
            }
            .splitFastq(record: true, limit: 1)
            .map { record -> record.readString.length() }
        
    // Report extracted read length:
    ch_read_length.view { "Read length extracted from samplesheet: $it" }

    GENERATE_PARAMS (
        params.marker_size_min, 
        params.minimum_overlap, 
        params.step_size, 
        ch_read_length,
        FW_primer_len,
        RV_primer_len,
        params.trunclenf_range,
        params.trunclenr_range 
    )

    // Create a channel with parameters as input to ampliseq (simplified from nf-core)
    ch_params = GENERATE_PARAMS.out.params_csv
    .splitCsv(header: true, sep: ',')
    .map { row -> 
           def is_best_run = false // initiallise is_best_run to false
           return tuple(row.runID, row.trunclenf, row.trunclenr, is_best_run) } 


    // Run simplified version of ampliseq
    AMPLISEQ_SIMPLIFIED(ch_samplesheet, ch_params)



 
    //
    // SUBWORKFLOW: Compare runs (moved from AMPLISEQ_SIMPLIFIED)
    //
    if (!params.skip_run_comparison) {
        COMPARE_RUNS ( 
            AMPLISEQ_SIMPLIFIED.out.runs_summary, 
            AMPLISEQ_SIMPLIFIED.out.runs_asv_table, 
            AMPLISEQ_SIMPLIFIED.out.runs_asv_tax 
        )

        // Process comparison results to identify best run
        COMPARE_RUNS.out
            .map{ stdout, report -> stdout }
            .splitJson()
            .flatten()
            .set{ ch_best_run }

        // Filter outputs to get best run data
        AMPLISEQ_SIMPLIFIED.out.runs_asv_table
            .map { meta, file -> [meta.runID, meta, file] }
            .join( ch_best_run.map { runID -> [runID, true] }, by: 0)
            .map { runID, meta, file, _ -> [meta, file] }
            .set { ch_best_tsv }

        AMPLISEQ_SIMPLIFIED.out.runs_asv_fasta
            .map { meta, file -> [meta.runID, meta, file] }
            .join( ch_best_run.map { runID -> [runID, true] }, by: 0)
            .map { runID, meta, file, _ -> [meta, file] }
            .set { ch_best_fasta }


    } else {
        // If comparison is skipped, output all runs
        ch_best_tsv = AMPLISEQ_SIMPLIFIED.out.runs_asv_table
        ch_best_fasta = AMPLISEQ_SIMPLIFIED.out.runs_asv_fasta
    }



    if (!params.publish_all_runs) {
        // Create updated metadata channel with is_best_run = true for best runs
        ch_best_run_formatted = ch_best_run
        .map { runID -> [runID, true] }        

        ch_params_best = ch_params
            .map {runID, trunclenf, trunclenr, is_best_run -> [runID, trunclenf, trunclenr] }
            .join( ch_best_run_formatted)

        ch_params_best.view()

        // Re-run AMPLISEQ_SIMPLIFIED with updated metadata (will use cached results but publish properly)
        AMPLISEQ_SIMPLIFIED_RERUN(ch_samplesheet, ch_params_best)
        
        // Update the output channels to use the second run's outputs
        ch_best_tsv = AMPLISEQ_SIMPLIFIED_RERUN.out.runs_asv_table
        ch_best_fasta = AMPLISEQ_SIMPLIFIED_RERUN.out.runs_asv_fasta
    
    }

    // Create links to best runs
    ch_best_tsv.collect().map { it[0] }
        .set { ch_best_for_link }
    CREATE_LINK ( ch_best_for_link )



    /*emit:
    runs_summary     = AMPLISEQ_SIMPLIFIED.out.runs_summary
    runs_asv_table   = AMPLISEQ_SIMPLIFIED.out.runs_asv_table
    runs_asv_tax     = AMPLISEQ_SIMPLIFIED.out.runs_asv_tax
    best_asv_table   = ch_best_tsv
    best_asv_fasta   = ch_best_fasta
    comparison_results = params.skip_run_comparison ? Channel.empty() : COMPARE_RUNS.out
    multiqc_report   = AMPLISEQ_SIMPLIFIED.out.multiqc_report
    versions         = AMPLISEQ_SIMPLIFIED.out.versions*/
}

