include { AMPLISEQ_SIMPLIFIED                                       } from '../ampliseq_simplified/main'
include { AMPLISEQ_SIMPLIFIED as AMPLISEQ_SIMPLIFIED_RERUN } from '../ampliseq_simplified/main'
include { GENERATE_PARAMS                                           } from '../../../modules/local/generate_params'
include { COMPARE_RUNS                  } from '../../../subworkflows/local/compare_runs/main'
include { CREATE_LINK                   } from '../../../modules/local/create_link'


// Input
if (params.metadata) {
    //ch_metadata = Channel.fromPath("${params.metadata}", checkIfExists: true)
    ch_metadata = Channel
        .fromPath(params.metadata, checkIfExists: true)
        .splitCsv(header: true)
    
    ch_metadata_samples = ch_metadata.filter { row -> row.condition == "sample" }
       // .map { it.sampleID }
        .map { row -> [row.sampleID, row.is_replicate == true || row.is_replicate == 'true'] } 

} else { 
    ch_metadata = Channel.empty() 
    ch_metadata_samples = Channel.empty()
}


workflow AMPLISEQ_SCREENING {
    take:
    ch_samplesheet
    
    main:
    // generate sets of parameters based on input ranges
    FW_primer_len = (params.FW_primer && !params.skip_cutadapt) ? params.FW_primer.size() : 0
    RV_primer_len = (params.RV_primer && !params.skip_cutadapt) ? params.RV_primer.size() : 0

    trunclenf_range = params.trunclenf_range ?: ""
    trunclenr_range = params.trunclenr_range ?: ""

    ch_read_length = ch_samplesheet
            .first()
            .map { meta, readfw, readrv -> 
                // Return the first FASTQ file (readfw)
                return readfw
            }
            .splitFastq(record: true, limit: 1)
            .map { record -> record.readString.length() }
        
    // Report extracted read length:
    ch_read_length.subscribe { read_len ->
        log.info "Read length extracted from samplesheet: ${read_len}"
    }

    GENERATE_PARAMS (
        params.marker_size_min, 
        params.minimum_overlap, 
        params.screen_interval, 
        ch_read_length,
        FW_primer_len,
        RV_primer_len,
        trunclenf_range,
        trunclenr_range 
    )

    ch_params = GENERATE_PARAMS.out.params_csv
    .splitCsv(header: true, sep: ',')
    .map { row -> 
           def is_best_run = false // initiallise is_best_run to false
           return tuple(row.runID, row.trunclenf, row.trunclenr, is_best_run) } 

    // Subset samples 
    ch_is_best_run = ch_params
        .map { runID, trunclenf, trunclenr, is_best_run -> is_best_run }
        .unique()

    subset_samples = params.subset_samples ?: false

    // restrict to samples only
    if (params.metadata) {
        ch_samplesheet_samples = ch_samplesheet
            .map { meta, read1, read2 -> [meta.id, [meta, read1, read2]] }
            .join (ch_metadata_samples)
            .map { id, tuple, is_replicate -> [tuple, is_replicate] }            
    }
    else { ch_samplesheet_samples = ch_samplesheet.map { tuple -> [tuple, false] } 
    }

    ch_samplesheet_subset = ch_samplesheet_samples.collect()
        .combine (ch_is_best_run)
        .map { tuple ->
            def is_best_run = tuple[-1]
            def all_items = tuple[0..-2].collate(2)

            // DEBUG: Check the structure
            //log.info "DEBUG: all_items size = ${all_items.size()}"
            //log.info "DEBUG: first item = ${all_items[0]}"
            //log.info "DEBUG: first item class = ${all_items[0].getClass()}"

            if (subset_samples && !is_best_run) {
                // Sort so replicates come first, then shuffle within each group
                def replicates = all_items.findAll { item -> item[1] }.collect { item -> item[0] }.shuffled()
                def non_replicates = all_items.findAll { item -> !item[1] }.collect { item -> item[0] }.shuffled()
                def sorted_samples = replicates + non_replicates
                
                if (subset_samples < sorted_samples.size()) {
                    log.info "Selecting ${subset_samples} samples (prioritizing ${replicates.size()} replicates) out of ${sorted_samples.size()} total"
                    return sorted_samples.take(subset_samples)
                } else {
                    log.info "Requested ${subset_samples} samples but only ${sorted_samples.size()} available - using all samples for screening"
                    return sorted_samples
                }
            } else {
                def all_samples = all_items.collect { item -> item[0]}
                log.info "Using all ${all_samples.size()} non-control samples for screening"
                return all_samples
            }
        }
        .flatMap { it }

    // Check the number of individual runs before running the workflow
    ch_samplesheet_subset
        .combine(ch_params)
        .count()
        .subscribe { count ->
            if (count > 10000) {
                log.warn """
        ================================================================================
        WARNING: Too many sample-runs 
        ================================================================================
        Sample numbers and parameter settings result in ${count} sample-runs.
        This can take a long time to finish, or take too much computing resources.
        
        Recommendations:
        - Consider subsetting your samples further with `--subset_samples`
        - Reduce the number of parameters to be screened by 
          1. increasing `--screen_interval` or
          2. adjusting `--trunclenf_range` and/or `--trunclenr_range`
        ================================================================================
                """
            } else {
                log.info "Number of sample-runs: ${count}"
            }
        }


    // Run simplified version of ampliseq
    AMPLISEQ_SIMPLIFIED(ch_samplesheet_subset, ch_params)

    ch_runs_summary = AMPLISEQ_SIMPLIFIED.out.runs_summary
    ch_runs_asv_table = AMPLISEQ_SIMPLIFIED.out.runs_asv_table
    ch_runs_asv_tax = AMPLISEQ_SIMPLIFIED.out.runs_asv_tax
 
    //
    // SUBWORKFLOW: Compare runs 
    //
    if (!params.skip_run_comparison) {
        COMPARE_RUNS ( 
            ch_runs_summary, 
            ch_runs_asv_table, 
            ch_runs_asv_tax 
        )

        // Process comparison results to identify best run
        COMPARE_RUNS.out
            .map{ stdout, report -> stdout }
            .splitJson()
            .flatten()
            .set{ ch_best_run }

        // Filter outputs to get best run data
        ch_runs_asv_table
            .map { meta, file -> [meta.runID, meta, file] }
            .join( ch_best_run.map { runID -> [runID, true] }, by: 0)
            .map { runID, meta, file, _ -> [meta, file] }
            .set { ch_best_tsv }

       /* AMPLISEQ_SIMPLIFIED.out.runs_asv_fasta
            .map { meta, file -> [meta.runID, meta, file] }
            .join( ch_best_run.map { runID -> [runID, true] }, by: 0)
            .map { runID, meta, file, _ -> [meta, file] }
            .set { ch_best_fasta }*/

    } else {
        // If comparison is skipped, output all runs
        ch_best_tsv = ch_runs_asv_table
        //ch_best_fasta = AMPLISEQ_SIMPLIFIED.out.runs_asv_fasta
        ch_best_run = Channel.empty()
    }


    //
    //
    //

    if (!params.skip_run_comparison && !params.subset_samples && params.publish_all_runs) { 

        // Create links to best runs
        ch_best_tsv.collect().map { it[0] }
            .set { ch_best_for_link }
        CREATE_LINK ( ch_best_for_link )

    } else if (!params.skip_run_comparison) {

        // Create updated metadata channel with is_best_run = true for best runs
        ch_best_run_annotated = ch_best_run
        .map { runID -> [runID, true] }        

        ch_params_best = ch_params
            .map {runID, trunclenf, trunclenr, is_best_run -> [runID, trunclenf, trunclenr] }
            .join( ch_best_run_annotated )

        // Re-run AMPLISEQ_SIMPLIFIED with updated metadata (will use cached results but publish properly)
        AMPLISEQ_SIMPLIFIED_RERUN(ch_samplesheet, ch_params_best)
        
        // Update the output channels to use the second run's outputs
        ch_best_tsv = AMPLISEQ_SIMPLIFIED_RERUN.out.runs_asv_table
        //ch_best_fasta = AMPLISEQ_SIMPLIFIED_RERUN.out.runs_asv_fasta
    
    }



    emit:
    best_run         = ch_best_run
    multiqc_report   = AMPLISEQ_SIMPLIFIED.out.multiqc_report
    versions         = AMPLISEQ_SIMPLIFIED.out.versions
    /*
    runs_summary     = AMPLISEQ_SIMPLIFIED.out.runs_summary
    runs_asv_table   = AMPLISEQ_SIMPLIFIED.out.runs_asv_table
    runs_asv_tax     = AMPLISEQ_SIMPLIFIED.out.runs_asv_tax
    best_asv_table   = ch_best_tsv
    best_asv_fasta   = ch_best_fasta
    comparison_results = params.skip_run_comparison ? Channel.empty() : COMPARE_RUNS.out
    */
}

