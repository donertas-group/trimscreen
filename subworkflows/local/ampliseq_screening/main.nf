include { AMPLISEQ_SIMPLIFIED                                       } from '../ampliseq_simplified/main'
include { AMPLISEQ_SIMPLIFIED as AMPLISEQ_SIMPLIFIED_RERUN } from '../ampliseq_simplified/main'
include { GENERATE_PARAMS                                           } from '../../../modules/local/generate_params'
include { SUMMARISE_RUNS                } from '../../../subworkflows/local/summarise_runs/main'
include { COMPARE_RUNS                  } from '../../../subworkflows/local/compare_runs/main'
include { CREATE_LINK                   } from '../../../modules/local/create_link'

// Input
if (params.metadata) {
    //ch_metadata = Channel.fromPath("${params.metadata}", checkIfExists: true)
    ch_metadata = Channel
        .fromPath(params.metadata, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> // Validation step using .map
            def required_columns = ['sampleID', 'condition', 'replicated']
            def missing_columns = required_columns.findAll { !row.containsKey(it) }
            
            if (missing_columns) {
                error "ERROR: The metadata file is missing mandatory column(s): ${missing_columns.join(', ')}. Please check your input."
            }
            return row
        }
    ch_metadata_samples = ch_metadata.filter { row -> row.condition == "sample" }
        .map { row -> 
            // Normalize the string: lowercase it and take the first character
            def val = row.replicated?.toString()?.toLowerCase()?.trim()
            def isReplicated = (val == 'true' || val == 't')
            
            return [row.sampleID, isReplicated]
        }

} else { 
    ch_metadata = Channel.empty() 
    ch_metadata_samples = Channel.empty()
}

if ( params.skip_trim_screening){ 
    if (params.trunclenf_range && params.trunclenr_range) {
    log.warn "`skip_trim_screening` is set to true. `trunclenf_range` and `trunclenr_range` will be ignored." 
}
}
if ( params.skip_trim_screening){ 
    if ( params.trunclenf && params.trunclenr){ 
        log.warn "`skip_trim_screening` is set to true, trimming will be performed as set by `trunclenf = ${params.trunclenf}` and `trunclenr = ${params.trunclenr}`."
    }
	else if (!params.skip_dada_quality) {
        log.warn "`skip_trim_screening` is set to true and `skip_dada_quality` is set to false. DADA2 quality rimming will be performed."
    }
    else {
        log.warn "Both `skip_trim_screening` and `skip_dada_quality` are set to true. No rear-end read trimming will be performed." 
    } 
}

ch_params = Channel.empty()

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

//
// Do trimscreening
//

    if (!params.skip_trim_screening) {


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
            def run_type = "screened" // initiallise run_type to "screened"
            return tuple(row.run, row.trunclenf, row.trunclenr, run_type) } 

    // Subset samples 
    ch_is_best_run = ch_params
        .map { run, trunclenf, trunclenr, run_type -> run_type }
        .unique()

    subset_samples = params.subset_samples ?: false

    // restrict to samples only
    if (params.metadata) {
        ch_samplesheet_samples = ch_samplesheet
            .map { meta, read1, read2 -> [meta.id, [meta, read1, read2]] }
            .join (ch_metadata_samples)
            .map { id, tuple, replicated -> [tuple, replicated] }            
    
    } else {     
        ch_samplesheet_samples = ch_samplesheet.map { tuple -> [tuple, false] } 
    }

    ch_samplesheet_subset = ch_samplesheet_samples.collect()
        .combine (ch_is_best_run)
        .map { tuple ->
            def run_type = tuple[-1]
            def all_items = tuple[0..-2].collate(2)

            // DEBUG: Check the structure
            //log.info "DEBUG: all_items size = ${all_items.size()}"
            //log.info "DEBUG: first item = ${all_items[0]}"
            //log.info "DEBUG: first item class = ${all_items[0].getClass()}"

            if (subset_samples && run_type == "screened" ) {
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

}


//
// Trimscreening is skipped
//

    else if ( params.trunclenf && params.trunclenr){
        AMPLISEQ_SIMPLIFIED(ch_samplesheet, ch_params) // ch_params set by params.trunclenf && params.trunclenr
    }
    else if (!params.skip_dada_quality) {
        AMPLISEQ_SIMPLIFIED(ch_samplesheet, ch_params) // ch_params set by dada default quality trimming
    }
    else {

        def run       = "no_trimming"
        def trunclenf = 0
        def trunclenr = 0
        def run_type  = "no_trimming"
        
        ch_params = Channel.of(tuple(run, trunclenf, trunclenr, run_type))

        AMPLISEQ_SIMPLIFIED(ch_samplesheet, ch_params) // ch_params has no trimming, just up- and downstream steps
    }


    ch_runs_summary   = AMPLISEQ_SIMPLIFIED.out.runs_summary
    ch_runs_asv_table = AMPLISEQ_SIMPLIFIED.out.runs_asv_table
    ch_runs_asv_tax   = AMPLISEQ_SIMPLIFIED.out.runs_asv_tax
    ch_run_qtrim      = AMPLISEQ_SIMPLIFIED.out.run_qtrim
 
    // Prepare structures that SUMMARISE needs regardless of whether we skip comparison
    ch_run_data = ch_runs_summary
        .combine(ch_runs_asv_table, by:0)
        .combine(ch_runs_asv_tax, by:0)

    if (params.metadata) {
        ch_metadata = Channel.fromPath("${params.metadata}", checkIfExists: true)
    } else {
        ch_metadata = Channel.value(file('NO_FILE'))
    }

    // 1. ALWAYS SUMMARIZE BY DEFAULT
    SUMMARISE_RUNS(ch_run_data, ch_metadata)
    def full_table = SUMMARISE_RUNS.out.full_table

    // 2. CONDITIONALLY RUN COMPARISONS
    if (!params.skip_run_comparison) {
        
        def run_qtrim_string = ch_run_qtrim
            .ifEmpty([[run: ''], []])
            .map { meta, file -> meta.run }
            .first()

        COMPARE_RUNS ( 
            full_table,
            ch_run_data,
            ch_runs_summary,
            ch_runs_asv_tax,
            run_qtrim_string,
            ch_metadata
        )

        // Process comparison results to identify best run
        COMPARE_RUNS.out.best_runs
            .map{ stdout, report -> stdout }
            .splitJson()
            .flatten()
            .set{ ch_best_run }

        // Filter outputs to get best run data
        ch_runs_asv_table
            .map { meta, file -> [meta.run, meta, file] }
            .join( ch_best_run.map { run -> [run, "suggested"] }, by: 0)
            .map { run, meta, file, _ -> [meta, file] }
            .set { ch_best_tsv }

    } else {
        // If comparison is skipped, output all runs or default channels empty
        ch_best_tsv = ch_runs_asv_table
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

        // Create updated metadata channel with run_type = "suggested" for best runs
        ch_best_run_annotated = ch_best_run
            .map { run -> [run, "suggested"] }        

        ch_params_best = ch_params
            .map {run, trunclenf, trunclenr, run_type -> [run, trunclenf, trunclenr] }
            .join( ch_best_run_annotated )

        // Re-run AMPLISEQ_SIMPLIFIED with updated metadata (will use cached results but publish properly)
        //AMPLISEQ_SIMPLIFIED_RERUN(ch_samplesheet, ch_params_best)
        
        // Update the output channels to use the second run's outputs
        //ch_best_tsv = AMPLISEQ_SIMPLIFIED_RERUN.out.runs_asv_table
    
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

