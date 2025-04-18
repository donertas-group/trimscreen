/*
 * Check input samplesheet or folder and get read channels
 */

include { CUTADAPT as CUTADAPT_BASIC                        } from '../../../modules/nf-core/cutadapt/main'
//include { CUTADAPT as CUTADAPT_READTHROUGH                  } from '../../../modules/nf-core/cutadapt/main'
//include { CUTADAPT as CUTADAPT_DOUBLEPRIMER                 } from '../../../modules/nf-core/cutadapt/main'
include { CUTADAPT_SUMMARY as CUTADAPT_SUMMARY_STD          } from '../../../modules/local/ampliseq/cutadapt_summary'
//include { CUTADAPT_SUMMARY as CUTADAPT_SUMMARY_DOUBLEPRIMER } from '../../../modules/local/ampliseq/cutadapt_summary'
include { CUTADAPT_SUMMARY_MERGE                            } from '../../../modules/local/ampliseq/cutadapt_summary_merge'

workflow CUTADAPT_WORKFLOW {
    take:
    ch_file
    illumina_pe_its
    double_primer

    main:
    ch_versions_cutadapt_workflow = Channel.empty()

    CUTADAPT_BASIC ( ch_file ).reads.set { ch_trimmed_reads }
    CUTADAPT_BASIC.out.log
        .map {
            info, log ->
                def meta = [:]
                meta.single_end = info.single_end
                meta.runID = info.runID // add runID to meta
                [ meta, log ] }
        .groupTuple(by: 0 )
        .set { ch_cutadapt_logs }
    ch_versions_cutadapt_workflow = ch_versions_cutadapt_workflow.mix( CUTADAPT_BASIC.out.versions )
    CUTADAPT_SUMMARY_STD ( "cutadapt_standard", ch_cutadapt_logs )
    ch_versions_cutadapt_workflow = ch_versions_cutadapt_workflow.mix( CUTADAPT_SUMMARY_STD.out.versions )

    if (illumina_pe_its) {
        CUTADAPT_READTHROUGH ( ch_trimmed_reads ).reads.set { ch_trimmed_reads }
    }

    if (double_primer) {
        CUTADAPT_DOUBLEPRIMER ( ch_trimmed_reads ).reads.set { ch_trimmed_reads }
        CUTADAPT_DOUBLEPRIMER.out.log
            .map {
                info, log ->
                    def meta = [:]
                    meta.single_end = info.single_end
                    [ meta, log ] }
            .groupTuple(by: 0 )
            .set { ch_cutadapt_doubleprimer_logs }
        CUTADAPT_SUMMARY_DOUBLEPRIMER ( "cutadapt_doubleprimer", ch_cutadapt_doubleprimer_logs )
        ch_summaries = CUTADAPT_SUMMARY_STD.out.tsv.combine( CUTADAPT_SUMMARY_DOUBLEPRIMER.out.tsv )
        ch_versions_cutadapt_workflow = ch_versions_cutadapt_workflow.mix( CUTADAPT_SUMMARY_DOUBLEPRIMER.out.versions )
        CUTADAPT_SUMMARY_MERGE ( "merge", ch_summaries )
        ch_versions_cutadapt_workflow = ch_versions_cutadapt_workflow.mix( CUTADAPT_SUMMARY_MERGE.out.versions )
    } else {
        CUTADAPT_SUMMARY_MERGE ( "copy", CUTADAPT_SUMMARY_STD.out.tsv )
    }

    //Filter empty files
    ch_trimmed_reads
        .branch {
            failed: it[0].single_end ? it[1].countFastq() < params.min_read_counts : it[1][0].countFastq() < params.min_read_counts || it[1][1].countFastq() < params.min_read_counts
            passed: true
        }
        .set { ch_trimmed_reads_result }


/// Added by Yi Wang to allow filtering runs that contain failed samples /////////////////////////////////////////
    ch_failed = ch_trimmed_reads_result.failed   ///////// currently ch_failed is not used because adding it causes error if no samples failed ////
        .map { meta, reads -> meta.run }
        .unique()
        .collect()
        .toList()

    // in case ch_failed is empty, create a dummy entry
    ch_failed_runs = Channel.value([[run:''],'']) 
        .map { meta, reads -> meta.run }
        .unique()
        .collect()
        .toList()
   // ch_failed_runs.mix(ch_failed ).set {ch_failed_runs}  ///// if run, error////
 
    // Filter out the failed runs from the passed samples
    ch_trimmed_reads_passed = ch_trimmed_reads_result.passed
       // .combine(ch_failed_runs)
       // .filter {meta, reads, failed_runs -> !(meta.run in failed_runs) }
       // .map { meta, reads, failed_runs -> [meta, reads] }

    // Log the failed runs
    ch_failed
        .subscribe { failed_runs ->
            if (failed_runs.size() > 0) {
                log.warn "The following runs contained samples with too few reads (<$params.min_read_counts) after trimming and will be filtered out:\n${failed_runs.join('\n')}"
            }
        }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ch_trimmed_reads_result.failed
        .map { meta, reads -> [ meta.id ] }
        .collect()
        .subscribe {
            samples = it.join("\n")
            if (params.ignore_failed_trimming) {
                log.warn "The following samples had too few reads (<$params.min_read_counts) after trimming with cutadapt:\n$samples\nIgnoring failed samples and continue!\n"
            } else {
                error("The following samples had too few reads (<$params.min_read_counts) after trimming with cutadapt:\n$samples\nPlease check whether the correct primer sequences for trimming were supplied. Ignore that samples using `--ignore_failed_trimming` or adjust the threshold with `--min_read_counts`.")
            }
        }


    emit:
    reads    = ch_trimmed_reads_passed
    logs     = CUTADAPT_BASIC.out.log
    summary  = CUTADAPT_SUMMARY_MERGE.out.tsv
    versions = ch_versions_cutadapt_workflow
}
