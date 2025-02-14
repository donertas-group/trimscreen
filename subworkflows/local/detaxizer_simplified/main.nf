//
// Subworkflow with functionality to do host/contaminant-removal with modules from the nf-core/detaxizer pipeline
//

workflow DETAXIZER_SIMPLIFIED {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    ch_short = ch_samplesheet.branch {
        shortReads: it[1]
        }.shortReads.map{
        meta, short_reads_fastq_1, short_reads_fastq_2, long_reads_fastq_1 ->
            if (short_reads_fastq_2){
                return [meta + [ single_end: false, long_reads: false , amount_of_files: 2 ], [ short_reads_fastq_1, short_reads_fastq_2 ] ]
            } else {
                return [meta + [ id: "${meta.id}_R1", single_end: true, long_reads: false, amount_of_files: 1 ], short_reads_fastq_1 ]
            }
    }

    ch_long = ch_samplesheet.branch {
        longReads: it[3]
    }.longReads.map {
        meta, short_reads_fastq_1, short_reads_fastq_2, long_reads_fastq_1 ->
            return [meta + [ id: "${meta.id}_longReads", single_end: true, long_reads: true, amount_of_files: 1 ], long_reads_fastq_1 ]
    }

    ch_short_long = ch_short.mix(ch_long)

    //
    // MODULE: Rename Fastq headers
    //
    RENAME_FASTQ_HEADERS_PRE(ch_short_long)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        RENAME_FASTQ_HEADERS_PRE.out.fastq
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_fastq_for_classification = RENAME_FASTQ_HEADERS_PRE.out.fastq



    //////////////////////////////////////////////////
    //  Classification
    //////////////////////////////////////////////////


    //
    // MODULE: Run bbduk
    //
    BBMAP_BBDUK (
        ch_fastq_for_classification,
        ch_fasta_bbduk.first()
    )
    ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions.first())

    //
    // MODULE: Run ISOLATE_BBDUK_IDS
    //
    ISOLATE_BBDUK_IDS(
        BBMAP_BBDUK.out.contaminated_reads
    )
    ch_versions = ch_versions.mix(ISOLATE_BBDUK_IDS.out.versions.first())


    //
    // MODULE: Merge IDs
    //
    MERGE_IDS(
        ISOLATE_BBDUK_IDS.out.classified_ids
    )
    ch_versions = ch_versions.mix(MERGE_IDS.out.versions.first())


    //
    // MODULE: Summarize the classification results
    //

    SUMMARY_CLASSIFICATION(
        MERGE_IDS.out.classified_ids
    )

    // Drop meta of kraken2_summary as it is not needed for the combination step of summarizer
    ch_classification_summary = SUMMARY_CLASSIFICATION.out.summary.map {
            meta, path -> [path]
    }
    ch_versions = ch_versions.mix(SUMMARY_CLASSIFICATION.out.versions.first())



    //////////////////////////////////////////////////
    //  Validation
    //////////////////////////////////////////////////

    if (params.validation_blastn) {

        //
        // MODULE: Extract the hits to fasta format
        //
        ch_combined = ch_fastq_for_classification
        .join(
            MERGE_IDS.out.classified_ids, by: [0]
        )

        PREPARE_FASTA4BLASTN (
            ch_combined
        )

        ch_versions = ch_versions.mix(PREPARE_FASTA4BLASTN.out.versions.first())

        //
        // MODULE: Run BLASTN
        //
        ch_reference_fasta = ch_fasta_blastn

        ch_reference_fasta_with_meta = ch_reference_fasta.map {
            item -> [['id': "id-fasta-for-makeblastdb"], item]
            }

        BLAST_MAKEBLASTDB (
                ch_reference_fasta_with_meta
        )
        ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

        ch_fasta4blastn = PREPARE_FASTA4BLASTN.out.fasta
            .flatMap { meta, fastaList ->
                if (fastaList.size() == 2) {
                return [
                    [ meta + [ id: "${meta.id}_R1" ], fastaList[0] ],
                    [ meta + [ id: "${meta.id}_R2" ], fastaList[1] ]
                ]

                } else {
                    return [
                        [ meta , fastaList ] ]
                }
            }

        ch_blastn_db = BLAST_MAKEBLASTDB.out.db.first()

        BLAST_BLASTN (
            ch_fasta4blastn,
            ch_blastn_db
        )

        ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

        ch_combined_blast = BLAST_BLASTN.out.txt.map {
            meta, path ->
                return [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
        }
        .map{
            meta, path -> tuple(groupKey(meta, meta.amount_of_files), path)
        }
        .groupTuple(
                by: [0]
            ).map {
                meta, paths -> [ meta, paths.flatten() ]
                }

        FILTER_BLASTN_IDENTCOV (
            BLAST_BLASTN.out.txt
        )
        ch_versions = ch_versions.mix(FILTER_BLASTN_IDENTCOV.out.versions.first())

        ch_filtered_combined = FILTER_BLASTN_IDENTCOV.out.classified.map {
            meta, path ->
                return [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
        }
        .map{
            meta, path -> tuple(groupKey(meta, meta.amount_of_files), path)
        }
        .groupTuple (by: [0])
        .map {
            meta, paths ->
                paths = paths.flatten()
                return [ meta, paths ]
        }

        ch_blastn_combined = ch_combined_blast.join(ch_filtered_combined, remainder: true).map{
            meta, blastn, filteredblastn ->
                if (blastn[0] == null){
                    blastn[0] = []
                }
                if (blastn[1] == null){
                    blastn[1] = []
                }
                if (filteredblastn[0] == null){
                    filteredblastn[0] = []
                }
                if (filteredblastn[1] == null){
                    filteredblastn[1] = []
                }
                return [ meta, blastn[0], blastn[1], filteredblastn[0], filteredblastn[1] ]
            }
        ch_blastn_summary = SUMMARY_BLASTN (
            ch_blastn_combined
        )
        ch_versions = ch_versions.mix(ch_blastn_summary.versions.first())

    // Drop meta of blastn_summary as it is not needed for the combination step of summarizer
        ch_blastn_summary = ch_blastn_summary.summary.map {
                meta, path -> [path]
            }
        }


    //
    // MODULE: Filter out the classified or validated reads
    //
    if ( !params.validation_blastn && params.enable_filter ) {

        ch_classification = RENAME_FASTQ_HEADERS_PRE.out.fastq
            .join(MERGE_IDS.out.classified_ids, by:[0])

        FILTER(
            ch_classification
        )

        ch_versions = ch_versions.mix(FILTER.out.versions.first())

    } else if ( params.enable_filter ) {

        ch_blastn2filter = FILTER_BLASTN_IDENTCOV.out.classified_ids.map {
            meta, path ->
                return [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
        }
        .map{
            meta, path -> tuple(groupKey(meta, meta.amount_of_files), path)
        }
        .groupTuple(by:[0])

        ch_combined_short_long_id = RENAME_FASTQ_HEADERS_PRE.out.fastq.map {
            meta, path ->
                return [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
        }

        ch_blastnfilter = ch_combined_short_long_id.join(
            ch_blastn2filter, by:[0]
        )

        FILTER(
            ch_blastnfilter
        )

        ch_versions = ch_versions.mix(FILTER.out.versions.first())
    
    }

    //
    // MODULE: Rename headers after filtering
    //
    if ( params.enable_filter ) {

    ch_headers = RENAME_FASTQ_HEADERS_PRE.out.headers.map {
        meta, path ->
            return [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
    }

    ch_filtered2rename = FILTER.out.filtered.map {
        meta, path ->
            return [ meta + [ id: meta.id.replaceAll("(_R1|_R2)", "") ], path ]
    }

    ch_removed2rename = Channel.empty()

    ch_rename_filtered = ch_filtered2rename.join(ch_headers, by:[0])

    ch_removed2rename = ch_removed2rename.ifEmpty(['empty', []])


    RENAME_FASTQ_HEADERS_AFTER(
        ch_rename_filtered,
        ch_removed2rename.first()
    )
    ch_versions = ch_versions.mix(RENAME_FASTQ_HEADERS_AFTER.out.versions.first())
    }


    //
    // MODULE: Summarize the classification process
    //
    if (params.validation_blastn){

    ch_summary = ch_classification_summary.mix(ch_blastn_summary).collect().map {
            item -> [['id': "summary_of_classification_and_blastn"], item]
        }
    } else {
        ch_summary = ch_classification_summary.collect().map {
            item -> [['id': "summary_of_classification"], item]
        }
    }

    ch_summary = SUMMARIZER (ch_summary)
    ch_versions = ch_versions.mix(ch_summary.versions)
    
    if (params.generate_downstream_samplesheet){
        GENERATE_DOWNSTREAM_SAMPLESHEETS( RENAME_FASTQ_HEADERS_AFTER.out.fastq )
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
