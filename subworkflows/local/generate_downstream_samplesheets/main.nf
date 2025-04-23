process WRITE_SAMPLESHEET {
    publishDir "${params.outdir}/downstream_samplesheets", mode: 'copy'

    input:
    val format
    val sep
    val header
    val sample_info_list

    output:
    path "samplesheet.${format}", emit: samplesheet

    script:
    def lines = [header.join(sep)] + sample_info_list.collect { it.values().join(sep) }
    def content = lines.join('\n')
    
    """
    echo '${content}' > samplesheet.${format}
    """
}

workflow GENERATE_DOWNSTREAM_SAMPLESHEETS {
    take:
    ch_reads

    main:
    def format = 'csv' // most common format in nf-core

    ch_list_for_samplesheet = ch_reads
        .map { meta, reads ->
            def out_path = file(params.outdir).toString() + '/filter/filtered/'
            def sampleID = meta.id
            def forwardReads = meta.single_end ? out_path + reads.getName() : out_path + reads[0].getName()
            def reverseReads = !meta.single_end ? out_path + reads[1].getName() : ""
            [sampleID: sampleID, forwardReads: forwardReads, reverseReads: reverseReads]
        }

    def header = ['sampleID', 'forwardReads', 'reverseReads']
    
    WRITE_SAMPLESHEET(
        format,
        ["csv":",", "tsv":"\t", "txt":"\t"][format],
        header,
        ch_list_for_samplesheet.collect()
    )

    emit:
    samplesheet = WRITE_SAMPLESHEET.out.samplesheet
}
