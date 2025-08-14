process FILTER_LEN {
    tag "${fasta},$meta.run"
    label 'process_single'

    conda "bioconda::bioconductor-biostrings=2.58.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-biostrings:2.58.0--r40h037d062_0' :
        'biocontainers/bioconductor-biostrings:2.58.0--r40h037d062_0' }"

    publishDir "${params.outdir}/runs/${meta.runID}/asv_length_filter",
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
        enabled:  "${meta.is_best_run}"



    input:
    tuple val(meta), path(fasta), path(table)

    output:
    tuple val(meta), path( "stats.len.tsv.gz" )      , emit: stats
    tuple val(meta), path( "ASV_table.len.tsv.gz" )  , emit: asv, optional: true
    tuple val(meta), path( "ASV_seqs.len.fasta.gz" ) , emit: fasta
    tuple val(meta), path( "ASV_len_orig.tsv.gz" )   , emit: len_orig
    tuple val(meta), path( "ASV_len_filt.tsv.gz" )   , emit: len_filt
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def min_len_asv = task.ext.min_len_asv ?: '1'
    def max_len_asv = task.ext.max_len_asv ?: '1000000'

    def read_table  = table ? "table <- read.table(file = '$table', sep = '\t', comment.char = '', header=TRUE)" : "table <- data.frame(matrix(ncol = 1, nrow = 0))"
    def asv_table_filtered  = table ? "ASV_table.len.tsv.gz" : "empty_ASV_table.len.tsv.gz"
    """
    #!/usr/bin/env Rscript

    #load packages
    suppressPackageStartupMessages(library(Biostrings))

    #read abundance file, first column is ASV_ID
    $read_table
    colnames(table)[1] <- "ASV_ID"

    #read fasta file of ASV sequences
    seq <- readDNAStringSet("$fasta")
    seq <- data.frame(ID=names(seq), sequence=paste(seq))

    #filter
    filtered_seq <- seq[nchar(seq\$sequence) %in% $min_len_asv:$max_len_asv,]
    list <- filtered_seq[, "ID", drop = FALSE]
    filtered_table <- merge(table, list, by.x="ASV_ID", by.y="ID", all.x=FALSE, all.y=TRUE)

    #report
    distribution_before <- table(nchar(seq\$sequence))
    distribution_before <- data.frame(Length=names(distribution_before),Counts=as.vector(distribution_before))
    distribution_after <- table(nchar(filtered_seq\$sequence))
    distribution_after <- data.frame(Length=names(distribution_after),Counts=as.vector(distribution_after))

    #write
    write.table(filtered_table, file = gzfile("$asv_table_filtered"), row.names=FALSE, sep="\t", col.names = TRUE, quote = FALSE, na = '')
    write.table(data.frame(s = sprintf(">%s\n%s", filtered_seq\$ID, filtered_seq\$sequence)), gzfile('ASV_seqs.len.fasta.gz'), col.names = FALSE, row.names = FALSE, quote = FALSE, na = '')
    write.table(distribution_before, file = gzfile("ASV_len_orig.tsv.gz"), row.names=FALSE, sep="\t", col.names = TRUE, quote = FALSE, na = '')
    write.table(distribution_after, file = gzfile("ASV_len_filt.tsv.gz"), row.names=FALSE, sep="\t", col.names = TRUE, quote = FALSE, na = '')

    #stats
    stats <- as.data.frame( t( rbind( colSums(table[-1]), colSums(filtered_table[-1]) ) ) )
    stats\$ID <- rownames(stats)
    colnames(stats) <- c("lenfilter_input","lenfilter_output", "sample")
    write.table(stats, file = gzfile("stats.len.tsv.gz"), row.names=FALSE, sep="\t")

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    Biostrings: ", packageVersion("Biostrings")) ), "versions.yml")
    """
}
