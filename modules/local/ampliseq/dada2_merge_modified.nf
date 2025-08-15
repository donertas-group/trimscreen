process DADA2_MERGE {
    tag "$meta.run"
    label 'process_single'

    // https://depot.galaxyproject.org/singularity/bioconductor-dada2=1.28.0--r43hf17093f_0 doesnt contain 'digest', so keep here v1.22.0
    conda "bioconda::bioconductor-dada2=1.22.0 conda-forge::r-digest=0.6.30"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    publishDir "${params.outdir}/runs/${meta.runID}/dada2",
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
        enabled: "${meta.is_best_run || params.publish_all_runs}"


    input:
    tuple val(meta), path(files)
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path( "DADA2_stats.tsv.gz" ), emit: dada2stats
    tuple val(meta), path( "DADA2_table.tsv.gz" ), emit: dada2asv
    tuple val(meta), path( "ASV_table.tsv.gz" ),   emit: asv
    tuple val(meta), path( "ASV_seqs.fasta.gz" ) , emit: fasta
    tuple val(meta), path( "DADA2_table.rds" ), emit: rds
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    suppressPackageStartupMessages(library(digest))

    #combine stats files
    for (data in sort(list.files(".", pattern = ".stats.tsv", full.names = TRUE))) {
        if (!exists("stats")){ stats <- read.csv(data, header=TRUE, sep="\\t") }
        if (exists("stats")){
            temp <-read.csv(data, header=TRUE, sep="\\t")
            stats <-unique(rbind(stats, temp))
            rm(temp)
        }
    }
    write.table( stats, file = gzfile("DADA2_stats.tsv.gz"), sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    #combine dada-class objects
    files <- sort(list.files(".", pattern = ".ASVtable.rds", full.names = TRUE))
    if ( length(files) == 1 ) {
        ASVtab = readRDS(files[1])
    } else {
        ASVtab <- mergeSequenceTables(tables=files, repeats = "error", orderBy = "abundance", tryRC = FALSE)
    }
    saveRDS(ASVtab, "DADA2_table.rds")

    df <- t(ASVtab)
    colnames(df) <- gsub('_1.filt.fastq.gz', '', colnames(df))
    colnames(df) <- gsub('.filt.fastq.gz', '', colnames(df))
    colnames(df) <- gsub('\\\\.run.*', '', colnames(df)) # added this line to remove run numbers in table
    df <- data.frame(sequence = rownames(df), df, check.names=FALSE)
    # Create an md5 sum of the sequences as ASV_ID and rearrange columns
    df\$ASV_ID <- sapply(df\$sequence, digest, algo='md5', serialize = FALSE)
    df <- df[,c(ncol(df),3:ncol(df)-1,1)]

    # file to publish
    write.table(df, file = gzfile("DADA2_table.tsv.gz"), sep = "\\t", row.names = FALSE, quote = FALSE, na = '')

    # Write fasta file with ASV sequences to file
    write.table(data.frame(s = sprintf(">%s\n%s", df\$ASV_ID, df\$sequence)), gzfile('ASV_seqs.fasta.gz'), col.names = FALSE, row.names = FALSE, quote = FALSE, na = '')

    # Write ASV file with ASV abundances to file
    df\$sequence <- NULL
    write.table(df, file = gzfile("ASV_table.tsv.gz"), sep="\\t", row.names = FALSE, quote = FALSE, na = '')

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}
