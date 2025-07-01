process WRITE_SINGLE_STATS {
    tag "$meta.run"
    label 'process_single'

    conda "bioconda::bioconductor-dada2=1.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    input:
    tuple val(meta), path('file.tsv')

    output:
    tuple val(meta), path("overall_summary.tsv") , emit: tsv
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript
    
    x <- read.table(\"file.tsv\", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    #write
    write.table(x, file = \"overall_summary.tsv\", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")) ), "versions.yml")
    """
}
