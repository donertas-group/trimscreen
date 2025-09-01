process NOVASEQ_ERR {
    tag "$meta.run"
    label 'process_single'

    conda "bioconda::bioconductor-dada2=1.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    publishDir "${params.outdir}/runs/${meta.runID}/dada2/QC",
        mode: params.publish_dir_mode,
        pattern: "*{.pdf,convergence.txt}",
        enabled: params.publish_all_runs

    publishDir "${params.outdir}/runs/${meta.runID}/dada2/QC/svg",
        mode: params.publish_dir_mode,
        pattern: "*.svg",
        enabled: params.publish_all_runs

    publishDir "${params.outdir}/best_run/${meta.runID}/dada2/QC",
        mode: params.publish_dir_mode,
        pattern: "*{.pdf,convergence.txt}",
        enabled: "${meta.is_best_run}"

    publishDir "${params.outdir}/best_run/${meta.runID}/dada2/QC/svg",
        mode: params.publish_dir_mode,
        pattern: "*.svg",
        enabled: "${meta.is_best_run}"


    input:
    tuple val(meta), path(errormodel)

    output:
    tuple val(meta), path("*.md.err.rds"), emit: errormodel
    tuple val(meta), path("*.md.err.pdf"), emit: pdf
    tuple val(meta), path("*.md.err.svg"), emit: svg
    tuple val(meta), path("*.md.err.convergence.txt"), emit: convergence
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (!meta.single_end) {
        """
        novaseq_err_pe.r ${errormodel[0]} ${errormodel[1]} ${meta.run}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
        END_VERSIONS
        """
    } else {
        """
        novaseq_err_se.r ${errormodel} ${meta.run}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
        END_VERSIONS
        """
    }
}
