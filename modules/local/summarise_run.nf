process SUMMARISE_RUN {
    tag "$meta.run"
    label "process_single"

    conda "conda-forge::python=3.12.0 biopython=1.81 numpy=1.26.3 pandas=1.1.5 scikit-bio=0.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scikit-bio:0.4.2--np112py36_0' :
        'oras://community.wave.seqera.io/library/scikit-bio:0.6.3--60b3440d8dded0f7' }"

    input:
    tuple val(meta), path(stats), path(asv), path(tax)
    each path(metadata)

    output:
    tuple val(meta), path("${meta.run}.samplerun_summary.csv"), path("${meta.run}.run_summary.csv"), emit: csv
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def metadata_arg = metadata.name != 'NO_FILE' ? "-m ${metadata}" : ""
    """
    summarise_run.py -i $stats $asv $tax ${meta.run} ${metadata_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """


}


