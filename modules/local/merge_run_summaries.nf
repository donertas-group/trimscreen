process MERGE_RUN_SUMMARIES {
    tag "merge"
    label "process_low"

    conda "conda-forge::python=3.12.0 biopython=1.81 numpy=1.26.3 pandas=1.1.5 scikit-bio=0.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scikit-bio:0.4.2--np112py36_0' :
        'oras://community.wave.seqera.io/library/scikit-bio:0.6.3--60b3440d8dded0f7' }"

    input:
    path stats

    output:
    path "full_table.csv" , emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat $stats > full_table.csv
    """
}
