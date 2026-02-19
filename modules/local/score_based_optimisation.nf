process SCORE_BASED_OPTIMISATION {
    tag "decide"
    label "process_single"

    conda "conda-forge::python=3.12.0 biopython=1.81 numpy=1.26.3 pandas=1.1.5 scikit-bio=0.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scikit-bio:0.4.2--np112py36_0' :
        'oras://community.wave.seqera.io/library/scikit-bio:0.6.3--60b3440d8dded0f7' }"

    input:
    path table
    path metadata

    output:
    tuple stdout, path("report.txt"), emit: info          

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    score_based_optimisation.py -i $table -m $metadata --metrics retained_reads_percent shannon_Genus rep_similarity --metric_directions + + + --metric_weights 0.5 1 0
    """

}
