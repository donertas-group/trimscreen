process FIND_BEST_RUN {
    tag "decide"
    label "process_low"

    conda "conda-forge::python=3.12.0 biopython=1.81 numpy=1.26.3 pandas=1.1.5 scikit-bio=0.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scikit-bio:0.4.2--np112py36_0' :
        'oras://community.wave.seqera.io/library/scikit-bio:0.6.3--60b3440d8dded0f7' }"

    input:
    path table
    path metadata

    output:
    path "best_runs.csv",   emit: id
    path "report.txt",      emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    find_best_run.py -i $table -m $metadata -t Phylum Genus 
    """



}
