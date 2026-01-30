process RAREFY_RUNS {
    tag "rarefy $meta.runID"
    label "process_single"

    conda "conda-forge::python=3.12.0 biopython=1.81 numpy=1.26.3 pandas=1.1.5 scikit-bio=0.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scikit-bio:0.4.2--np112py36_0' :
        'oras://community.wave.seqera.io/library/scikit-bio:0.6.3--60b3440d8dded0f7' }"

    publishDir "${params.outdir}/runs/${meta.runID}/dada2",
        mode: params.publish_dir_mode,
        pattern: "*.tsv.gz"


    input:
    tuple val(meta), path(asv_table), val(depth)

    output:
    tuple val(meta), path("ASV_table_rarefied.tsv.gz"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    rarefy_runs.py -i $asv_table -d $depth
    """
}

