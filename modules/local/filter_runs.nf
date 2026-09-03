process FILTER_RUNS {
    tag "filter"
    label "process_single"

    conda "conda-forge::python=3.12.0 biopython=1.81 numpy=1.26.3 pandas=1.1.5 scikit-bio=0.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scikit-bio:0.4.2--np112py36_0' :
        'oras://community.wave.seqera.io/library/scikit-bio:0.6.3--60b3440d8dded0f7' }"

    input:
    path table
    val suffix
    val min_reads
    val expected_samples   // comma-separated String of all expected sample names

    output:
    // stdout now carries three JSON lines:
    //   line 1: [[runID, depth], ...]   (good runs, as before)
    //   line 2: [sample, ...]           (samples that never clear min_reads
    //                                    in any run - empty list if none;
    //                                    error-level, stops the pipeline)
    //   line 3: [sample, ...]           (samples missing from the table
    //                                    entirely, e.g. dropped upstream by
    //                                    ASV length filtering - empty list
    //                                    if none; warning-level only)
    tuple stdout, path("samplerun_summaries${suffix}.csv"), emit: filtered

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    filter_runs.py -i $table -o samplerun_summaries${suffix}.csv -n $min_reads -s "$expected_samples"
    """
}

