process CREATE_LINK {
    tag "$meta.runID"
    publishDir "$params.outdir/best_run", mode: 'link', overwrite: true

    input:
    val(meta)

    output:
    path "${meta.runID}", emit: dir

    script:
    """
    ln -s ${workflow.projectDir}/../$params.outdir/runs/${meta.runID} ${meta.runID}
    """
}
