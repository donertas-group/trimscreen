process CREATE_LINK {
    tag "$meta.run"
    publishDir "$params.outdir/best_run", mode: 'link', overwrite: true

    input:
    val(meta)

    output:
    path "${meta.run}", emit: dir

    script:
    """
    if [[ "$params.outdir" = /* ]]; then
        ln -s $params.outdir/runs/${meta.run} ${meta.run}
    else
        ln -s ${workflow.projectDir}/../$params.outdir/runs/${meta.run} ${meta.run}
    fi
    """
}
