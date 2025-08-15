process PICRUST {
    tag "${seq},${abund},${meta.runID}"
    label 'process_single'

    conda "bioconda::picrust2=2.5.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picrust2:2.5.3--pyhdfd78af_0' :
        'biocontainers/picrust2:2.5.3--pyhdfd78af_0' }"

    publishDir "${params.outdir}/runs/${meta.runID}/picrust",
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
        enabled: "${meta.is_best_run || params.publish_all_runs}"

    input:
    tuple val(meta), path(seq)
    tuple val(meta), path(abund)
    val(source)
    val(message)

    output:
    tuple val(meta), path("all_output/*")    , emit: outfolder
    tuple val(meta), path("*_descrip.tsv")   , emit: pathways
    path "versions.yml"                      , emit: versions
    tuple val(meta), path("*.args.txt")      , emit: args
    tuple val(meta), path("message.txt")  , emit: message

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    #If input is QIIME2 file, than (1) the first line and (2) the first character (#) of the second line need to be removed
    if [ "$source" == 'QIIME2' ]
    then
        tail -n +2 "$abund" > "${abund}.tmp" && mv "${abund}.tmp" "$abund"
    fi

    picrust2_pipeline.py \\
        $args \\
        -s $seq \\
        -i $abund \\
        -o all_output \\
        -p $task.cpus \\
        --in_traits EC,KO \\
        --verbose

    #Add descriptions to identifiers
    add_descriptions.py -i all_output/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o EC_pred_metagenome_unstrat_descrip.tsv
    add_descriptions.py -i all_output/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
                    -o KO_pred_metagenome_unstrat_descrip.tsv
    add_descriptions.py -i all_output/pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o METACYC_path_abun_unstrat_descrip.tsv

    echo "$message" > "message.txt"
    echo -e "picrust\t$args" > "picrust.args.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        picrust2: \$( picrust2_pipeline.py -v | sed -e "s/picrust2_pipeline.py //g" )
    END_VERSIONS
    """
}
