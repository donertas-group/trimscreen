process DADA2_ERR {
    tag "$meta.run"
    label 'process_low'

    conda "bioconda::bioconductor-dada2=1.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    publishDir "${params.outdir}/runs/${meta.run}/dada2/QC",
        mode: params.publish_dir_mode,
        pattern: "*{.pdf,convergence.txt}",
        enabled: params.publish_all_runs

    publishDir "${params.outdir}/runs/${meta.run}/dada2/QC/svg",
        mode: params.publish_dir_mode,
        pattern: "*.svg",
        enabled: params.publish_all_runs

    publishDir "${params.outdir}/runs/${meta.run}/dada2/args",
        mode: params.publish_dir_mode,
        pattern: "*.args.txt",
        enabled: params.publish_all_runs

    publishDir "${params.outdir}/runs/${meta.run}/dada2/log",
        mode: params.publish_dir_mode,
        pattern: "*.log",
        enabled: params.publish_all_runs

    publishDir "${params.outdir}/best_run/${meta.run}/dada2/QC",
        mode: params.publish_dir_mode,
        pattern: "*{.pdf,convergence.txt}",
        enabled: "${meta.run_type=='suggested'}"

    publishDir "${params.outdir}/best_run/${meta.run}/dada2/QC/svg",
        mode: params.publish_dir_mode,
        pattern: "*.svg",
        enabled: "${meta.run_type=='suggested'}"

    publishDir "${params.outdir}/best_run/${meta.run}/dada2/args",
        mode: params.publish_dir_mode,
        pattern: "*.args.txt",
        enabled: "${meta.run_type=='suggested'}"

    publishDir "${params.outdir}/best_run/${meta.run}/dada2/log",
        mode: params.publish_dir_mode,
        pattern: "*.log",
        enabled: "${meta.run_type=='suggested'}"



    input:
    tuple val(meta), path(reads)

    output:
    // RDS is mandatory for downstream, but we handle its "emptiness" in the workflow
    tuple val(meta), path("*.err.rds")            , emit: errormodel
    // All these must be optional so Nextflow doesn't crash if the R block is skipped
    tuple val(meta), path("*.err.pdf")            , emit: pdf, optional: true
    tuple val(meta), path("*.err.svg")            , emit: svg, optional: true
    tuple val(meta), path("*.err.log")            , emit: log, optional: true
    tuple val(meta), path("*.err.convergence.txt"), emit: convergence, optional: true
    path "versions.yml"                           , emit: versions
    path "*.args.txt"                             , emit: args

    script:
    def prefix = task.ext.prefix ?: "${meta.run}"
    def args   = task.ext.args   ?: 'nbases = 1e8'
    def seed   = task.ext.seed   ?: '100'
    
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    set.seed($seed)

    safe_learn <- function(files, out_rds) {
        err <- tryCatch({
            learnErrors(files, $args, multithread = $task.cpus, verbose = TRUE)
        }, error = function(e) {
            message("Error in learnErrors: ", e\$message)
            return(NULL) 
        })
        saveRDS(err, out_rds)
        return(err)
    }

    ${meta.single_end ? 
        """
        fnFs <- sort(list.files(".", pattern = ".filt.fastq.gz", full.names = TRUE))
        errF <- safe_learn(fnFs, "${prefix}.err.rds")
        
        if(!is.null(errF)){
            # Generate PDF
            pdf("${prefix}.err.pdf"); plotErrors(errF, nominalQ = TRUE); dev.off()
            # Generate SVG
            svg("${prefix}.err.svg"); plotErrors(errF, nominalQ = TRUE); dev.off()
            # Generate Convergence
            sink("${prefix}.err.convergence.txt"); dada2:::checkConvergence(errF); sink(NULL)
        }
        """ : 
        """
        fnFs <- sort(list.files(".", pattern = "_1.filt.fastq.gz", full.names = TRUE))
        fnRs <- sort(list.files(".", pattern = "_2.filt.fastq.gz", full.names = TRUE))
        
        errF <- safe_learn(fnFs, "${prefix}_1.err.rds")
        errR <- safe_learn(fnRs, "${prefix}_2.err.rds")
        
        if(!is.null(errF) && !is.null(errR)){
            pdf("${prefix}_1.err.pdf"); plotErrors(errF, nominalQ = TRUE); dev.off()
            svg("${prefix}_1.err.svg"); plotErrors(errF, nominalQ = TRUE); dev.off()
            
            pdf("${prefix}_2.err.pdf"); plotErrors(errR, nominalQ = TRUE); dev.off()
            svg("${prefix}_2.err.svg"); plotErrors(errR, nominalQ = TRUE); dev.off()
            
            sink("${prefix}_1.err.convergence.txt"); dada2:::checkConvergence(errF); sink(NULL)
            sink("${prefix}_2.err.convergence.txt"); dada2:::checkConvergence(errR); sink(NULL)
        }
        """
    }

    write.table('learnErrors\t$args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}
