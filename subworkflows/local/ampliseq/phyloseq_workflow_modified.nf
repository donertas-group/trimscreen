/*
 * Create phyloseq objects
 */

include { PHYLOSEQ                                } from '../../../modules/local/ampliseq/phyloseq_modified'
include { PHYLOSEQ_INASV                          } from '../../../modules/local/ampliseq/phyloseq_inasv'

workflow PHYLOSEQ_WORKFLOW {
    take:
    ch_tax
    ch_tsv
    ch_meta
    ch_tree
    run_qiime2

    main:
    if ( params.metadata ) {
        ch_phyloseq_inmeta = ch_meta.first() // The .first() is to make sure it's a value channel
    } else {
        ch_phyloseq_inmeta = []
    }

    if ( params.pplace_tree ) {
        ch_phyloseq_intree = ch_tree.map { it = it[1] }.first()
    } else {
        ch_phyloseq_intree = []
    }

    if ( run_qiime2 ) {
        if ( params.exclude_taxa != "none" || params.min_frequency != 1 || params.min_samples != 1 ) {
            ch_phyloseq_inasv = PHYLOSEQ_INASV ( ch_tsv ).tsv
        } else {
            ch_phyloseq_inasv = ch_tsv
        }
    } else {
        ch_phyloseq_inasv = ch_tsv
    }

    PHYLOSEQ ( ch_tax.combine(ch_phyloseq_inasv, by:0).map { it.flatten() }, ch_phyloseq_inmeta, ch_phyloseq_intree )// flatten it as meta is included in ch_phyloseq_inasv

    emit:
    rds      = PHYLOSEQ.out.rds
    versions = PHYLOSEQ.out.versions
}
