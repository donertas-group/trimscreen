/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                                                    } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                                          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                    } from '../subworkflows/local/utils_nfcore_trimscreen_pipeline'
include { getGenomeAttribute                                        } from '../subworkflows/local/utils_nfcore_trimscreen_pipeline'
include { GENERATE_DOWNSTREAM_SAMPLESHEETS                          } from '../subworkflows/local/generate_downstream_samplesheets/main.nf'// detaxizer has original subworkflow
include { DETAXIZER_SIMPLIFIED                                      } from '../subworkflows/local/detaxizer_simplified/main.nf'
include { AMPLISEQ_SCREENING                                        } from '../subworkflows/local/ampliseq_screening/main.nf'
include { AMPLISEQ_SIMPLIFIED                                       } from '../subworkflows/local/ampliseq_simplified/main.nf'

include { BBMAP_BBDUK                                               } from '../modules/nf-core/bbmap/bbduk/main'
include { BLAST_BLASTN                                              } from '../modules/nf-core/blast/blastn/main'
include { BLAST_MAKEBLASTDB                                         } from '../modules/nf-core/blast/makeblastdb/main'

include { ISOLATE_BBDUK_IDS                                         } from '../modules/local/detaxizer/isolate_bbduk_ids'
include { MERGE_IDS                                                 } from '../modules/local/detaxizer/merge_ids'
include { RENAME_FASTQ_HEADERS_PRE                                  } from '../modules/local/detaxizer/rename_fastq_headers_pre'
include { PREPARE_FASTA4BLASTN                                      } from '../modules/local/detaxizer/prepare_fasta4blastn'
include { FILTER_BLASTN_IDENTCOV                                    } from '../modules/local/detaxizer/filter_blastn_identcov'
include { FILTER                                                    } from '../modules/local/detaxizer/filter'
include { RENAME_FASTQ_HEADERS_AFTER                                } from '../modules/local/detaxizer/rename_fastq_headers_after'
include { SUMMARY_CLASSIFICATION                                    } from '../modules/local/detaxizer/summary_classification'
include { SUMMARY_BLASTN                                            } from '../modules/local/detaxizer/summary_blastn'
include { SUMMARIZER                                                } from '../modules/local/detaxizer/summarizer'


//include { GENERATE_PARAMS                                           } from '../modules/local/generate_params'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// specify the ch_fasta_blastn channel if it is not provided via --fasta_blastn
def ch_fasta_blastn = Channel.empty()

if ( !params.fasta_blastn && params.validation_blastn ) {
    ch_fasta_blastn = Channel.fromPath(getGenomeAttribute('fasta'))
} else if ( params.validation_blastn ){
    // If params.fasta_blastn is there, use it for the creation of the blastn database
    ch_fasta_blastn = Channel.fromPath(params.fasta_blastn)
}

// specify the ch_fasta_bbduk channel if it is not provided via --fasta_bbduk

def ch_fasta_bbduk = Channel.empty()

if ( !params.fasta_bbduk && params.classification_bbduk ) {
    ch_fasta_bbduk = Channel.fromPath(getGenomeAttribute('fasta'))
} else if ( params.classification_bbduk ){
    // If params.fasta_bbduk is there, use it for the creation of the blastn database
    ch_fasta_bbduk = Channel.fromPath(params.fasta_bbduk)
}


// This is for new TRIMSCREEN after making detaxizer optional
workflow TRIMSCREEN {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    if (params.run_host_removal) {

        DETAXIZER_SIMPLIFIED(ch_samplesheet)
        downstream_samplesheet = DETAXIZER_SIMPLIFIED.out.ch_samplesheet

        AMPLISEQ_SCREENING (downstream_samplesheet)

        multiqc_report = AMPLISEQ_SCREENING.out.multiqc_report

    } else {

        AMPLISEQ_SCREENING (ch_samplesheet)

        //multiqc_report = AMPLISEQ_SCREENING.out.multiqc_report
    }
    









    //emit:
    //multiqc_report
}




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
