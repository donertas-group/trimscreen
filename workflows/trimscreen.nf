/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samplesheetToList                                         } from 'plugin/nf-schema'
include { getGenomeAttribute                                        } from '../subworkflows/local/utils_nfcore_trimscreen_pipeline'
include { DETAXIZER_SIMPLIFIED                                      } from '../subworkflows/local/detaxizer_simplified/main.nf'
include { AMPLISEQ_SCREENING                                        } from '../subworkflows/local/ampliseq_screening/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
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
*/

// This is for new TRIMSCREEN after making detaxizer optional
workflow TRIMSCREEN {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    if (!params.skip_host_removal) {

        DETAXIZER_SIMPLIFIED(ch_samplesheet)
        new_samplesheet = DETAXIZER_SIMPLIFIED.out.new_samplesheet

        ch_new_samplesheet = new_samplesheet.map { sheet ->
            def listified = samplesheetToList(sheet, "${projectDir}/assets/schema_input.json")
            return listified.collect { sample ->
                def meta = sample[0]  // Assuming the first element is the sample ID
                def files = sample[1..2]    // Assuming the next two elements are the file paths
                return tuple(meta, files)
            }
        }.flatten().collate(3)

        AMPLISEQ_SCREENING (ch_new_samplesheet)

       // multiqc_report = AMPLISEQ_SCREENING.out.multiqc_report

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
