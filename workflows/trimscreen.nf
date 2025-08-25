/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samplesheetToList                                         } from 'plugin/nf-schema'
include { DETAXIZER_SIMPLIFIED                                      } from '../subworkflows/local/detaxizer_simplified/main.nf'
include { AMPLISEQ_SCREENING                                        } from '../subworkflows/local/ampliseq_screening/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// This is for new TRIMSCREEN after making detaxizer optional
workflow TRIMSCREEN {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    if (params.enable_filter) {

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
