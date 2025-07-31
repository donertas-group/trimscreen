include { AMPLISEQ_SIMPLIFIED                                       } from '../ampliseq_simplified/main'
include { GENERATE_PARAMS                                           } from '../../../modules/local/generate_params'

workflow AMPLISEQ_SCREENING {
    take:
    ch_samplesheet
    
    main:
    // generate sets of parameters based on input ranges
    
    FW_primer_len = params.FW_primer ? params.FW_primer.size() : 0
    RV_primer_len = params.RV_primer ? params.RV_primer.size() : 0

    GENERATE_PARAMS (
        params.marker_size_min, 
        params.minimum_overlap, 
        params.step_size, 
        params.read_length,
        FW_primer_len,
        RV_primer_len
    )

    // create a channel with parameters as input to ampliseq (simplified from nf-core)
    ch_params = GENERATE_PARAMS.out.params_csv
    .splitCsv(header: true, sep: ',')
    .map { row -> tuple(row.runID, row.trunclenf, row.trunclenr) }

    AMPLISEQ_SIMPLIFIED(ch_samplesheet, ch_params)
}
