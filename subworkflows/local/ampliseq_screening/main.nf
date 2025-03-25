include { AMPLISEQ_SIMPLIFIED                                       } from '../ampliseq_simplified/main'


workflow AMPLISEQ_SCREENING {
    take:
    ch_samplesheet
    
    main:
    // generate sets of parameters based on input ranges
    GENERATE_PARAMS(params.trunclenf_range, params.trunclenr_range)

    // create a channel with parameters as input to ampliseq (simplified from nf-core)
    ch_params = GENERATE_PARAMS.out.params_csv
    .splitCsv(header: true, sep: ',')
    .map { row -> tuple(row.runID, row.trunclenf, row.trunclenr) }

    /*
    ch_samplesheet.combine(ch_params).map{ meta, reads, runID, trunclenf, trunclenr ->
        def new_meta = meta.clone()  // Clone meta to avoid mutating the original object
        new_meta.sample = meta.id
        new_meta.id = "run_${trunclenf}${trunclenr}.${meta.id}"  // Create new sample field with concatenated value
        tuple(new_meta + [runID: runID], reads, ['FW', trunclenf], ['RV', trunclenr])
    }.set{ ch_samplesheet_ }*/

    ch_samplesheet.combine(ch_params)
    //.map{meta, read1, read2, _, runID, trunclenf, trunclenr -> tuple([meta+[runID:runID], read1, read2, _], trunclenf, trunclenr)}
    //.set{ ch_samplesheet_w_params }
    .map { meta, read1, read2, _, runID, trunclenf, trunclenr -> 
    def new_meta = meta + [ sample: meta.id, id: "${meta.id}.${runID}", runID: runID, run: runID, trunclenf: trunclenf, trunclenr: trunclenr]  
    tuple(new_meta, read1, read2, _)}
    .set{ ch_samplesheet_w_params }
    

    AMPLISEQ_SIMPLIFIED(ch_samplesheet_w_params)
    /* AMPLISEQ_SIMPLIFIED(
         ch_samplesheet_w_params.map{ sample, trunclenf, trunclenr -> sample},
         ch_samplesheet_w_params.map{ sample, trunclenf, trunclenr -> trunclenf},
         ch_samplesheet_w_params.map{ sample, trunclenf, trunclenr -> trunclenr}) 
    */
}

process GENERATE_PARAMS {

    def out_path = file(params.outdir).toString() + '/generate_params/'

    publishDir "$out_path", mode: 'copy'

    input:
    val trunclenf_range
    val trunclenr_range

    output:
    path "summary_params_settings.csv", emit: params_csv

    script:
    """
    generate_params.py -f $trunclenf_range -r $trunclenr_range -o .
    """
}

