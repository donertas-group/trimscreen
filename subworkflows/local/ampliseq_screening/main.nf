workflow AMPLISEQ_SCREENING {
    
    GENERATE_PARAMS(params.trunclenf_range, params.trunclenr_range)

    // Create a channel from your parameter file
    //ch_generate_params = GENERATE_PARAMS.out.params_csv
   // ch_params = file("${params.outdir}/generate_params/summary_params_settings.csv")
   // .splitCsv(header: true, sep: ',')
   // .map { row -> tuple(row.runID, row.trunclenf, row.trunclenr) }

/*        .splitText()
        .map { line ->
            def fields = line.trim().split(',')  // Adjust the split method based on your file format
            return [
                runID: fields[0],
                trunclenf: fields[1],
                trunclenr: fields[2]
            ]*/
        


   // AMPLISEQ_SIMPLIFIED(ch_params)
}

// Define your subworkflow
workflow AMPLISEQ_SIMPLIFIED {
    take:
    ch_params

    main:
    // Your subworkflow logic here
    PROCESS_1(ch_params)
    // PROCESS_2(PROCESS_1.out)
    // ... more processes

    //emit:
    // Define your outputs here
}

process PROCESS_1 {
    input:
    tuple val(runID), val(trunclenf), val(trunclenr) // ... more parameters

    output:
    // Define your outputs

    script:
    """
    echo "$runID $trunclenf $trunclenr"
    """
}

process GENERATE_PARAMS {

    def out_path = file(params.outdir).toString() + '/generate_params/'

    publishDir "$out_path", mode: 'copy'
    //publishDir "${params.outdir}/generate_params", mode: 'copy'

    input:
    val trunclenf_range
    val trunclenr_range

    output:
    //path "${out_path}/summary_params_settings.csv", emit: params_csv
    path "summary_params_settings.csv"

    script:
    """
    generate_params.py -f $trunclenf_range -r $trunclenr_range -o .
    """
}

// ... more process definitions
