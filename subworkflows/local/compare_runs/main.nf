
workflow COMPARE_RUNS {
    take:
    runs_summary 
    runs_asv_table
    runs_asv_tax

    main:
    runs_summary.view()
    runs_asv_table.view()
    runs_asv_tax.view()


   // emit:

}
