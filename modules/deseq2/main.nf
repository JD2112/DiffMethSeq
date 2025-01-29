process DESeq2{
    label 'analysis'
    tag "$compare_str"

    publishDir "${params.outdir}/DESeq2", mode: 'link', overwrite: true

    input:
    // file reads_count_input from htseq_count_input.collect()
    // file formatted_designfile from formatted_designfile.collect()
    // val compare_str from compareLines_for_DESeq2    
     path(formatted_designfile)
     path(compare_str)
     path(reads_count_input1)
     path(reads_count_input2)

    output:
    file "DESeq2*.csv" 
    
    // when:
    // !params.skip_deseq2 && !params.skip_expression && params.comparefile
    
    script:
    // println LikeletUtils.print_purple("Differential expression analysis performed by DESeq2 ($compare_str)")
    """
    Rscript $baseDir/bin/DESeq2.R ${formatted_designfile} ${compare_str} ${reads_count_input1} ${reads_count_input2}
    """ 
}