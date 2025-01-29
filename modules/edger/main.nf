//Rscript $baseDir/bin/edgeR.R $formatted_designfile $compare_str $params.aligners
process EdgeR{
    label 'analysis'
    tag "$compare_str"
    publishDir "${params.outdir}/edgeR", mode: 'link', overwrite: true

    input:
    // file reads_count_input from htseq_count_input.collect()
    // file formatted_designfile from formatted_designfile.collect()
    // val compare_str from compareLines_for_edgeR
    path(formatted_designfile)
    path(compare_str)
    path(cov_dir)
    //path(reads_count_input1)
    //path(reads_count_input2) 
    
 
    output:
    file "edgeR*.csv" 
    
    // when:
    // !params.skip_edger && !params.skip_expression && params.comparefile

    script:
    // println LikeletUtils.print_purple("Differential expression analysis performed by EdgeR ($compare_str)")
    //Rscript $baseDir/bin/edgeR.R ${formatted_designfile} ${compare_str} ${reads_count_input1} ${reads_count_input2}
    """
    Rscript $baseDir/bin/edgeR.R ${formatted_designfile} ${compare_str} ${cov_dir}
    """ 
}