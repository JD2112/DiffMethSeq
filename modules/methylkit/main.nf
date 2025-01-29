// DNA Methylation data analysis from Illumina, Nanopore, Pacific Biosciences and Twist methylation panel.

// Author: Jyotirmoy Das

// Date Created: January 17, 2024
// Date Last Modified: February 11, 2024

process MethylKit {
    label 'analysis'
    tag "$compare_str"
    publishDir "${params.outdir}/MethylKit", mode: 'link', overwrite: true

    input:
    path(formatted_designfile)
    path(compare_str)
    path(bam_dir)
    path(phenoGroup)    
 
    output:
    file "methylkit*.csv" 
    
    script:

    """
    Rscript $baseDir/bin/methylKit.R ${formatted_designfile} ${compare_str} ${bam_dir} ${phenoGroup}
    """ 
}