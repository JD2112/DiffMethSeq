// enable dsl2
nextflow.enable.dsl=2

// DNA Methylation data analysis from Illumina, Nanopore, Pacific Biosciences and Twist methylation panel.

// Author: Jyotirmoy Das

// Date Created: January 17, 2024
// Date Last Modified: 

params.formatted_designfile = false
params.compare_str = false
params.cov_dir = false
params.bam_dir = false
params.reads_count_input1 = false
params.reads_count_input2 = false
params.phenoGroup = false

log.info """\
 DNAm-NF PIPELINE
 ===================================
 reads          : ${params.formatted_designfile}
 comparefile    : ${params.compare_str} 
 coverage_dir   : ${params.cov_dir}
 bam_dir        : ${params.bam_dir}
 pheno          : ${params.phenoGroup}
 reads_count1   : ${params.reads_count_input1}
 reads_count2   : ${params.reads_count_input2}
 outdir         : ${params.outdir}
 """

include { EdgeR } from './modules/edger/main'
include { MethylKit } from './modules/methylkit/main'
//include { DESeq2 } from './modules/deseq2/main'
// include { BISMARK_ALIGNMENT } from '../modules/bismark/alignment/main'
// include { BISMARK_DEDUP } from '../modules/bismark/deduplicate/main'
// include { BISMARK_METHYLATIONEXTRACTOR } from '../modules/bismark/methylextractor/main'
// include { BISMARK_REPORT } from '../modules/bismark/report/main'


workflow {        
        design_file = Channel.fromPath( params.formatted_designfile ).collect()
        comp_file = Channel.fromPath( params.compare_str ).collect()        
        cover_file = Channel.fromPath( params.cov_dir)
        bamDir = Channel.fromPath( params.bam_dir )
        phenoGr = Channel.fromPath( params.phenoGroup ).collect() 
        //reads1 = Channel.fromPath( params.reads_count_input1 ).collect()
        //reads2 = Channel.fromPath( params.reads_count_input2 ).collect()
   
        EdgeR(design_file, comp_file, cover_file)         
        //EdgeR(design_file, comp_file, reads1, reads2)

        //DESeq2(design_file, comp_file, reads1, reads2)
        MethylKit(design_file, comp_file, bamDir, phenoGr)
  
}