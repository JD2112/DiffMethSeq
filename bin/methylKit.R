#!/bin/Rscript
## Author: Jyotirmoy Das
## Date Created: January 17, 2024
## Date Last Modified: February 11, 2024
## Rscript edgeR.R <designfile> <compare_str> aligner eg. Rscript edgeR.R designfile_single.txt T_vs_N star
### designfile: Sample_id, Input_filename, IP_filename, group_id
### compare_str: Compairision design (eg: A_vs_B)

suppressPackageStartupMessages(library(methylKit))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(genomation))

args<-commandArgs(T) 

targets <- read_csv(args[1])
compare_str <- read_file(args[2])
bam_dir = args[3]
phenoGroup <- noquote(read_file(args[4]))

## add files
#pheno <- read_file(phenoGroup)

# Running edgeR by compare.file
## while there are only 2 groups, running edgeR without compare.file
if(length(unique(targets$Group)) < 2){
  stop( "The count of Group is less than two, please check your designfile.")
}else if( compare_str == "two_groups" ){
  # Running edgeR without compare_str beacause of only two groups
  ## Combine expression matrix
  group_id_1 <- unique(targets$Group)[1]
  group_id_2 <- unique(targets$Group)[2]
}else{
  # Running edgeR with compare_str 
  ## Combine expression matrix
  group_id_1 <- strsplit(as.character(compare_str), "_vs_")[[1]][1]
  group_id_2 <- strsplit(as.character(compare_str), "_vs_")[[1]][2]
}

# using regrex for files

#cpgFiles <- dir("./data", full.name=TRUE, pattern = "*.bismark.cov.gz")
cpgFiles <- dir(bam_dir, full.name=TRUE, pattern = "*.bam")
#cpgFiles_bam <- list.files(paste0(bamfile.list, cpgFiles))

sampleNames <- basename(unname(sapply(cpgFiles, function(x) str_extract(x, ".+?(?=_)"))))

myObj = processBismarkAln(
    location = as.list(cpgFiles),
    sample.id = as.list(sampleNames),
    assembly = "hg38",
    #treatment = phenoGroup,
    treatment = c(0,1,0,1,0,1,0,1),
    mincov = 5,
    minqual = 20
)

# filtering the data
filtered.myObj = filterByCoverage(myObj, 
                                    lo.count = 10, # default is 10
                                    lo.perc = NULL, 
                                    hi.count = NULL, hi.perc = 99.9)

# Comparative analysis
## Normalization
myobj.filt.norm <- normalizeCoverage(filtered.myObj, method = "median")

## merge data
meth1 <- unite(myobj.filt.norm, destrand=FALSE)

# get percent methylation matrix
pm=percMethylation(meth1)

# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)

# keep only CpG with standard deviations larger than 2%
meth1 <- meth1[sds > 2]

#calculate diff meth bases

myDiff.filt <- calculateDiffMeth(meth1,
                            overdispersion = "MN",
                            adjust="BH",
                            mc.cores = 12)


# get hyper methylated bases - 5%
myDiff5p.hyper.filt=getMethylDiff(myDiff.filt,difference=5,qvalue=0.01,type="hyper")
# get hypo methylated bases
myDiff5p.hypo.filt=getMethylDiff(myDiff.filt,difference=5,qvalue=0.01,type="hypo")
# get all differentially methylated bases
myDiff5p.filt=getMethylDiff(myDiff.filt,difference=5,qvalue=0.01)

# remove rows with unknown chromosome
myDiff5p.filt.clean <- myDiff5p.filt[!grepl("_random", myDiff5p.filt$chr),]

# download refseq file on the fly
url <- "https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz/download"
tmpfile <- tempfile(tmpdir=getwd())
file.create(tmpfile)
download.file(url, tmpfile)
file.rename(tmpfile, "hg38_RefSeq.bed.gz")
gene.obj=readTranscriptFeatures("hg38_RefSeq.bed.gz")

# annotate filtered results
myDiff5p.filt.annot <- annotateWithGeneParts(as(myDiff5p.filt.clean,"GRanges"),gene.obj)

test.Annot.methDiff <- cbind(myDiff5p.filt.clean, getAssociationWithTSS(myDiff5p.filt.annot))

# Get the annotation
anno <- AnnotationDbi::select(org.Hs.eg.db,keys=gsub("\\..*","", test.Annot.methDiff$feature.name),
              columns=c("SYMBOL","GENENAME"),
              keytype="REFSEQ")

# remove the refseq extra . from names and get the data
DMC.final.annot <- cbind(test.Annot.methDiff, anno)

output_name <- paste0("methylkit_group_",group_id_1, "_",group_id_2)
write.csv(DMC.final.annot, file = paste0(output_name,".csv"), quote = FALSE)