#!/bin/Rscript
## Rscript edgeR.R <designfile> <compare_str> aligner eg. Rscript edgeR.R designfile_single.txt T_vs_N star
### designfile: Sample_id, Input_filename, IP_filename, group_id
### compare_str: Compairision design (eg: A_vs_B)
library("edgeR")
library(readr)
library(org.Hs.eg.db)

args<-commandArgs(T) 
#tabledata <- args[1]
#compare_str <- as.character(args[2])
targets <- read_csv(args[1])
compare_str <- read_file(args[2])


#targets = read.csv(tabledata)
Sample <- targets$Sample_Name

files <- targets$File
#files <- args[3]

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

# read bismark data to deg list
tt <- readBismark2DGE(files, sample.names = Sample, readr = TRUE, verbose =TRUE)

# cleaning of unwanted chromosomes
keep<- rep(TRUE,nrow(tt))
Chr<- as.character(tt$genes$Chr)
keep[grep("random",Chr) ]<-FALSE
keep[grep("chrUn",Chr) ]<-FALSE
#keep[Chr=="chrY"]<-FALSE ##need to ask
keep[Chr=="chrM"]<-FALSE
#keep[Chr=="chrX"]<-FALSE ##need to ask
tt1<-tt[keep,,keep.lib.sizes=FALSE]

tt1_data_rem <- cbind(tt1$counts, tt1$genes)

# Assign chromosome names
ChrNames<- paste0("chr",c(1:22, "X", "Y"))
tt1$genes$Chr<- factor(tt1$genes$Chr,levels=ChrNames)
o<- order(tt1$genes$Chr, tt1$genes$Locus)
tt1<-tt1[o,]
tt1_data_rem1 <- cbind(tt1$counts, tt1$genes)

# annotation 
# TSS<- nearestTSS(tt1$genes$Chr, tt1$genes$Locus, species="Hs")

# tt1$genes$EntrezID<-TSS$gene_id
# tt1$genes$Symbol<-TSS$symbol
# tt1$genes$Strand<-TSS$strand
# tt1$genes$Distance<-TSS$distance
# tt1$genes$Width<-TSS$width

# Filtering and normalizing
Methylation <- gl(2,1,ncol(tt1), labels=c("Me","Un"))
Me <- tt1$counts[, Methylation == "Me"]
Un <- tt1$counts[, Methylation=="Un"]
Coverage <- Me + Un

## keeping at least 3 cpgs per sample
keep<- rowSums(Coverage >= 3) == 8
HasBoth <- rowSums(Me) > 0 & rowSums(Un) > 0

y <-tt1[keep & HasBoth,,keep.lib.sizes=FALSE]
TotalLibSize <- 0.5 * y$samples$lib.size[Methylation=="Me"] + 0.5 * y$samples$lib.size[Methylation=="Un"]
y$samples$lib.size <- rep(TotalLibSize, each=2)

# Data Exploration
Me <- y$counts[, Methylation=="Me"]
Un <- y$counts[, Methylation=="Un"]
M <- log2(Me + 2) - log2(Un + 2)
colnames(M) <- targets$Sample_Name

y1 <- cbind(y$counts, y$genes)
My1 <- cbind(M, y1)

# design the model.matrix
designSL<- model.matrix(~ 0 + factor(Group), data=targets) # check Sample_sheet.csv if matrix shows different groups.
colnames(designSL) <- c(group_id_1, group_id_2)

design<- modelMatrixMeth(designSL)

## by passing with estimateGLMCommonDisp
yy <- estimateGLMCommonDisp(y, design, verbose = TRUE)
# calculate toptags
fit<- glmFit(yy, design)


#contr<- makeContrasts(DvsH=group_id_2 - group_id_1,levels=design)
#lrt<- glmLRT(fit, contrast=contr)
lrt<- glmLRT(fit, coef=2)
topTags(lrt)

#topME<- topTags(lrt,n=Inf,p=0.05)$table
# save table
lrt$table$padj <- p.adjust(lrt$table$PValue, "BH")
lrt.res <- lrt$table[order(lrt$table$padj),]

output_name <- paste0("edgeR_group_",group_id_1, "_",group_id_2)
write.csv(lrt.res, file = paste0(output_name, ".csv"), quote = FALSE)
#write.table(topME, "/mnt/SD3/CFFMHS-KW-/edgeR_results/DMClist_DvsH.txt", sep = "\t", quote = FALSE)