#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(data.table)

#==========================================================================================
#       Read data
##=========================================================================================

# column descrptions
colData <- read.table("colData",header=T,sep="\t")

# RNA-seq reads per sample
countData <- read.table("countData",sep="\t",header=T,row.names=1) 

# Identified vitamin pathway genes
pathways <- fread("../vitamin_pathways/pathways.txt",sep="\t")

#==========================================================================================
#       DGE analysis
##=========================================================================================

# create DES object
dds <- DESeqDataSetFromMatrix(countData,colData,~1)

# Create a grouping column 
dds$groupby <- paste(dds$condition,dds$sample,sep="_")

# Combine technical replicates by grouping column
dds <- collapseReplicates(dds,groupby=dds$groupby)

# Calcute library size normalisation factors
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds))

# DES model to calculate
design=~condition

# add model to DES object
design(dds) <- design 

# calculate model
dds <- DESeq(dds,parallel=T)

# BH adjust cut-off
alpha <- 0.05

# extract results from DES oject (fold changes relate to the first condition in the model)
res <- lapply(seq(2,8), function(i) results(dds,alpha=alpha,contrast=c("condition",levels(dds$condition)[i],"RH1"),parallel=T))

# extract results for MWT vs MKO ( as above - fold change relates to MWT ("RH5"))             
res[[8]] <- results(dds,alpha=alpha,contrast=c("condition","RH5","RH7"),parallel=T)

# join pathways to results
res.merged <- lapply(res,function(x) left_join(pathways,rownames_to_column(as.data.frame(x)),by=c("Gene_ID"="rowname")))

# create an output file  - 1:7 is gene name and descriptions
out <- res.merged[[1]][,c(1:7)]

# get the fold change and adj p for each contrast in the results list                     
invisible(lapply(res.merged,function(o) out<<-cbind(out,o[,c(8,12)])))

# rename the columns in the output file
colnames(out)[8:23] <- c("FC_02793","P_02793","FC_F55","P_F55","FC_10170","P_10170","FC_MWT","P_MWT","FC_MOL","P_MOL","FC_MKO","P_MKO","FC_TJ","P_TJ","FC_MWT_MKO","P_MWT_MKO")

# write the output file                  
write.table(out,"../vitamin_pathways/pathways.merged.tsv",sep="\t",quote=F,na="",row.names=F)
