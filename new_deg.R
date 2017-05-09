#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(ggplot2)
library(Biostrings)
library(devtools)
load_all("myfunctions") # this contains various R scripts for plotting graphs
library(data.table)
library(naturalsort)
library(tibble)

#===============================================================================
#       Load features counts data 
#===============================================================================

#load tables into a list of data tables
qq <- lapply(list.files(".",".*.txt$",full.names=T,recursive=F),function(x) fread(x))

# rename the sample columns (7th column in a feature counts table, saved as the path to the BAM file)
# in the below I'm saving the 8th ([[1]][8]) path depth (which was the informative folder name containg the BAM file)
invisible(lapply(seq(1:length(qq)), function(i) colnames(qq[[i]])[7]<<-strsplit(colnames(qq[[i]])[7],"\\/")[[1]][8]))

#merge the list of data tables into a single data table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)

#output "countData"
write.table(m[,c(1,7:(ncol(m))),with=F],"countData",sep="\t",na="",quote=F,row.names=F) 
#output gene details
write.table(m[,1:6,with=F],"genes.txt",sep="\t",quote=F,row.names=F) 

#==========================================================================================
#       Read pre-prepared colData,  countData and annotations
##=========================================================================================

colData <- read.table("colData",header=T,sep="\t")
# colData$condition <- rep(c("02780","02793","F55","10170","MWT","MOL","MKO","TJ"),3) # need to test this - will set columns to numbers 
countData <- read.table("countData",sep="\t",header=T,row.names=1) # produced above, could just subset the data table countData <- m[,c(1,7:length(m),with=F]	
countData <- countData[,colData$SampleID] # reorder countData columns to same order as colData rows

annotations <- fread("WT_annotation.tsv")
	
	
#===============================================================================
#       DESeq2 analysis
#		Set alpha to the required significance level. This also effects how
#		DESeq calculated FDR - setting to 0.05 and then extracting results with a
#		significance below 0.01 will give slightly different results form setting
#		alpha to 0.01
#================================================================================

dds <- 	DESeqDataSetFromMatrix(countData,colData,~1) 
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds))
dds$groupby <- paste(dds$condition,dds$sample,sep="_")
dds <- collapseReplicates(dds,groupby=dds$groupby)
design=~condition
design(dds) <- design # could just replace the ~1 in the first step with the design, if you really wanted to...
dds <- DESeq(dds,parallel=T)

alpha <- 0.01

# calculate the differences - uses the "levels" of the condition factor as the third term for the contrast
# res is a list object containing the DESeq results objects for each contrast
# contrast=c("condition","RH1","RH2") etc. (the below just runs through all of the different sample types (excluding RH1))
res <- lapply(seq(2,8), function(i) results(dds,alpha=alpha,contrast=c("condition","RH1",levels(dds$condition)[i])))
names(res) <- c("02793","F55","10170","MWT","MOL","MKO","TJ")
		# RH2,   RH3,  RH4,    RH5,  RH6,  RH7,  RH8
		#lapply(res,function(x) strsplit(x@elementMetadata@listData$description[2]," ")[[1]][8])
	
# merge results with annotations
res.merged <- lapply(res,function(x) left_join(rownames_to_column(as.data.frame(x)),annotations,by=c("rowname"="query_id")))	
	
# get, then order the significant results
sig.res <- lapply(res.merged, function(x) subset(x,padj<=alpha))
sig.res <- lapply(sig.res,function(x) x[order(x$padj),])

# sig in all conditions	
test <- sig.res[[7]]	
newtest <- lapply(seq(1:6), function(i) test<<-inner_join(sig.res[[i]],test,by="rowname"))
newtest[[6]]$rowname
all.sig <- lapply(sig.res,function(o) o[which(o$rowname%in%newtest[[6]]$rowname),])

# write tables of results, and significant results
lapply(seq(1:7),function(x) {
	write.table(res.merged[[x]],paste(names(res.merged)[x],"merged.txt",sep="."),quote=F,na="",row.names=F,sep="\t")
	write.table(sig.res[[x]],paste(names(sig.res)[x],"sig.merged.txt",sep="."),quote=F,na="",row.names=F,sep="\t")
	write.table(all.sig[[x]],paste(names(all.sig)[x],"all.sig.merged.txt",sep="."),quote=F,na="",row.names=F,sep="\t")
})	
	
	
#===============================================================================
#       FPKM
#===============================================================================

rowRanges(dds) <- GRangesList(apply(m,1,function(x) GRanges(x[[1]],IRanges(1,as.numeric(x[[6]])),"+")))
myfpkm <- data.table(GeneID=m[,1],length=m[,6],fpkm(dds,robust=T))
write.table(myfpkm,"fpkm.txt",quote=F,na="",sep="\t")
	
#===============================================================================
#       Heirachical clustering
#===============================================================================

clus <- function(X,clusters=10,m=1,name="hclust.pdf") {
	if (m==1) {d <- dist(X, method = "manhattan")}
	else if (m==2) {d <- dist(X, method = "euclidean")}
	else if (m==3) {d <- dist(X, method = "maximum")}
	else if (m==4) {d <- dist(X, method = "canberra")}
	else if (m==5) {d <- dist(X, method = "binary")}
	else d <- {dist(X, method = "minkowski")}
	hc <- hclust(d, method="ward")
	groups <- cutree(hc, k=clusters) # cut tree into n clusters
	pdf(name,height=8,width=8)
	plot(hc)
	rect.hclust(hc,k=clusters)
	dev.off()
	return(list(hc,groups,d))
}

#===============================================================================
#       Graphs
#===============================================================================

rld <- varianceStabilizingTransformation(dds,blind=F,fitType="local")
rld$label <- dds$sample
rld$condition <- c("02780","02793","F55","10170","MWT","MOL","MKO","TJ","02780","02793","F55","10170","MWT","MOL","MKO","TJ","02780","02793","F55","10170","MWT","MOL","MKO","TJ")
pdf("quorn.pca_2.pdf",height=10,width=10)
plotPCAWithLabels(rld)
dev.off()
