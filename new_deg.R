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

#load tables into a list of data tables - "." should point to counts directory, e.g. "counts/."
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

	
# merged  merged
out <- res.merged[[1]][,c(1:2)]
lapply(res.merged,function(o) out<<-cbind(out,o[,c(3,7)]))
out <- cbind(out,res.merged[[1]][,8:16])
colnames(out)[3:16] <- c("FC_02793","P_02793","FC_F55","P_F55","FC_10170","P_10170","FC_MWT","P_MWT","FC_MOL","P_MOL","FC_MKO","P_MKO","FC_TJ","P_TJ")
write.table(out,"all.merged.csv",sep=",",quote=F,na="",row.names=F)

# sig all
write.table(subset(out,P_02793<=0.05&P_F55<=0.05&P_10170<=0.05&P_MWT<=0.05&P_MOL<=0.05&P_MKO<=0.05&P_TJ<=0.05),"all.sig.csv",sep=",",quote=F,na="",row.names=F)
	
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
	
# PCA 1 vs 2 plot
vst <- varianceStabilizingTransformation(dds,blind=F,fitType="local")
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH1"] <- "02780"
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH2"] <- "02793"
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH3"] <- "F55"
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH4"] <- "10170"
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH5"] <- "MWT"
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH6"] <- "MOL"
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH7"] <- "MKO"
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH8"] <- "TJ"
mypca <- prcomp(t(assay(vst)))
mypca$percentVar <- mypca$sdev^2/sum(mypca$sdev^2)
df <- t(data.frame(t(mypca$x)*mypca$percentVar))

pdf("quorn.pca.pdf",height=8,width=8)
plotOrd(df,vst@colData,design="condition",xlabel=PC1,ylabel=PC2, pointSize=3)
dev.off()
	
# MA plots	
pdf("MA_plots.pdf")
lapply(res.merged,function(obj) {
	with(obj,plot(log2FoldChange,log10(baseMean),pch=20, xlim=c(-6,6),bty="n",
		xlab=expression("Log"[2]*" Fold Change"),ylab=expression("Log"[10]*" Mean Expression")))
	with(subset(obj, padj<0.5 ), points(log2FoldChange, log10(baseMean), pch=20, col="#E69F00"))
	with(subset(obj, abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="#56B4E9"))
	with(subset(obj, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, log10(baseMean), pch=20, col="#009E73"))
})
dev.off()
