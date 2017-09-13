#===============================================================================
#       Load libraries
#===============================================================================

library(DEXSeq)
library("BiocParallel")
BPPARAM=MulticoreParam(12)
register(BPPARAM)
library(data.table)

#===============================================================================
#       Load featureCounts data 
#===============================================================================

# load tables into a list of data tables - "." should point to counts directory, e.g. "counts/."
qq <- lapply(list.files(".",".*.txt$",full.names=T,recursive=F),function(x) fread(x,skip=1)) 

# rename the sample columns (7th column in a feature counts table, saved as the path to the BAM file)
# in the below this removes everything after the first dot in the 7th column
invisible(lapply(seq(1:length(qq)), function(i) colnames(qq[[i]])[7]<<-sub("\\..*","",colnames(qq[[i]])[7])))

# merge the list of data tables into a single data table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)

# output "countData"
write.table(m[,c(1,7:(ncol(m))),with=F],"countData",sep="\t",na="",quote=F,row.names=F) 
	    
# output gene details
write.table(m[,1:6,with=F],"genes.txt",sep="\t",quote=F,row.names=F) 

#==========================================================================================
#      Read in data 
##=========================================================================================

# Sample data 	    
colData <- read.table("colData",sep="\t",row.names=1,header=T)

# extract counts from m or load them from a file (countData will need to converted to a data frame before calling DEXSeqDataSet)	    
countData <- m[order(Chr,Start,Geneid),c(1,7:length(m)),with=F]r

# reorder countData columns to same order as colData rows 
countData <- countData[,c("Geneid",row.names(colData)),with=F]

# get the gene data	    
geneData <-  m[order(Chr,Start,Geneid),c(1,2:6),with=F]
	    
# add an exon column to geneData by counting occurance of each gene (ordered by position)
geneData[, Exonid := paste0("exon_",seq_len(.N)), by = Geneid]

# set GRanges for each exon (optional - but set featureRanges to NULL otherwise)	    
featureRanges <- GRanges(geneData$Chr,IRanges(geneData$Start,as.numeric(geneData$End)),geneData$Strand)	    

# retrieve transcipts (optional - but set transcripts to NULL otherwise)	    
features <- fread("grep exon featureCounts.gtf")	  
# get the gene id
features$Geneid <- gsub("\"|.*gene_id \"","",features$V9)	    
# get transcript id
features$Transcriptid <- gsub("\";.*|transcripts \"","",features$V9)	    
# left join geneData and features
geneData[features[,Geneid,Transcriptid],Transcriptid:=i.Transcriptid,on="Geneid"]

transcripts <- geneData$Transcriptid    

#==========================================================================================
#      DEXSeq analysis simple (one factor)
##=========================================================================================

##### NOTE #####
# estimateDispersions and testForDEU can both be parallelsed by adding 
# BPPARAM=BPPARAM
# BUT parallisation did not work when I was tesing using MulticoreParam 		    
# This is unfortunate as dispersion estimates may take hours to calculate

# combine technical replicates	    
dds <- 	DESeqDataSetFromMatrix(countData[,-1],colData,~1)
	    
# add grouping factor to identify technical replicates where sample contains the replicate info   
dds$groupby <- paste(dds$condition,dds$sample,sep="_")

# sum (collapse) replicates 
dds <- collapseReplicates(dds,groupby=dds$groupby)
	    
# create DEXSeq object     
dxd <- DEXSeqDataSet(assay(dds),
		     as.data.frame(colData(dds)),
		     featureID=geneData$Exonid,
		     groupID=geneData$Geneid,
		     featureRanges=featureRanges,
		     transcripts=transcripts
)
	    
# calculate model 
dxr = DEXSeq(dxd) # BPPARAM=BPPARAM   
	    
#==========================================================================================
#      Plots
##=========================================================================================	    
	    
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")	  

pdf("test.pdf",width=8)
plotDEXSeq( dxr, "g6103", displayTranscripts=F, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )	    

# plot everything (=< fdr)
DEXSeqHTML( dxr, FDR=0.05, color=cbbPalette,path=".")
	    
	    
