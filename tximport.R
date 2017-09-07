library(tximport)
install.packages("readr")
library(DESeq2)

tx2gene <- read.table("trans2gene.txt",header=T,sep="\t")
txi <-tximport("quant.sf",type="salmon",tx2gene=tx2gene,txOut=T)
txi <- summarizeToGene(test,tx2gene)

dds<-DESeqDataSetFromTximport(txi,data.frame(S="S1",C="H"),design=~1)
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds)

## for DESeq < v1.12  (R 3.3 or greater)
DESeqDataSetFromTximport <- function(txi, colData, design, ...) 
{
  counts <- round(txi$counts)
  mode(counts) <- "integer"
  object <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design=design, ...)
  stopifnot(txi$countsFromAbundance %in% c("no","scaledTPM","lengthScaledTPM"))
  if (txi$countsFromAbundance %in% c("scaledTPM","lengthScaledTPM")) {
    message("using just counts from tximport")
  } else {
    message("using counts and average transcript lengths from tximport")
    lengths <- txi$length
    stopifnot(all(lengths > 0))
    dimnames(lengths) <- dimnames(object)
    assays(object)[["avgTxLength"]] <- lengths
  }
  return(object)
}   

nm <- assays(dds)[["avgTxLength"]]
nm <- nm / exp(rowMeans(log(nm)))
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds,normMatrix=nm))

