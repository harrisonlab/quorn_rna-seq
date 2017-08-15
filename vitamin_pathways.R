
dds <- DESeqDataSetFromMatrix(countData,colData,~1)

dds$groupby <- paste(dds$condition,dds$sample,sep="_")

dds <- collapseReplicates(dds,groupby=dds$groupby)

sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds))

design=~condition

design(dds) <- design 

dds <- DESeq(dds,parallel=T)

alpha <- 0.05

res <- lapply(seq(2,8), function(i) results(dds,alpha=alpha,contrast=c("condition",levels(dds$condition)[i],"RH1"),parallel=T))

res[[8]] <- results(dds,alpha=alpha,contrast=c("condition","RH5","RH7"),parallel=T)

pathways <- fread("../vitamin_pathways/pathways.txt",sep="\t")

res.merged <- lapply(res,function(x) left_join(pathways,rownames_to_column(as.data.frame(x)),by=c("Gene_ID"="rowname")))

out <- res.merged[[1]][,c(1:7)]

invisible(lapply(res.merged,function(o) out<<-cbind(out,o[,c(8,12)])))

colnames(out)[8:23] <- c("FC_02793","P_02793","FC_F55","P_F55","FC_10170","P_10170","FC_MWT","P_MWT","FC_MOL","P_MOL","FC_MKO","P_MKO","FC_TJ","P_TJ","FC_MWT_MKO","P_MWT_MKO")

write.table(out,"../vitamin_pathways/pathways.merged.tsv",sep="\t",quote=F,na="",row.names=F)
