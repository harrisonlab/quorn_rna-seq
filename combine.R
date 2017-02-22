#########################################################
#							#
#   concatenates a set of featureCounts files 		#
#   into a single countData file for use with DESeq	#
#							#
#########################################################

#install and load libraries		
require(data.table)

#load tables into a list of data tables
qq <- lapply(list.files(".",".*.txt$",full.names=T,recursive=T),function(x) fread(x))

#rename the sample columns 
lapply(qq,function(x) {names(x)[7]<-sub("\\/.*","",sub(".*\\/WTCHG","WTCHG",names(x)[7]))})

#merge the list of data tables into a single data table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)
	
#output "countData"
write.table(m[,c(1,7:(ncol(m))),with=F],"countData",sep="\t",na="",quote=F) 
#output gene details
write.table(m[,1:6,with=F],"genes.txt",sep="\t",quote=F,row.names=F) 
