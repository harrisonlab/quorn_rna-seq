library(sqldf)

res2 <- res.2
res3 <- res.3
res4 <- res.4
res5 <- res.5
res6 <- res.6
res7 <- res.7
res8 <- res.8

# perform multiple joins with sqldf is easier than with merge (proably could do the same with plyr which is faster than sqldf, but I don't know the syntax)
all_results <- sqldf("
	select 
		r2.log2FoldChange as r2_FC,
		r3.log2FoldChange as r3_FC,
		r4.log2FoldChange as r4_FC,
		r5.log2FoldChange as r5_FC,
		r6.log2FoldChange as r6_FC,
		r7.log2FoldChange as r7_FC,
		r8.log2FoldChange as r8_FC,
		r2.padj as r2_p,
		r3.padj as r3_p,
		r4.padj as r4_p,
		r5.padj as r5_p,
		r6.padj as r6_p,
		r7.padj as r7_p,
		r8.padj as r8_p
	from res2 r2
		join res3 r3 on r2.row_names=r3.row_names
		join res4 r4 on r2.row_names=r4.row_names
		join res5 r5 on r2.row_names=r5.row_names
		join res6 r6 on r2.row_names=r6.row_names
		join res7 r7 on r2.row_names=r7.row_names
		join res8 r8 on r2.row_names=r8.row_names "
   , row.names=T
)
rownames(all_results) <- rownames(res2) # I'm pretty certain this could have been added to the abpve (r2.row_names)

# find genes which are differentially expressed in at least one condition (FC>2, p<=.05)
all_results$test <- apply(all_results,1,function(x) {
	as.logical(sum(
		sqrt(x[1]^2)>1&x[8]<=0.05,
		sqrt(x[2]^2)>1&x[9]<=0.05,
		sqrt(x[3]^2)>1&x[10]<=0.05,
		sqrt(x[4]^2)>1&x[11]<=0.05,
		sqrt(x[5]^2)>1&x[12]<=0.05,
		sqrt(x[6]^2)>1&x[13]<=0.05,
		sqrt(x[7]^2)>1&x[14]<=0.05
	   )) 
	}
)

myfpkm <- myfpkm[,c(1,9,17,2,10,18,3,11,19,4,12,20,5,13,21,6,14,22,7,15,23,8,16,24)]


## calculate pearson corr matrix and return the mean value of the lower left triangle (excluding corr to self)
p_corr <- function(X) {
	t1 <- suppressWarnings(cor(t(X),use="all.obs",method="pearson"))
	t1[!lower.tri(t1)] <- NA
	t1 <- sqrt(t1^2)
	sum(t1,na.rm=T)/(0.5*(nrow(X)^2-nrow(X)))
}

## 1. get random set of three genes, 2. find average correlation 3. replicate 1000 times 4. get 95th quantile
## not certain whether to use a cutoff per experiment or a global value - global value is going to be a lot easier
cut_off_1 <-quantile(#4
	replicate(#3
		1000,
		p_corr( #2 
			myfpkm[sample(1:length(myfpkm[,1]),3, replace=FALSE),1:3] #1 
			)
		)
	,.95
)
cut_off_2 <- quantile(replicate(1000,p_corr(myfpkm[sample(1:length(myfpkm[,1]),3, replace=FALSE),4:6])),.95)
cut_off_3 <- quantile(replicate(1000,p_corr(myfpkm[sample(1:length(myfpkm[,1]),3, replace=FALSE),7:9])),.95)
cut_off_4 <- quantile(replicate(1000,p_corr(myfpkm[sample(1:length(myfpkm[,1]),3, replace=FALSE),10:12])),.95)
cut_off_5 <- quantile(replicate(1000,p_corr(myfpkm[sample(1:length(myfpkm[,1]),3, replace=FALSE),13:15])),.95)
cut_off_6 <- quantile(replicate(1000,p_corr(myfpkm[sample(1:length(myfpkm[,1]),3, replace=FALSE),16:18])),.95)
cut_off_7 <- quantile(replicate(1000,p_corr(myfpkm[sample(1:length(myfpkm[,1]),3, replace=FALSE),19:21])),.95)
cut_off_8 <- quantile(replicate(1000,p_corr(myfpkm[sample(1:length(myfpkm[,1]),3, replace=FALSE),22:24])),.95)
cut_off <- quantile(replicate(1000,p_corr(myfpkm[sample(1:length(myfpkm[,1]),3, replace=FALSE),])),.95)

get_seeds <- function(X,cut_off=0.573) {
	sig <- X[,ncol(X)]
	X <- X[,-ncol(X)]
	if(sum(sig)<=1) {return(FALSE)}	
	if(p_corr(X)>cut_off) {return(TRUE)}
	return(FALSE)
}

myfpkm <- myfpkm[order(rownames(myfpkm)),]
all_results <- all_results[order(rownames(all_results)),]
myfpkm$test <- all_results$test
myfpkm$test[is.na(myfpkm$test)] <- FALSE

test<-c()
for (i in (1:(nrow(myfpkm)-2))) {
	test <- c(test,get_seeds(myfpkm[i:(i+2),],cut_off))
}
myfpkm$seed <- c(test,FALSE,FALSE)




### this is all cool, but it didn't do what I wanted...
p_corr <- function(x) {
	X <- cbind(
		t(x[1:3]),
		t(x[4:6]),
		t(x[7:9]),
		t(x[10:12]),
		t(x[13:15]),
		t(x[16:18]),
		t(x[19:21]),
		t(x[22:24])
	)
	t1 <- suppressWarnings(cor(X,use="all.obs",method="pearson"))
	t1[!lower.tri(t1)] <- NA
	t1 <- sqrt(t1^2)
	sum(t1,na.rm=T)/(0.5*(length(x)^2-length(x)))
}

#set.seed=0.24

quantile(
	apply(
		myfpkm[sample(1:length(myfpkm[,1]), 1000, replace=FALSE),]
		,1
		,function(x) {
			p_corr(as.data.frame(t(x)))
		}
	),.95
)
