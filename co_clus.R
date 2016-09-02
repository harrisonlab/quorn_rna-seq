res2 <- res.2
res3 <- res.3
res4 <- res.4
res5 <- res.5
res6 <- res.6
res7 <- res.7
res8 <- res.8

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
