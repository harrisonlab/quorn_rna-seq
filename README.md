# quorn_rna-seq

Set quorn project file to:
```shell
QUORN=~/projects/quorn
```

## QC
Qualtiy checking with fastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

?? This needs updating - copied from another project  
```shell
for FILE in $ARDERI/data/$RUN/fastq/*; do 
	$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c qcheck $FILE $ARDERI/data/$RUN/quality
done
```

## Trim data
Trimming was performed with Trimmomatic (trim.sh, submit_trim.sh and truseq.fa should all be in same directory)
Around 25% - 30% of reverse reads were discarded due to adapter contamination. Trimmomatic was set to capture (rather than dump) unpaired forward reads. SE workflow to follow...


```shell
counter=0
for f in $QUORN/raw_rna/*.gz
do
	counter=$((counter+1))
	if (( $counter % 2 == 0 )) 
	then
		R2=$f
		$QUORN/scripts/trim.sh $R1 $R2 $QUORN/scripts 20 50
	fi
	R1=$f
done
```

## Filter data
```shell
counter=0
for f in $QUORN/trimmed/*.gz
do
	counter=$((counter+1))
	if (( $counter % 2 == 0 )) 
	then
		R2=$f
		S=$(echo $f|awk -F"/" '{print $NF}'|awk -F"_" '{print $1,$2,$3}' OFS="_")
		$QUORN/scripts/bowtie.sh $R1 $R2 $QUORN/filtered/phix/phix $QUORN/filtered/$S 200 400
	fi
	R1=$f
done
```

## Align to ref with Tophat 
(maybe update to Hisat2/Tophat3 (when available))
Unpaired reads can be added to the tophat workflow
```shell

for f in ../filtered/*[1].fq; 
do 
	S=$(echo $f|awk -F"/" '{print $NF}'|awk -F"." '{print $1}')
	$QUORN/scripts/tophat.sh $QUORN/filtered/${S}.1.fq $QUORN/filtered/${S}.2.fq $QUORN/filtered/${S}_1_SE.fq $QUORN/ref/venenatum $QUORN/filtered/$S 200 400  
done

```

## NEW - Align to ref with STAR - 
I prefer STAR now - it's performance is not so depedent on choice of input parameters.
```shell
# Create star index
STAR --runMode genomeGenerate --genomeDir your_out_dir --genomeFastaFiles redgauntlet.fa --sjdbGTFfile redgauntlet.gff
# align 
STAR --genomeDir your_out_dir --outFileNamePrefix something --readFilesIn fastq_F fastq_R --outSAMtype SAM --runThreadN 16

```
## Count features
Using featureCounts. 

This can be done in R, but is slow and either imports all files or fails. Outside R each sample is counted seperately.

```shell
featureCounts -o output_file -a gff_file sam_files

```

Can be useful to produce a SAF file from an input gff (as they are not always consistent)
The below will extract exon annotations and output the ninth column stipped of ID= and anything after the first ".", then the first column (chromosome) and etc.
```
grep exon final_genes_appended.gff3|awk -F"\t" '{gsub(/ID=/,"",$NF);gsub(/\..*/,"",$NF);print $NF,$1,$4,$5,$7}' OFS="\t" > $QUORN/counts/exons.SAF
```

The RNA-seq pipeline can be used to run featureCounts:
PIPELINE.sh -c counts annotaions output_dir output_file sam/bam(s) [options]

The below runs with 12 threads (-T 12), counts all multimapping reads (-M) and uses a SAF file for input (-F SAF)
```
for D in $QUORN/aligned/treatment/WTCHG*; do
OUTFILE=$(echo $D|awk -F"/" '{print $(NF)}').counts
$QUORN/RNA-seq_pipeline/scripts/PIPELINE.sh -c counts \
$QUORN/counts/exons.SAF \
$QUORN/counts \
$OUTFILE \
$D/star_aligmentAligned.sortedByCoord.out.bam -T 12 -M -F SAF
done
```

### DESeq2 analysis
Follow script new_DEG.R

## Functional analysis
GO funtional analysis with BinGO
Requires GO.obo for gene ontology (downloaded 22/06/17)
IPR_2_GO (downloaded version 17/06/17)

```shell
# produce simplfied ipr2go with GO_ID\tIPR_ID
#tail -n +7 iprtogo.txt|awk -F" " '{gsub(/InterPro:/, "",$1);print $1,$NF}' OFS="\t" > ipr2go.txt # remove all annotation
grep -oP "IPR\d*|GO:.*" iprtogo.txt|awk '{if($1~/^I/){i=$1;}else{print i,$0}}' OFS=";" > ipr2go.txt # retain go annotation
# produce simplified annotawion with GENE_ID\tIPR_ID
grep -oP "g\d*\.t1|IPR\d*" annotation.txt|awk '{if($1~/^g/){i=$1;}else{print i,$1}}' OFS="\t"|sort|uniq >gene_ipr.txt
```

```R
# merge ipr2go with gene_ipr to output GENE_ID\tGO_ID
library(data.table)
library(dplyr)
gipr <- fread("gene_ipr.txt",header=F)
iprgo <- fread("ipr2go.txt",header=F,sep=";")
output <- inner_join(gipr,iprgo,by=c("V2"="V1"))
colnames(output) <- c("GENE_ID","IPR_ID","IPR_DESC","GO_DESC","GO_ID")
output$GENE_ID <- sub("\\.t1.*","",output$GENE_ID)
output$GO_DESC <- sub("GO:","",output$GO_DESC)
# output<-cbind(a="xx",output[,1:2],b="",output[,4],c="","ISS","UNKNOWN","C",output[,3],"gene","taxon:5555","210617","GD")
# output <- output[complete.cases(output),]
# output$unique<-paste(output[,2],output[,5],sep="_")
# output <- output[!duplicated(output$unique),]
write.table(output,"gene_IPR_GO.txt",row.names=F,quote=F,sep="\t",na="")
# write.table(output[,1:14],"genes_go.txt",sep="\t",row.names=F,col.names=F,quote=F)
```

BinGO requires a gene2go file in a specific format, with

## Clusters
co_clusters.R will find groups of n consecutive genes with expression correlation higher than the .95 quantile of n random genes (n set to 3 by default). Still in development, but will work after a fashion.
For n = 3 and using correlation between all datasets (probably better to limit it to each experiment) 1406 clusters, just over 11% of total possible clusters have correlation higher than the .95% quantile

cluster length|total
---|---
13|1
12|1
11|1
10|4
9|6
8|14
7|26
6|57
5|104
4|171


## promoter motif finder
get_primers.pl will return n nucleotides upstream of a feature in a gff file.

```shell
# get 1000  nucleotides upstream of each gene in final_genes_Braker.gff
./get_primer.pl contigs_unmasked.fa final_genes_Braker.gff 1000 gene >brake_up1000.fa
#motif is GTGA...GTGA seperated by at most 8 nucleotides.
grep -E [TGA].AGGCC brake_up1000.fa wc -l # 4797 
grep -E '(GTGA|TCAC)'.\{0,8\}'(GTGA|TCAC)' brake_up1000.fa|wc -l # 5579
grep -E [TGA].AGGCC brake_up1000.fa|grep -E  '(GTGA|TCAC)'.\{0,8\}'(GTGA|TCAC)' |wc -l # 2092
grep -E -B1 [TGA].AGGCC brake_up1000.fa|grep -E -B1 '(GTGA|TCAC)'.\{0,8\}'(GTGA|TCAC)'|grep ">" >set_both.txt
sed -i -e 's/>//' set_both.txt
```
The number of differentially expressed genes containing both promters per media is given below. 

Media|DGE|set_both
---|---|---
2793|1616|39
F55|2280|153
10170|4016|520
MWT|6999|1058
MOL|5012|635
MKO|6988|1056
TJ|4338|476

Maybe of interest, given the number of genes in this set all the media show lower numbers of DE genes compared to what might be expected.

##Not implemented

### Cufflinks
second argument to cufflinks.sh is no. processors

This will need editing or cufflinks will output the same file name multiple times.... 
```shell
for f in $QUORN/tophat/*.gz
do
	$QUORN/scripts/cufflinks.sh $f 8
done
```
