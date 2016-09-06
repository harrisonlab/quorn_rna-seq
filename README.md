# quorn_rna-seq

Set quorn project file to:
```shell
QUORN=~/projects/quorn
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

## DESeq2 analysis
Using braker gene models (cufflinks is still running after a couple of weeks)

The method in dge_deseq.R was followed to produce list of diffrentially expressed braker gene models.

## Clusters
co_clusters.R will find groups of n consecutive genes with expression correlation higher than the .95 quantile of n random genes (n set to 3 by default). Still in development, but will work after a fashion.

## promoter motif finder
get_primers.pl will return n nucleotides upstream of a feature in a gff file.

```shell
# get 1000  nucleotides upstream of each gene in final_genes_Braker.gff
./get_primer.pl contigs_unmasked.fa final_genes_Braker.gff 1000 gene >brake_up1000.fa
#motif is GTGA...GTGA seperated by at most 8 nucleotides.
grep -E -o -B1 GTGA.\{0,8\}GTGA brake_up1000.fa
```

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
