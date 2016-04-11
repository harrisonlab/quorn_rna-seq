# quorn_rna-seq

Set quorn project file to:
```shell
QUORN=~/projects/quorn
```

## Trim data
Trimming was performed with Trimmomatic (trim.sh, submit_trim.sh and illumina_full_adapters.fa should all be in same directory)

```shell
counter=0
for f in $QUORN/raw_rna/*.gz
do
	counter=$((counter+1))
	if (( $counter % 2 == 0 )) 
	then
		R2=$f
		$QUORN/scripts/trim.sh $R1 $R2 $QUORN/scripts 25 150
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
		$QUORN/scripts/bowtie.sh $R1 $R2 $QUORN/filtered/phix/phix $QUORN/filtered 200 400
	fi
	R1=$f
done
```
## Align to ref with Tophat 
(maybe update to Hisat2/Tophat3 (when available))
```shell
counter=0
for f in $QUORN/filtered/*.gz
do
	counter=$((counter+1))
	if (( $counter % 2 == 0 )) 
	then
		R2=$f
		$QUORN/scripts/tophat.sh $R1 $R2 REF 200 400
	fi
	R1=$f
done
```
## Cufflinks

## Cuffmerge

## Cuffquant

## Cuffdiff
