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
```shell
counter=0
for f in $QUORN/filtered/*.*
do
    counter=$((counter+1))
    if (( $counter % 2 == 0 )) 
    then
        R2=$f
        $QUORN/scripts/tophat.sh $R1 $R2 $QUORn/ref/venenatum 200 400
    fi
    R1=$f
done
```
## Cufflinks
second argument to cufflinks.sh is no. processors

This will need editing or cufflinks will output the same file name multiple times.... 
```shell
for f in $QUORN/tophat/*.gz
do
	$QUORN/scripts/cufflinks.sh $f 8
done
```
## Cuffmerge

## Cuffquant

## Cuffdiff
