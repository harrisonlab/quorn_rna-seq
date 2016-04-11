# quorn_rna-seq

Set quorn project file to:
```shell
QUORN=~/projects/quorn
```

## Trim data
Trimming was performed with Trimmomatic 

```shell
counter=0
for f in $QUORN/raw_rna/*
do
	counter=$((counter+1))
	if (( $counter % 2 == 0 )) 
	then
		R2=$f
	    $QUORN/scripts/trim.sh $R1 $R2 $QUORN/scripts 25 100
	fi
	R1=$f
done
```

## Filter data

## Align to ref with Tophat

## Cufflinks

## Cuffmerge

## Cuffquant

## Cuffdiff
