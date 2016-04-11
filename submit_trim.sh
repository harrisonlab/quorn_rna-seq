#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

FORWARD=$1
REVERSE=$2
TRIMLOC=$3
QUALITY=$4
MINLEN=$5

java -jar $TRIMLOC/trimmomatic-0.33.jar PE -phred33 $1 $2 $1.trimmed.fq /dev/null $2.trimmed.fq /dev/null ILLUMINACLIP:$TRIMLOC/illumina_full_adapters.fa:2:30:10 SLIDINGWINDOW:8:$QUALITY MINLEN:$MINLEN

