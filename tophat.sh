#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

R1=$1
R2=$2
REF=$3
OUT=$4
I=$5
X=$6

qsub $SCRIPT_DIR/submit_tophat.sh $R1 $R2 $REF $OUT $I $X 
