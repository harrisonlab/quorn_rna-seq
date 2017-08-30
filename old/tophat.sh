#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

R1=$1
R2=$2
REF=$3
I=$4
X=$5

qsub $SCRIPT_DIR/submit_tophat.sh $R1 $R2 $REF $I $X 
