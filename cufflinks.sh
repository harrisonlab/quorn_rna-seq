#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})
BAM=$1
OUT=$2
PROC=$3



qsub $SCRIPT_DIR/submit_cufflinks.sh $BAM $OUT $PROC
