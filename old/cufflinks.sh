#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})
BAM=$1
PROC=$2

qsub $SCRIPT_DIR/submit_cufflinks.sh $BAM $PROC
