#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

R1=$1
R2=$2
TRIMLOC=$3
QUALITY=$4
MINLEN=$5

		qsub $SCRIPT_DIR/submit_trim.sh $R1 $R2 $TRIMLOC $QUALITY $MINLEN
