#!/bin/bash

#Assemble short reads with Tophat
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

R1=$1
R2=$2
REF=$3
I=$4
X=$5

echo "Running Tophat with the following in= REF IS '$REF' READ 1 '$R1' READ 2 ' $R2 ' DEST is '$DEST' OUTNAME IS '$OUTNAME' MIN is '$I' MAX is '$X'"

tophat $REF $R1 $R2 -p 8 --fusion-search --b2-I $I --b2-X $X


