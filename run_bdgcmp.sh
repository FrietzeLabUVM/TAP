#!/usr/bin/bash
#$ -o bdgcmp."$JOB_ID".out
#$ -e bdgcmp."$JOB_ID".error


#required inputs:
#TREATMENT - control lambda bdg file from macs2
#CONTROL - treatment pileup bdg file from macs2
#METHOD - should probably be logFE, see macs2 bdgcmp help for more options
WD=$WD
t=$TREATMENT
c=$CONTROL
met=$METHOD
if [ -z $t ]
	then t=$1
	c=$2
	met=$3
	WD=$(pwd)
fi
echo treat file - $t
echo control file - $c
echo method - $met
OUT=$WD/$(basename $t)
OUT="${OUT/_treat_pileup.bdg/}"_"$met".bdg
if [ -f $OUT ] || [ -f ${OUT/.bdg/.bw} ]; then
	echo skipping bdgcmp for "$t", file "$OUT" exists
else
	macs2 bdgcmp -t $t -c $c -m $met -o $OUT
fi
