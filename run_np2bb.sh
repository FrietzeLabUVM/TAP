#!/usr/bin/bash
#$ -o np2bb."$JOB_ID".out
#$ -e np2bb."$JOB_ID".error


#required inputs:
if [ -z $inputNarrowPeak ]
then
	inputNarrowPeak=$1
	inputChromSizes=$2
	outputBigBed=$3
fi
echo input - $inputNarrowPeak
echo chrSizes - $inputChromSizes
echo output - $outputBigBed
if [ -f $outputBigBed ]; then
  echo file $outputBigBed exists so np2bb not necessary for $inputBegGraph
  echo nothing done.
else
  f_bed=${inputNarrowPeak}.2bb.bed
  sort -k1,1 -k2,2n $inputNarrowPeak | awk 'BEGIN {FS="\t"; OFS="\t"} {if ($5 > 1000) $5 = 1000; print $1,$2,$3,$4,$5}' > $f_bed
  bedToBigBed $f_bed $inputChromSizes $outputBigBed
  rm $f_bed
fi

