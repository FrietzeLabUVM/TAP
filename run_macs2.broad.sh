#!/bin/bash
#$ -o macs2."$JOB_ID".out
#$ -e macs2."$JOB_ID".error


echo TREAT_BAM is $TREAT_BAM
echo INPUT_BAM is $INPUT_BAM
echo OUTDIR is $OUTDIR
echo PREFIX is $PREFIX
echo PVAL is $PVAL
echo QVAL is $QVAL
echo BROADCUTOFF is $BROADCUTOFF
echo GEN is $GEN
echo macs2 on $(basename $TREAT_BAM):$(basename $INPUT_BAM)

if [ -z "$TREAT_BAM" ] #|| [ ! -e $TREAT_BAM ]
then
	echo TREAT_BAM $TREAT_BAM arg required! stop
	exit 1
fi
for tb in $TREAT_BAM; do
  if [ ! -e $tb ]; then
    echo TREAT_BAM $tb not found! stop
    exit 1
  fi
done

if [ -z "$INPUT_BAM" ] #|| [ ! -e $INPUT_BAM ]
then
	echo INPUT_BAM $INPUT_BAM arg required! stop
	exit 1
fi
for ib in $INPUT_BAM; do
  if [ ! -e $ib ]; then
    echo INPUT_BAM $ib not found! stop
    exit 1
  fi
done

if [ -z $OUTDIR ] || [ ! -d $OUTDIR ]
then
	echo OUTDIR $OUTDIR not found! stop
	exit 1
fi
if [ -z $PREFIX ]
then
	echo PREFIX $PREFIX missing! stop
	exit 1
fi


if [ -z $PVAL ]
then
  if [ -z $QVAL ]; then
    echo PVAL and QVAL missing, need 1! stop
    exit 1
  else
    stat="q"
    stat_val=$QVAL
  fi
else
  stat="p"
  stat_val=$PVAL
fi

if [ $GEN == "mm10" ] || [ $GEN == "mm9" ]; then
	g=mm
elif [ $GEN == "hg38" ] || [ $GEN == "hg19" ] || [ $GEN == "U13369" ]; then
	g=hs
fi

#model is disabled for testing on small files!
if [ -f $PREFIX"_peaks.broadPeak" ]; then
	echo $PREFIX"_peaks.broadPeak" exists. macs2 has already been run for $TREAT_BAM.
	echo delete $PREFIX"_peaks.broadPeak" if you want to rerun macs2.
elif [[ $GEN == *"simGenome"* ]]; then
        cmd="macs2 callpeak -t $TREAT_BAM -c $INPUT_BAM -g $g --outdir $OUTDIR -n $PREFIX -$stat $stat_val --broad --broad-cutoff $BROADCUTOFF --nomodel --extsize 200"
else
        cmd="macs2 callpeak -t $TREAT_BAM -c $INPUT_BAM -g $g --outdir $OUTDIR -n $PREFIX -$stat $stat_val --broad --broad-cutoff $BROADCUTOFF"
fi
echo cmd is:
echo $cmd
$cmd

