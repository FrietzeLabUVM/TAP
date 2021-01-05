#!/bin/bash
#$ -cwd
#$ -N exactSNP
#py=/slipstream/home/dbgap/anaconda3/bin/python3
bam=$1
fasta=$2
if [ -z $bam ]; then echo must have bam file as arg1. quit.; exit 1; fi
if [ ! -s $bam ]; then echo could not find bam file $bam. quit.; exit 1; fi

if [ ! -s ${bam/.bam/""}.all.exactSNP.vcf ]; then
  exactSNP -b -s 20 -C ~/tmp -i $bam -g $fasta -o ${bam/.bam/""}.all.exactSNP.vcf
else
  echo found ${bam/.bam/""}.exactSNP.vcf. no rerun.
fi
