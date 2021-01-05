#!/bin/bash
#$ -cwd
#$ -N featureCounts
bam=$1
gtf=$2
if [ -z $bam ]; then echo must have bam file as arg1. quit.; exit 1; fi
if [ ! -f $bam ]; then echo could not find bam file $bam. quit.; exit 1; fi

if [ ! -f ${bam/.bam/""}.featureCounts.txt ]; then
  featureCounts -p -B -s 2 -t exon -g gene_id -a ${gtf} -o ${bam/.bam/""}.featureCounts.txt $bam
else
  echo found ${bam/.bam/""}.featureCounts.txt. no rerun.
fi
