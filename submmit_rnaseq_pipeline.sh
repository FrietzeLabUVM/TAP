#!/bin/bash
wd=~/lab_shared/data/Frietze_SK_iLabs_15333_Cahan05_12092020
idx=~/dbgap_scripts/vacc_scripts/DM6_STAR_INDEX
gtf=~/lab_shared/indexes/DM6/GTF/dm6.ensGene.gtf
fasta=~/lab_shared/indexes/DM6/FASTA/dm6.fa

for f in $wd/*/*_R1_001.fastq.gz; do
  echo $f;
  bash rnaseq_pipeline.sh -f1 $f -idx $idx -g $gtf -fa $fasta
done
