#!/bin/bash
input=~/lab_shared/data/Frietze_SK_iLabs_15333_Cahan05_12092020
ref=~/lab_shared/indexes/DM6

for f in $input/*/*_R1_001.fastq.gz; do
  echo $f;
  bash rnaseq_pipeline.sh -f1 $f -ref $ref
done
