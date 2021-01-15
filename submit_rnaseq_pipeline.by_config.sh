#!/bin/bash
#input=~/lab_shared/data/Frietze_SK_iLabs_15333_Cahan05_12092020
cfg=DM6_AA_config.csv
ref=~/lab_shared/indexes/DM6
out=~/lab_shared/AA_DM6_by_config

todo=$(cat $cfg)
for l in $todo; do 
  f1=$(echo $l | awk -v FS="," '{print $1}'); 
  root=$(echo $l | awk -v FS="," '{print $2}'); 
  echo fastq $(basename $f1) to root $root in $out
  bash rnaseq_pipeline.sh -f1 $f1 -ref $ref -o $out -p $root
done

