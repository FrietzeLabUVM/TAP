#!/bin/bash
#$ -cwd
#$ -N exactSNP
echo $0 $@

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -b|--bam) BAM="$2"; shift ;;
        -fa|--fasta) FASTA="$2"; shift ;;
        -docker|--docker) docker="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ -z $BAM ]; then echo must have bam file as arg1. quit.; exit 1; fi
if [ ! -s $BAM ]; then echo could not find bam file $BAM. quit.; exit 1; fi
OUT_SNP_FILE=${BAM/.bam/""}.all.exactSNP.vcf
OUT_SNP_FILE_local=$OUT_SNP_FILE

run="no"
if [ ! -s $OUT_SNP_FILE ]; then
  run="yes"
fi

TMP_DIR=~/tmp
mkdir -p $TMP_DIR

# docker for samtools and UCSC tools v1.0
echo docker is $docker
if [ -n "$docker" ]; then
  dBAM=/input_bam/$(basename $BAM)
  dFASTA=/input_fasta/$(basename $FASTA)
  dOUT_SNP_FILE=/output/$(basename $OUT_SNP_FILE)
  dTMP_DIR=/user_tmp

  base_cmd="docker run \
    -u $(id -u):$(id -g) \
    -v $(dirname $BAM):$(dirname $dBAM) \
    -v $(dirname $FASTA):$(dirname $dFASTA) \
    -v $(dirname $OUT_SNP_FILE):$(dirname $dOUT_SNP_FILE) \
    -v $TMP_DIR:$dTMP_DIR \
    --entrypoint"
    
  cmd_exactSNP="$base_cmd exactSNP $docker"
  
  BAM=$dBAM
  FASTA=$dFASTA
  OUT_SNP_FILE=$dOUT_SNP_FILE
  TMP_DIR=$dTMP_DIR
else
  cmd_exactSNP=exactSNP
fi


if [ $run = "yes" ]; then
  run_cmd="$cmd_exactSNP -b -s 20 -C $TMP_DIR -i $BAM -g $FASTA -o $OUT_SNP_FILE"
  echo exactSNP command:
  echo $run_cmd
  $run_cmd
else
  echo found ${OUT_SNP_FILE_local}. no rerun.
fi
