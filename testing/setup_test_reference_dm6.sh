#!/bin/bash


while [[ "$#" -gt 0 ]]; do
  case $1 in
    -docker|--docker) docker="$2"; shift ;;
    -singularity|--singularity) singularity="$2"; shift ;;
    -p|--runThreadN) CPUS="$2"; shift ;;
    *) echo "Unknown parameter passed: $1"; exit 1 ;;
  esac
  shift
done

SCRIPT_PATH=$(dirname "$(readlink -f "$0")")
if [ ! -f dm6.fa ]; then
 wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
 gunzip dm6.fa.gz
fi
if [ ! -f dm6.ensGene.gtf.gz ]; then
 wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/genes/dm6.ensGene.gtf.gz
fi

extra_args=""
if [ -n "$docker" ]; then extra_args=$extra_args"--docker $docker "; fi
if [ -n "$singularity" ]; then extra_args=$extra_args"--singularity $singularity "; fi
if [ -n "$CPUS" ]; then extra_args=$extra_args"--runThreadN $CPUS "; fi
bash "$SCRIPT_PATH"/../setup_scripts/setup_new_reference.sh -o "$SCRIPT_PATH"/references/dm6 --gtf_ensemble dm6.ensGene.gtf.gz -f dm6.fa --genomeSAindexNbases 12 $extra_args
