#!/bin/bash
#SLURM pipeline for RNAseq

in_param="$@"

set -- -f1 test.1.fastq.gz -p my_prefix -ref my_ref
echo "$#"
echo "$@"

#set -a -b -e -f arg1
#echo "$#"
#echo "$@"

#exit 0

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f1|--fastq1) F1="$2"; shift ;;
        -f2|--fastq2) F2="$2"; shift ;;
        -p|--prefix) root="$2"; shift ;;
        -o|--outdir) align_path="$2"; shift ;;
        -ref|--reference) ref="$2"; shift ;;
        -idx|--starindex) star_index="$2"; shift ;;
        -s|--supparef) suppa_ref="$2"; shift ;;
        -g|--gtf) gtf="$2"; shift ;;
        -fa|--fasta) fasta="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

echo from config
echo $F1
echo $root
echo $ref

set -- ${in_param}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f1|--fastq1) F1="$2"; shift ;;
        -f2|--fastq2) F2="$2"; shift ;;
        -p|--prefix) root="$2"; shift ;;
        -o|--outdir) align_path="$2"; shift ;;
        -ref|--reference) ref="$2"; shift ;;
        -idx|--starindex) star_index="$2"; shift ;;
        -s|--supparef) suppa_ref="$2"; shift ;;
        -g|--gtf) gtf="$2"; shift ;;
        -fa|--fasta) fasta="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

echo from command
echo $F1
echo $root
echo $ref

