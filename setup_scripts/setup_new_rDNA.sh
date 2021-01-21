#!/bin/bash
#SBATCH --nodes=1                               # Request one core
#SBATCH --ntasks-per-node=1                               # Request one node (if you request more than one core with -n, also using
#SBATCH --cpus-per-task=1                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=8000                        # Memory total in MB (for all cores)
#SBATCH -o ref_build_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e ref_build_%j.err                 # File to which STDERR will be written, including job ID


while [[ "$#" -gt 0 ]]; do
  case $1 in
    -o|--output) out="$2"; shift ;;
    -f|--fasta) fa="$2"; shift ;;
    --genomeSAindexNbases) saIndex="$2"; shift ;;
    *) echo "Unknown parameter passed: $1"; exit 1 ;;
  esac
  shift
done

if [ -z $out ] || [ -z $fa ]; then echo output \(-o\|--output\) and fasta \(-f\|--fasta\) are all required. quit.; exit 1; fi

echo out will be $out
echo fasta will be $fa

#if [ -d $out ]; then echo "output location $out exists! delete to or change location. quit."; exit 1; fi
if [ ! -f $fa ]; then echo fasta file $fa not found! quit.; exit 1; fi
if [ -z $saIndex ]; then saIndex=5; fi

echo saIndex will be $saIndex \(set below 14 for smaller genomes\)


fa=$(readlink -f $fa)

mkdir -p $out
out=$(readlink -f $out)
cd $out

bfa=$(basename $fa)

echo staging gtf and fasta...

#fasta
mkdir FASTA
cd FASTA
cp $fa $bfa
ln -sf $bfa genome.fa
cd ..

echo making STAR index...

#star index
STAR --runThreadN 1 \
--runMode genomeGenerate \
--genomeDir $out/STAR_INDEX \
--genomeFastaFiles FASTA/genome.fa \
--genomeSAindexNbases $saIndex
