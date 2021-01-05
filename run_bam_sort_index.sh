#!/bin/bash
#$ -N sort_index

#SBATCH --nodes=1                               # Request one core
#SBATCH --ntasks-per-node=1                               # Request one node (if you request more than one core with -n, also using
#SBATCH --cpus-per-task=1                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=8000                        # Memory total in MB (for all cores)
#SBATCH -o index_bam_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e index_bam_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-user=jrboyd@med.uvm.edu   # Email to which notifications will be sent


b=$1
if [ -z $b ]; then echo bam file expected as arg1, quit.; exit 1; fi

ob=${b/.Aligned.out.bam/""}.Aligned.sortedByCoord.out.bam
if [ -f $ob ]; then echo sorted bam file $ob found! rm to rerun; exit 0; fi

if [ ! -f $b ]; then echo bam file $b not found! quit.; exit 1; fi

if [ -f ${ob}.bai ]; then rm ${ob}.bai; fi
mkdir -p ~/tmp
samtools sort -T ~/tmp -o $ob $b
samtools index $ob

if [ -f ${ob}.bai ]; then rm $b; fi
echo $ob

