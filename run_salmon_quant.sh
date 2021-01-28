#!/bin/bash
#$ -cwd
#$ -V
#$ -N salmon
#$ -pe threads 4
#umask 077
#SBATCH --cpus-per-task=4                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p bluemoon                           # Partition to run in
#SBATCH --mem=8000   

eval "$(/gpfs2/pi-sfrietze/anaconda3/bin/conda shell.bash hook)"
conda activate /gpfs2/pi-sfrietze/envs/salmon_env 

echo $0 $@

bam=$1
gtf=$2
if [ -z $bam ]; then echo need bam as arg1, quit.; exit 1; fi
if [ ! -f $bam ]; then echo can not find bam $bam, quit.; exit 1; fi

out=${bam/.Aligned.toTranscriptome.out.bam/""}.salmon_quant

if [ $out == ${bam}.salmon_quant ]; then
  echo input bam does not seem valid, should end in .Aligned.toTranscriptome.out.bam.
  exit 1
fi

#if [ -d $out ]; then echo output $out already exists, will not rerun.; exit 0; fi
salmon quant -p 2 -t $gtf -l A -a $bam -o $out --gencode
