#!/bin/bash
#SBATCH --nodes=1                               # Request one core
#SBATCH --ntasks-per-node=1                               # Request one node (if you request more than one core with -n, also using
#SBATCH --cpus-per-task=8                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p bluemoon                           # Partition to run in
#SBATCH --mem=32000                        # Memory total in MB (for all cores)
#SBATCH -o star_pe_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e star_pe_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-user=jrboyd@med.uvm.edu   # Email to which notifications will be sent



#alignment script for dbgap data
#assumes data is paired end fastq.gz and files end in _R1_001.fastq.gz

echo $0 $@

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f1|--fastq1) F1="$2"; shift ;;
        -f2|--fastq2) F2="$2"; shift ;;
        -o|--out) B="$2"; shift ;;
        -wd|--workdir) wd="$2"; shift ;;
        -idx|--starindex) star_idx="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

#umask 077 #limits rwx permission to user only
#wd=/slipstream/home/dbgap/data/fetched_fastq_RNA-Seq/alignment2
if [ -z $wd ]; then wd=$(dirname $F1)/alignment; fi
mkdir -p $wd

if [ -z $F1 ]; then echo "expect read1 fastq as arg1, stop"; exit 1; fi
if [ ! -f $F1 ]; then echo "read1 fastq file, $F1, not found! stop"; exit 1; fi
#F2=$2if [ -z $F2 ]; then F2=${F1/_R1/_R2}; if [ ! -f $F2 ]; then echo "attempt at guessing read2 fastq as $F2 failed, file not found. try supplying as arg2. stop"; exit 1; fi; fi
if [ -z $F2 ]; then F2=${F1/_R1/_R2}; if [ ! -f $F2 ]; then echo "attempt at guessing read2 fastq as $F2 failed, file not found. try supplying as -f2. stop"; exit 1; fi; fi
if [ ! -f $F2 ]; then echo "read2 fastq file, $F2, not found! stop"; exit 1; fi

#B=$3
if [ -z $B ]; then B=$(basename ${F1/_R1_001.fastq.gz/""}); fi
B=$wd/$(basename $B)


echo read1 fastq is $F1
echo read2 fastq is $F2
echo output prefix is $B

if [ ! -d $(dirname $B) ]; then
  mkdir -p $(dirname $B)
fi

if [ -f ${B}.Aligned.sortedByCoord.out.bam ]; then # && [ -f ${B}.Aligned.toTranscriptome.out.bam ]; then
 echo output bam ${B}.Aligned.sortedByCoord.out.bam exists! will not rerun.
 exit 0
fi

if [ -f ${B}.Aligned.out.bam ]; then # && [ -f ${B}.Aligned.toTranscriptome.out.bam ]; then
 echo output bam ${B}.Aligned.out.bam exists! will not rerun.
 exit 0
fi

STAR \
--runThreadN 8 \
--readFilesIn $F1 $F2 \
--readFilesCommand gunzip -c \
--genomeDir $star_idx \
--runDirPerm User_RWX \
--outFilterType BySJout `#encode params` \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFileNamePrefix $B"." `#output params` \
--outSAMtype BAM Unsorted \
--quantMode TranscriptomeSAM GeneCounts \
--twopassMode Basic

#samtools index ${B}.Aligned.sortedByCoord.out.bam
#samtools index ${B}.Aligned.toTranscriptome.out.bam
