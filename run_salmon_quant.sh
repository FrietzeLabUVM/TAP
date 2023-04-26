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

echo $0 $@

BAM=$1
GTF=$2
container=$4

if [ -z $BAM ]; then echo need bam as arg1, quit.; exit 1; fi
if [ ! -f $BAM ]; then echo can not find bam $BAM, quit.; exit 1; fi

OUT=${BAM/.Aligned.toTranscriptome.out.bam/""}.salmon_quant

if [ $OUT == ${BAM}.salmon_quant ]; then
  echo input bam does not seem valid, should end in .Aligned.toTranscriptome.out.bam.
  exit 1
fi

mkdir -p $OUT

# container for suppa2 v1.1
echo container is $container
container_type=""
if [ -n "$container" ]; then container_type="docker"; fi
if [[ "$container" == *.sif ]]; then container_type="singularity"; fi

echo container_type is $container_type

if [ -n "$container_type" ]; then
  dBAM=/input_bam/$(basename $BAM)
  dGTF=/input_gtf/$(basename $GTF)
  dOUT=/output_salmon

  if [ $container_type = "docker" ]; then
    cmd_salmon="docker run \
      -u $(id -u):$(id -g) \
      -v $(dirname $BAM):$(dirname $dBAM) \
      -v $(dirname $GTF):$(dirname $dGTF) \
      -v $OUT:$dOUT \
      --entrypoint \
      salmon $container"  
  elif [ $container_type = "singularity" ]; then
    cmd_salmon="singularity exec \
      --bind $(dirname $BAM):$(dirname $dBAM),$(dirname $GTF):$(dirname $dGTF),$OUT:$dOUT \
      $singularity \
      $container salmon"
  else
      echo "Unrecognized container_type $container_type";
      exit 1;
  fi
  BAM=$dBAM
  GTF=$dGTF
  OUT=$dOUT
else
  cmd_salmon=salmon
fi

#if [ -d $OUT ]; then echo output $OUT already exists, will not rerun.; exit 0; fi
$cmd_salmon quant -p 2 -t $GTF -l A -a $BAM -o $OUT --gencode
