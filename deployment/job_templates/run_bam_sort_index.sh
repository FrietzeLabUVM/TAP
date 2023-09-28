#!/bin/bash
#SBATCH --nodes=1                         # Request one core
#SBATCH --ntasks-per-node=1               # Request one node (if you request more than one core with -n, also using
#SBATCH --cpus-per-task=1                 # -N 1 means all cores will be on the same node)
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p short                          # Partition to run in
#SBATCH --mem=8000                        # Memory total in MB (for all cores)
#SBATCH -o index_bam_%j.out               # File to which STDOUT will be written, including job ID
#SBATCH -e index_bam_%j.err               # File to which STDERR will be written, including job ID

echo $0 $@

BAM=$1
container=$3

if [ -z $BAM ]; then echo bam file expected as arg1, quit.; exit 1; fi
echo BAM is $BAM
OUT_BAM=${BAM/.Aligned.out.bam/""}.Aligned.sortedByCoord.out.bam
echo OUT_BAM is $OUT_BAM
if [ -f $OUT_BAM ]; then echo sorted bam file $OUT_BAM found! rm to rerun; exit 0; fi

if [ ! -f $BAM ]; then echo bam file $BAM not found! quit.; exit 1; fi

if [ -f ${OUT_BAM}.bai ]; then rm ${OUT_BAM}.bai; fi
TMP_DIR=~/tmp
mkdir -p $TMP_DIR

lOUT_BAM=$OUT_BAM
lBAM=$BAM

# container for suppa2 v1.1
echo container is $container
container_type=""
if [ -n "$container" ]; then container_type="docker"; fi
if [[ "$container" == *.sif ]]; then container_type="singularity"; fi

echo container_type is $container_type

if [ -n "$container_type" ]; then
  #derive container paths
  dBAM=/input_bam/$(basename $BAM)
  dOUT_BAM=/output_bam/$(basename $OUT_BAM)
  dTMP_DIR=/user_tmp

  echo $container_type BAM is $dBAM
  echo $container_type OUT_BAM is $dOUT_BAM
  echo $container_type TMP_DIR is $dTMP_DIR

  if [ $container_type = "docker" ]; then
    cmd_samtools="docker run \
      -u $(id -u):$(id -g) \
      -v $(dirname $BAM):$(dirname $dBAM) \
      -v $(dirname $OUT_BAM):$(dirname $dOUT_BAM) \
      -v $TMP_DIR:$dTMP_DIR \
      --entrypoint \
      samtools $container"
  elif [ $container_type = "singularity" ]; then
    cmd_samtools="singularity exec \
      --bind $(dirname $BAM):$(dirname $dBAM),$(dirname $OUT_BAM):$(dirname $dOUT_BAM),$TMP_DIR:$dTMP_DIR \
      $container samtools"
  else
      echo "Unrecognized container_type $container_type";
      exit 1;
  fi

  #overwrite with container paths
  BAM=$dBAM
  OUT_BAM=$dOUT_BAM
  TMP_DIR=$dTMP_DIR
else
  cmd_samtools=samtools
fi

$cmd_samtools sort -T $TMP_DIR -o $OUT_BAM $BAM
$cmd_samtools index $OUT_BAM

if [ -f ${lOUT_BAM}.bai ]; then rm $lBAM; fi
echo $lOUT_BAM

