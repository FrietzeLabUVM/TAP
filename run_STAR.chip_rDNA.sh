#!/bin/bash
#SBATCH --nodes=1                               # Request one core
#SBATCH --ntasks-per-node=1                               # Request one node (if you request more than one core with -n, also using
#SBATCH --cpus-per-task=12                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 1-06:00                         # Runtime in D-HH:MM format
#SBATCH -p bluemoon                           # Partition to run in
#SBATCH --mem=32000                        # Memory total in MB (for all cores)
#SBATCH -o star_pe_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e star_pe_%j.err                 # File to which STDERR will be written, including job ID

#alignment script developed for dbgap data, uses ENCODE standards
#assumes data is paired end fastq.gz and files end in _R1_001.fastq.gz

echo $0 $@

mode=PE

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f1|--fastq1) F1="$2"; shift ;;
        -f2|--fastq2) F2="$2"; shift ;;
        -f1s|--f1_suffix) F1_suff="$2"; shift ;;
        -f2s|--f2_suffix) F2_suff="$2"; shift ;;
        -o|--out) B="$2"; shift ;;
        -wd|--workdir) wd="$2"; shift ;;
        -idx|--starindex) star_idx="$2"; shift ;;
        -docker|--docker) docker="$2"; shift ;;
        -SE|--SE) mode=SE ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ -z $F1_suff ]; then F1_suff=_R1_001.fastq.gz; fi
if [ -z $F2_suff ]; then F2_suff=_R2_001.fastq.gz; fi

if [ -z $wd ]; then wd=$(dirname $F1)/alignment; fi
mkdir -p $wd

if [ -z $F1 ]; then echo "expect read1 fastq as arg1, stop"; exit 1; fi
if [ -z $F2 ]; then F2=${F1//$F1_suff/$F2_suff}; fi

F1=${F1//"&"/" "}
for f in $F1; do if [ ! -f $f ]; then echo fastq1 $f could not be found! quit; exit 1; fi; done
F1=${F1//" "/","}

if [ $mode = PE ]; then
  F2=${F2//"&"/" "}
  for f in $F2; do 
    if [ ! -f $f ]; then 
      echo fastq2 $f could not be found! Are you supplying as -f2 or deriving from read1? quit; 
      exit 1; 
    fi; 
  done
  F2=${F2//" "/","}
else
  F2=""
fi

if [ -z $B ]; then B=$(basename ${F1/$F1_suff/""}); fi
B=$wd/$(basename $B)


echo read1 fastq is $F1
if [ $mode = PE ]; then echo read2 fastq is $F2; fi
echo output prefix is $B

if [ ! -d $(dirname $B) ]; then
  mkdir -p $(dirname $B)
fi

#if [ -f ${B}.Aligned.sortedByCoord.out.bam ]; then # && [ -f ${B}.Aligned.toTranscriptome.out.bam ]; then
# echo output bam ${B}.Aligned.sortedByCoord.out.bam exists! will not rerun.
# exit 0
#fi

#if [ -f ${B}.Aligned.out.bam ]; then # && [ -f ${B}.Aligned.toTranscriptome.out.bam ]; then
# echo output bam ${B}.Aligned.out.bam exists! will not rerun.
# exit 0
#fi


echo STAR index is $star_idx

# docker RNA v1.0
echo docker is $docker
if [ ! -z $docker ]; then
  # Derive mount points for input files inside docker
  # Fastq files may be "," delim lists of fastqs. This requires special handling.
  F1=${F1//","/" "}
  F2=${F2//","/" "}
  fastq_dir=$(for f in $F1 $F2; do echo $(dirname $f); done | sort | uniq)
  arr=($fastq_dir)
  num_uniq=${#arr[@]}
  if [ $num_uniq != 1 ]; then 
    echo "For docker usage, all fastqs must be in the same directory! Found $num_uniq different directories. Quit!"
    exit 1;
  fi
  F1=${F1//" "/","}
  F2=${F2//" "/","}
  dF=/input
  dF1=${F1//$fastq_dir/"/input"}
  dF2=${F2//$fastq_dir/"/input"}
  dB=/output/$(basename $B)
  dstar_idx=/reference/$(basename $star_idx)

  echo docker F1 is $dF1
  echo docker F2 is $dF2
  echo docker B is $dB
  echo docker star_idx is $dstar_idx

  cmd="docker run \
    -u $(id -u):$(id -g) \
    -v $fastq_dir:/input \
    -v $(dirname $B):$(dirname $dB) \
    -v $star_idx:$dstar_idx \
    --entrypoint STAR\
    $docker \
    "

  #update STAR command inputs to use docker paths
  F1=$dF1
  F2=$dF2
  B=$dB
  star_idx=$dstar_idx
  #echo $cmd
else 
  cmd=STAR
fi

cmd="$cmd \
  --runThreadN 12 \
  --readFilesIn $F1 $F2 \
  --readFilesCommand gunzip -c \
  --genomeDir $star_idx \
  --alignIntronMax 1 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --alignMatesGapMax 1000000 \
  --outFileNamePrefix ${B}. `#output params` \
  --outSAMtype BAM Unsorted
"
echo "$cmd"
$cmd
