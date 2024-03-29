#!/bin/bash
set -e
libType=SE

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -b|--bam) BAM="$2"; shift ;;
        -s|--chrSizes) CHR_SIZES="$2"; shift ;;
        -o|--outDir) OUT_DIR="$2"; shift ;;
        -pe|--pe) libType=PE ;;
        -docker|--docker) container="$2"; shift ;;
        -singularity|--singularity) container="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ -z $BAM ] || [ -z $CHR_SIZES ] || [ -z $OUT_DIR ]; then
  echo bam file \(-b\|--bam\), chrSizes file \(-s\|--chrSizes\), and output directory \(-o\|--outDir\) are all required. quit.
  exit 1
fi

if [ ! -f $BAM ]; then echo bam file $BAM not found! quit.; exit 1; fi
BAM=$(readlink -f $BAM)
echo BAM is $BAM
if [ ! -f $BAM ]; then
  echo BAM does not exist, stop.
  echo add BAM= in qsub or supply as first arg.
  exit 1
fi

if [ ! -f $CHR_SIZES ]; then echo chrSizes file $CHR_SIZES not found! quit.; exit 1; fi
CHR_SIZES=$(readlink -f $CHR_SIZES)
echo CHR_SIZES is $CHR_SIZES
if [ ! -f $CHR_SIZES ]; then
  echo CHR_SIZES does not exist, stop.
  echo add CHR_SIZES= in qsub or supply as second arg
  exit 1
fi

mkdir -p $OUT_DIR
OUT_DIR=$(readlink -f $OUT_DIR)
echo output directory is $OUT_DIR
name=${BAM/.bam/""}
name=$(basename $name)
tmpdir=$OUT_DIR/tmp_bam2bw.${name}

echo TMPDIR is $tmpdir

mkdir -p $tmpdir

OUT_DIR_local=$OUT_DIR

echo container is $container
container_type=""
if [ -n "$container" ]; then container_type="docker"; fi
if [[ "$container" == *.sif ]]; then container_type="singularity"; fi
echo container_type is $container_type

if [ -n "$container" ]; then
  dBAM=/input_bam/$(basename $BAM)
  dCHR_SIZES=/input_chr_sizes/$(basename $CHR_SIZES)
  dOUT_DIR=/output_bigwigs/$(basename $OUT_DIR)

  if [ $container_type = "docker" ]; then
    base_cmd="docker run \
      -u $(id -u):$(id -g) \
      -v $(dirname $BAM):$(dirname $dBAM) \
      -v $(dirname $CHR_SIZES):$(dirname $dCHR_SIZES) \
      -v $(dirname $OUT_DIR):$(dirname $dOUT_DIR) \
      --entrypoint"
  
    cmd_samtools="$base_cmd samtools $container"
    cmd_genomeCoverageBed="$base_cmd genomeCoverageBed $container"
    cmd_bedGraphToBigWig="$base_cmd bedGraphToBigWig $container"
    cmd_sortBed="$base_cmd sort-bed $container"
  elif [ $container_type = "singularity" ]; then
    base_cmd="singularity exec \
        --bind $(dirname $BAM):$(dirname $dBAM),$(dirname $CHR_SIZES):$(dirname $dCHR_SIZES),$(dirname $OUT_DIR):$(dirname $dOUT_DIR)"
        
    cmd_samtools="$base_cmd $container samtools"
    cmd_genomeCoverageBed="$base_cmd $container genomeCoverageBed"
    cmd_bedGraphToBigWig="$base_cmd $container bedGraphToBigWig"
    cmd_sortBed="$base_cmd $container sort-bed"
  else
      echo "Unrecognized container_type $container_type";
      exit 1;
  fi
  BAM=$dBAM
  CHR_SIZES=$dCHR_SIZES
  OUT_DIR=$dOUT_DIR
else
  cmd_samtools=samtools
  cmd_genomeCoverageBed=genomeCoverageBed
  cmd_bedGraphToBigWig=bedGraphToBigWig
  cmd_sortBed=sort-bed
fi

echo "samtools command is: $cmd_samtools"
echo "genomeCoverageBed command is: $cmd_genomeCoverageBed"
echo "bedGraphToBigWig command is: $cmd_bedGraphToBigWig"
echo "sort-bed command is: $cmd_sortBed"

#for PE need to filter for read 1
if [ libType = PE ]; then
  $cmd_samtools view -hb -f 64 $BAM > read1.bam
  BAM=read1.bam
fi


F_FILE=${OUT_DIR_local}/${name}.factor
if [ -f $F_FILE ]; then
  echo skip factor calc
  FACTOR=$(cat $F_FILE)
else
  echo calc factor
  FACTOR=$(echo "scale=5; 1000000/$($cmd_samtools view -c $BAM | awk '{print $1}')" | bc)
  echo $FACTOR > $F_FILE
fi
echo FACTOR is $FACTOR

for strand in unstranded; do
for norm in raw normalized; do
for splice in hide; do
  sdir=unstranded
  if [ $strand != unstranded ]; then sdir=stranded; fi
  bwdir=$OUT_DIR/$sdir/$norm
  bwdir_local=$OUT_DIR_local/$sdir/$norm
  mkdir -p $bwdir_local

  strand_arg=""
  if [ $strand = positive ]; then 
    strand_arg="-strand +";
  elif [ $strand = negative ]; then
    strand_arg="-strand -";
  fi
  scale_arg=""
  if [ $norm = normalized ]; then
    scale_arg="-scale $FACTOR"
  fi
  suff=bdg
  splice_arg=""
  if [ $splice = "show" ]; then
    suff=showSplice.bdg
  else
    splice_arg="-split"
  fi
    
  BDG=$bwdir/${name}_${norm}_${strand}.${suff}
  BDG_local=$bwdir_local/${name}_${norm}_${strand}.${suff}
  echo make bedgraph $BDG

  cmd=""
  if [ -f $BDG_local ]; then 
    echo skip $BDG_local, delete to rerun; 
  else 
    cmd="$cmd_genomeCoverageBed -bg $splice_arg $scale_arg $strand_arg -ibam $BAM"
    echo $cmd to $BDG_local
    $cmd > $BDG_local
    cmd_sort="$cmd_sortBed --max-mem 19G $BDG"
    echo $cmd_sort to $BDG_local
    $cmd_sort > ${BDG_local}.tmp
    mv ${BDG_local}.tmp $BDG_local
    #$cmd_genomeCoverageBed -bg $splice_arg $scale_arg $strand_arg -ibam $BAM > $BDG
    #$cmd_sortBed -i $BDG > $BDG
  fi

  BW=${BDG/%.bdg/.bw}
  BW_local=${BDG_local/%.bdg/.bw}
  echo make bigwig $BW_local
  cmd2=""
  if [ -f $BW_local ]; then 
    echo skip bigwig $BW_local, delete to rerun
  else 
    cmd2="$cmd_bedGraphToBigWig $BDG $CHR_SIZES $BW"
    echo $cmd2
    $cmd2
    #$cmd_bedGraphToBigWig $BDG $CHR_SIZES $BW
  fi

  mv "$BW_local" "$OUT_DIR_local"/
done; done; done

cd ..
#echo delete "$tmpdir" if bigwigs look good
rm -r $tmpdir
#rm ${cell_line}_${extraction}_${replicate}_normalized.bedgraph
echo FINISHED
