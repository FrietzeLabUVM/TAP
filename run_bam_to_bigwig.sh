#!/bin/bash
#SBATCH --mem=16000
#SBATCH -o bigwig_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e bigwig_%j.err                 # File to which STDERR will be written, including job ID

echo script:
echo $0
echo args:
echo $@

libType=SE

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -b|--bam) BAM="$2"; shift ;;
        -s|--chrSizes) CHR_SIZES="$2"; shift ;;
        -o|--outDir) OUT="$2"; shift ;;
        -pe|--pe) libType=PE ;;
        -docker|--docker) docker="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

echo BAM is $BAM
echo CHR_SIZES is $CHR_SIZES
echo OUT is $OUT

if [ -z $BAM ] || [ -z $CHR_SIZES ] || [ -z $OUT ]; then
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

mkdir -p $OUT
O=$(readlink -f $OUT)

name=${BAM/.bam/""}
name=$(basename $name)
tmpdir=$OUT/tmp_bam2bw.${name}

echo TMPDIR is $tmpdir

mkdir -p $tmpdir
cd $tmpdir

# docker for samtools and UCSC tools v1.0
echo docker is $docker
if [ -n $docker ]; then
  dBAM=/input_bam/$(basename $BAM)
  dCHR_SIZES=/input_chr_sizes/$(basename $CHR_SIZES)
  dO=/output_bigwigs/$(basename $O)

  base_cmd="docker run \
    -u $(id -u):$(id -g) \
    -v $(dirname $BAM):$(dirname $dBAM) \
    -v $(dirname $CHR_SIZES):$(dirname $dCHR_SIZES) \
    -v $(dirname $O):$(dirname $dO) \
    --entrypoint"
    
  cmd_samtools="$base_cmd samtools $docker"
  cmd_genomeCoverageBed="$base_cmd genomeCoverageBed $docker"
  cmd_bedGraphToBigWig="$base_cmd bedGraphToBigWig $docker"
  cmd_bedSort="$base_cmd bedSort $docker"

  BAM=$dBAM
  CHR_SIZES=$dCHR_SIZES
  O=$dO
else
  cmd_samtools=samtools
  cmd_genomeCoverageBed=genomeCoverageBed
  cmd_bedGraphToBigWig=bedGraphToBigWig
  cmd_bedSort=bedSort
fi

#for PE need to filter for read 1
if [ libType = PE ]; then
  cmd0="$cmd_samtools view -hb -f 64 $BAM > read1.bam"
  echo running:
  echo $cmd0
  $cmd0
  BAM=read1.bam
fi


F_FILE=${name}.factor
if [ -f $F_FILE ]; then
  echo skip factor calc, read from $F_FILE
  FACTOR=$(cat $F_FILE)
else
  echo calc factor, save to $F_FILE
  echo $cmd_samtools view -c $BAM \| awk \'{print \$1}\'
  FACTOR=$(echo "scale=5; 1000000/$($cmd_samtools view -c $BAM | awk '{print $1}')" | bc)
  echo $FACTOR > $F_FILE
fi
echo FACTOR is $FACTOR

for strand in unstranded positive negative; do
for norm in raw normalized; do
for splice in show hide; do
  sdir=unstranded
  if [ $strand != unstranded ]; then sdir=stranded; fi
  bwdir=$sdir/$norm
  mkdir -p $bwdir

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
  echo make bedgraph $BDG
  cmd=""
  if [ -f $BDG ]; then 
    echo skip $BDG, delete to rerun; 
  else 
    cmd1="$cmd_genomeCoverageBed -bg $splice_arg $scale_arg $strand_arg -ibam $BAM -g $CHR_SIZES > $BDG"; 
    cmd2="$cmd_bedSort $BDG $BDG"
    echo running: 
    echo $cmd1
    $cmd1
    echo $cmd2
    $cmd2 
$cmd_bedSort $BDG $BDG; fi

  BW=${BDG/%.bdg/.bw}
  echo make bigwig $BW
  cmd3=""
  if [ -f $BW ]; then 
    echo skip bigwig $BW, delete to rerun; 
  else 
    cmd3="$cmd_bedGraphToBigWig $BDG $CHR_SIZES $BW"; 
    echo $cmd3; 
    $cmd3; 
  fi

  ln $BW ..
done; done; done

cd ..
echo delete $tmpdir if bigwigs look good
#rm -r $tmpdir
#rm ${cell_line}_${extraction}_${replicate}_normalized.bedgraph

