#!/bin/sh

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -b|--bam) BAM="$2"; shift ;;
        -s|--chrSizes) CHR_SIZES="$2"; shift ;;
        -o|--outDir) O="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ -z $BAM ] || [ -z $CHR_SIZES ] || [ -z $O ]; then
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

mkdir -p $O
O=$(readlink -f $O)

name=${BAM/.bam/""}
name=$(basename $name)
tmpdir=$O/tmp_bam2bw.${name}
mkdir $tmpdir
cd $tmpdir

F_FILE=${name}.factor
if [ -f $F_FILE ]; then
  echo skip factor calc
  FACTOR=$(cat $F_FILE)
else
  echo calc factor
  FACTOR=$(echo "scale=5; 1000000/$(samtools view -c $BAM | awk '{print $1}')" | bc)
  echo $FACTOR > $F_FILE
fi
echo FACTOR is $FACTOR

for strand in unstranded postive negative; do
for norm in raw normalized; do
for splice in show hide; do
  sdir=unstranded
  if [ $strand != unstranded ]; then sdir=stranded; fi
  bwdir=$sdir/$norm
  mkdir -p $bwdir
  suff=bdg
  if [ $splice = "show" ]; then suff=showSplice.bdg; fi
  BDG=$bwdir/${name}_${norm}_${strand}.${suff}
  echo make bedgraph $BDG
  sarg=""
  if [ $strand = positive ]; then 
    sarg="-strand +";
  elif [ $strand = negative ]; then
    sarg="-strand -";
  fi
  farg=""
  if [ $norm = normalized ]; then
    farg="-scale $FACTOR"
  fi
  sparg=""
    if [ $norm = normalized ]; then
    sparg="-split"
  fi

  if [ -f $BDG ]; then echo skip $BDG, delete to rerun; else genomeCoverageBed -bg $sparg $farg $sarg -ibam $BAM -g $CHR_SIZES > $BDG; bedSort $BDG $BDG; fi
  BW=${BDG/%.bdg/.bw}
  echo make bigwig $BW
  if [ -f $BW ]; then echo skip bigwig $BW, detete to rerun; else bedGraphToBigWig $BDG $CHR_SIZES $BW; fi
  ln $BW ..
done; done; done

cd ..
echo delete $tmpdir if bigwigs look good
rm -r $tmpdir
#rm ${cell_line}_${extraction}_${replicate}_normalized.bedgraph

