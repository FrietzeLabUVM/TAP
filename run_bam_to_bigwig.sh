#/bin/sh

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

BDGp=${name}_norm_pos.bdg
BDGn=${name}_norm_neg.bdg
BDGp_showSplice=${name}_norm_pos.showSplice.bdg
BDGn_showSplice=${name}_norm_neg.showSplice.bdg

echo make bedgraph
if [ -f $BDGp ]; then echo skip $BDGp, delete to rerun; else genomeCoverageBed -bg -split -scale $FACTOR -strand - -ibam $sBAM -g $CHR_SIZES > $BDGp; fi
if [ -f $BDGn ]; then echo skip $BDGn, delete to rerun; else genomeCoverageBed -bg -split -scale $FACTOR -strand + -ibam $sBAM -g $CHR_SIZES > $BDGn; fi
if [ -f $BDGp_showSplice ]; then echo skip $BDGp_showSplice, delete to rerun; else genomeCoverageBed -bg -scale $FACTOR -strand - -ibam $sBAM -g $CHR_SIZES > $BDGp_showSplice; fi
if [ -f $BDGn_showSplice ]; then echo skip $BDGn_showSplice, delete to rerun; else genomeCoverageBed -bg -scale $FACTOR -strand + -ibam $sBAM -g $CHR_SIZES > $BDGn_showSplice; fi

BWp=${name}_norm_pos.bw
BWn=${name}_norm_neg.bw
BWp_showSplice=${name}_norm_pos.showSplice.bw
BWn_showSplice=${name}_norm_neg.showSplice.bw

if [ -f $BWp ] && [ -f $BWn ] && [ -f $BWp_showSplice ] && [ -f $BWn_showSplice ]; then
  echo skip making bigwig
else
  echo make bigwig
  bedGraphToBigWig $BDGp $CHR_SIZES $BWp
  bedGraphToBigWig $BDGn $CHR_SIZES $BWn
  bedGraphToBigWig $BDGp_showSplice $CHR_SIZES $BWp_showSplice
  bedGraphToBigWig $BDGn_showSplice $CHR_SIZES $BWn_showSplice

fi

mv $BWp ..
mv $BWn ..
mv $BWp_showSplice ..
mv $BWn_showSplice ..
cd ..
echo delete $tmpdir if $BWp and $BWn and $BWp_showSplice and $BWn_showSplice are correct
#rm ${cell_line}_${extraction}_${replicate}_normalized.bedgraph

