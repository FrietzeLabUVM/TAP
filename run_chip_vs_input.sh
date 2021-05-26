#!/bin/bash

### runs MACS2 (tight, loose, and broad) generates FE bigwigs

BDG="--bdg"
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -t|--treat) TREAT_BAM="$2"; shift ;;
        -i|--input) INPUT_BAM="$2"; shift ;;
        -o|--out_dir) OUTDIR="$2"; shift ;;
        -p|--prefix) PREFIX="$2"; shift ;;
        -pval_narrow|--pval_narrow) PVAL="$2"; shift ;;
        -qval_narrow|--qval_narrow) QVAL="$2"; shift ;;
	-pval_loose|--pval_loose) PVAL_LOOSE="$2"; shift ;;
        -qval_loose|--qval_loose) QVAL_LOOSE="$2"; shift ;;
	-pval_broad|--pval_broad) PVAL_BROAD="$2"; shift ;;
        -qval_broad|--qval_broad) QVAL_BROAD="$2"; shift ;;
	-broad_cutoff|--broad_cutoff) BROADCUTOFF="$2"; shift ;;
        -g|--gen) GEN="$2"; shift ;;
	-s|--chrSizes) CHR_SIZES="$2"; shift ;;
        -no_bdg|--no_bdg) BDG="";;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done



echo TREAT_BAM is $TREAT_BAM
echo INPUT_BAM is $INPUT_BAM
echo OUTDIR is $OUTDIR
echo PREFIX is $PREFIX
echo PVAL is $PVAL
echo QVAL is $QVAL
echo PVAL is $PVAL_LOOSE
echo QVAL is $QVAL_LOOSE
echo PVAL is $PVAL_BROAD
echo QVAL is $QVAL_BROAD
echo BROADCUTOFF is $BROADCUTOFF
echo GEN is $GEN
echo CHR_SIZES is $CHR_SIZES
echo macs2 on $(basename $TREAT_BAM):$(basename $INPUT_BAM)

if [ -z "$TREAT_BAM" ] #|| [ ! -e $TREAT_BAM ]
then
	echo TREAT_BAM $TREAT_BAM arg required! stop
	exit 1
fi
for tb in $TREAT_BAM; do
  if [ ! -e $tb ]; then
    echo TREAT_BAM $tb not found! stop
    exit 1
  fi
done

if [ -z "$INPUT_BAM" ] #|| [ ! -e $INPUT_BAM ]
then
	echo INPUT_BAM $INPUT_BAM arg required! stop
	exit 1
fi
for ib in $INPUT_BAM; do
  if [ ! -e $ib ]; then
    echo INPUT_BAM $ib not found! stop
    exit 1
  fi
done

if [ ! -z $OUTDIR ]; then
  mkdir -p $OUTDIR
fi

if [ -z $OUTDIR ]; then
        OUTDIR=$(pwd)
	echo using current directory as output location
#	exit 1
fi
if [ -z $PREFIX ]; then
	echo PREFIX $PREFIX missing! stop
	exit 1
fi
#tight peak cutoff
if [ -z $PVAL ]; then
  if [ -z $QVAL ]; then
    stat="q"
    stat_val="1e-2"
  else
    stat="q"
    stat_val=$QVAL
  fi
else
  stat="p"
  stat_val=$PVAL
fi
#broad peak cutoff
if [ -z $PVAL_BROAD ]; then
  if [ -z $PVAL_BROAD ]; then
    stat_broad="q"
    stat_val_broad="1e-2"
  else
    stat_broad="q"
    stat_val_broad=$PVAL_BROAD
  fi
else
  stat_broad="p"
  stat_val_broad=$PVAL_BROAD
fi
#broadcutuff
if [ -z $BROADCUTOFF ]; then
  BROADCUTOFF="1e-1"
fi
#loose peak cutoff
if [ -z $PVAL_LOOSE ]; then
  if [ -z $PVAL_LOOSE ]; then
    stat_loose="p"
    stat_val_loose="1e-2"
  else
    stat_loose="q"
    stat_val_loose=$PVAL_LOOSE
  fi
else
  stat_loose="p"
  stat_val_loose=$PVAL_LOOSE
fi


if [ $GEN = "mm10" ] || [ $GEN = "mm9" ] || [ $GEN = "MM9" ] || [ $GEN = "MM10" ]; then
        g=mm
elif [ $GEN = "hg38" ] || [ $GEN = "hg19" ] || [ $GEN = "U13369" ] || [ $GEN = "HG38" ] || [ $GEN = "HG19" ]; then
        g=hs
else
        echo unrecognized genome ${g} !; exit 1
fi

extra="--nomodel --extsize 147"

#narrow tight
if [ -f $OUTDIR/$PREFIX"_peaks.narrowPeak" ]; then
  echo $PREFIX"_peaks.narrowPeak" exists. macs2 has already been run for $TREAT_BAM.
  echo delete $PREFIX"_peaks.narrowPeak" if you want to rerun macs2.
else
  cmd_tight="macs2 callpeak -t $TREAT_BAM -c $INPUT_BAM -g $g --outdir $OUTDIR -n $PREFIX -$stat $stat_val --bdg $extra"
  echo cmd_tight is:
  echo $cmd_tight
  $cmd_tight
fi

#narrow loose
if [ -f $OUTDIR/$PREFIX"_loose_peaks.narrowPeak" ]; then
  echo $PREFIX"_loose_peaks.narrowPeak" exists. macs2 has already been run for $TREAT_BAM.
  echo delete $PREFIX"_loose_peaks.narrowPeak" if you want to rerun macs2.
else
  cmd_loose="macs2 callpeak -t $TREAT_BAM -c $INPUT_BAM -g $g --outdir $OUTDIR -n ${PREFIX}_loose -$stat $stat_val $extra"
  echo cmd_loose is:
  echo $cmd_loose
  $cmd_loose

fi

#broad
if [ -f $OUTDIR/$PREFIX"_peaks.broadPeak" ]; then
  echo $PREFIX"_peaks.broadPeak" exists. macs2 has already been run for $TREAT_BAM.
  echo delete $PREFIX"_peaks.broadPeak" if you want to rerun macs2.
else
  cmd_broad="macs2 callpeak -t $TREAT_BAM -c $INPUT_BAM -g $g --outdir $OUTDIR -n $PREFIX -$stat $stat_val --broad --broad-cutoff $BROADCUTOFF $extra"
  echo cmd_broad is:
  echo $cmd_broad
  $cmd_broad
fi

if [ $BDG = "--bdg" ]; then
	if [ -z $CHR_SIZES ]; then
		echo "-s|--chrSizes is required for bedgraph comparison"
		exit 1
	fi
	#required inputs:
	#TREATMENT - control lambda bdg file from macs2
	#CONTROL - treatment pileup bdg file from macs2
	#METHOD - should probably be FE, see macs2 bdgcmp help for more options
	TREATMENT=$OUTDIR/$PREFIX"_treat_pileup.bdg"
	CONTROL=$OUTDIR/$PREFIX"_control_lambda.bdg"
	METHOD=FE

	echo treat file - $TREATMENT
	echo control file - $CONTROL
	echo method - $METHOD
	CMP_BDG=$OUTDIR/$(basename $TREATMENT)
	CMP_BDG="${CMP_BDG/_treat_pileup.bdg/}"_"$METHOD".bdg
	if [ -f $CMP_BDG ] || [ -f ${CMP_BDG/.bdg/.bw} ]; then
		echo skipping bdgcmp for "$TREATMENT", file "$CMP_BDG" exists
	else
		macs2 bdgcmp -t $TREATMENT -c $CONTROL -m $METHOD -o $CMP_BDG
                bedSort $CMP_BDG $CMP_BDG
	fi

	CMP_BW=$OUTDIR/$PREFIX"_FE.bw"

	#required inputs:
	echo input - $CMP_BDG
	echo chrSizes - $CHR_SIZES
	echo output - $CMP_BW
	if [ -f $CMP_BW ]; then
		echo file $CMP_BW exists so bdg2bw not necessary for $inputBegGraph
		echo nothing done.
	else
		bedGraphToBigWig $CMP_BDG $CHR_SIZES $CMP_BW
	fi
fi
