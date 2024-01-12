#!/bin/bash
set -e
### runs MACS2 (tight, loose, and broad) generates FE bigwigs
BDG="--bdg"
extra=""
read_mode=SE

echo script:
echo $0
echo args:
echo $@

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -t|--treat) TREAT_BAM="$2"; shift ;;
        -i|--input) INPUT_BAM="$2"; shift ;;
        -o|--out_dir) OUT_DIR="$2"; shift ;;
        -p|--prefix) PREFIX="$2"; shift ;;
        -pval_narrow|--pval_narrow) PVAL="$2"; shift ;;
        -qval_narrow|--qval_narrow) QVAL="$2"; shift ;;
        -pval_loose|--pval_loose) PVAL_LOOSE="$2"; shift ;;
        -qval_loose|--qval_loose) QVAL_LOOSE="$2"; shift ;;
        -pval_broad|--pval_broad) PVAL_BROAD="$2"; shift ;;
        -qval_broad|--qval_broad) QVAL_BROAD="$2"; shift ;;
        -broad_cutoff|--broad_cutoff) BROADCUTOFF="$2"; shift ;;
        -PE|--PE) read_mode=PE ;;
        -SE|--SE) read_mode=SE ;;
        -g|--gen) GEN="$2"; shift ;;
        -s|--chrSizes) CHR_SIZES="$2"; shift ;;
        -noModel|--noModel) extra="--nomodel --extsize 147" ;;
        -no_bdg|--no_bdg) BDG="";;
        -docker|--docker) container="$2"; shift ;;
        -singularity|--singularity) container="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done



echo TREAT_BAM is $TREAT_BAM
echo INPUT_BAM is $INPUT_BAM
echo OUT_DIR is $OUT_DIR
echo PREFIX is $PREFIX
echo PVAL is $PVAL
echo QVAL is $QVAL
echo PVAL is $PVAL_LOOSE
echo QVAL is $QVAL_LOOSE
echo PVAL is $PVAL_BROAD
echo QVAL is $QVAL_BROAD
echo BROADCUTOFF is $BROADCUTOFF
echo read_mode is $read_mode
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

if [ ! -z $OUT_DIR ]; then
  mkdir -p $OUT_DIR
fi

if [ -z $OUT_DIR ]; then
        OUT_DIR=$(pwd)
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
    stat_val="5e-3"
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
  if [ -z $QVAL_BROAD ]; then
    stat_broad="q"
    stat_val_broad="5e-3"
  else
    stat_broad="q"
    stat_val_broad=$QVAL_BROAD
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
  if [ -z $QVAL_LOOSE ]; then
    stat_loose="p"
    stat_val_loose="2e-2"
  else
    stat_loose="q"
    stat_val_loose=$QVAL_LOOSE
  fi
else
  stat_loose="p"
  stat_val_loose=$PVAL_LOOSE
fi

g=$(cat $CHR_SIZES | awk -v  total=0 '{total = total + $2}; END {print total}')

OUT_DIR_local=$OUT_DIR

# container for macs2 v1.1
echo container is $container
container_type=""
if [ -n "$container" ]; then container_type="docker"; fi
if [[ "$container" == *.sif ]]; then container_type="singularity"; fi
echo container_type is $container_type

if [ -n "$container" ]; then
  dTREAT_BAM=/input_treat_bam/$(basename $TREAT_BAM)
  dINPUT_BAM=/input_input_bam/$(basename $INPUT_BAM)
  dCHR_SIZES=/input_chr_sizes/$(basename $CHR_SIZES)
  dOUT_DIR=/output_bigwigs/$(basename $OUT_DIR)

  if [ $container_type = "docker" ]; then
    base_cmd="docker run \
      -u $(id -u):$(id -g) \
      -v $(dirname $TREAT_BAM):$(dirname $dTREAT_BAM) \
      -v $(dirname $INPUT_BAM):$(dirname $dINPUT_BAM) \
      -v $(dirname $CHR_SIZES):$(dirname $dCHR_SIZES) \
      -v $OUT_DIR:$dOUT_DIR \
      --entrypoint"
      
    cmd_macs2="$base_cmd macs2 $container"
    cmd_bedGraphToBigWig="$base_cmd bedGraphToBigWig $container"
    cmd_sortBed="$base_cmd sort-bed $container"
  elif [ $container_type = "singularity" ]; then
    base_cmd="singularity exec \
        --bind $(dirname $TREAT_BAM):$(dirname $dTREAT_BAM),$(dirname $INPUT_BAM):$(dirname $dINPUT_BAM),$(dirname $CHR_SIZES):$(dirname $dCHR_SIZES),$OUT_DIR:$dOUT_DIR"
        
    cmd_macs2="$base_cmd $container macs2"
    cmd_bedGraphToBigWig="$base_cmd $container bedGraphToBigWig"
    cmd_sortBed="$base_cmd $container sort-bed"
  else
      echo "Unrecognized container_type $container_type";
      exit 1;
  fi

  TREAT_BAM=$dTREAT_BAM
  INPUT_BAM=$dINPUT_BAM
  CHR_SIZES=$dCHR_SIZES
  OUT_DIR=$dOUT_DIR
else
  cmd_macs2=macs2
  cmd_bedGraphToBigWig=bedGraphToBigWig
  cmd_sortBed=sort-bed
fi

#narrow tight
if [ -f $OUT_DIR/$PREFIX"_peaks.narrowPeak" ]; then
  echo $PREFIX"_peaks.narrowPeak" exists. macs2 has already been run for $TREAT_BAM.
  echo delete $PREFIX"_peaks.narrowPeak" if you want to rerun macs2.
else
  cmd_tight="$cmd_macs2 callpeak -t $TREAT_BAM -c $INPUT_BAM -g $g --outdir $OUT_DIR -n $PREFIX -$stat $stat_val --bdg $extra"
  echo cmd_tight is:
  echo $cmd_tight
  $cmd_tight
fi

#narrow loose
if [ -f $OUT_DIR/$PREFIX"_loose_peaks.narrowPeak" ]; then
  echo $PREFIX"_loose_peaks.narrowPeak" exists. macs2 has already been run for $TREAT_BAM.
  echo delete $PREFIX"_loose_peaks.narrowPeak" if you want to rerun macs2.
else
  cmd_loose="$cmd_macs2 callpeak -t $TREAT_BAM -c $INPUT_BAM -g $g --outdir $OUT_DIR -n ${PREFIX}_loose -$stat_loose $stat_val_loose $extra"
  echo cmd_loose is:
  echo $cmd_loose
  $cmd_loose

fi

#broad
if [ -f $OUT_DIR/$PREFIX"_peaks.broadPeak" ]; then
  echo $PREFIX"_peaks.broadPeak" exists. macs2 has already been run for $TREAT_BAM.
  echo delete $PREFIX"_peaks.broadPeak" if you want to rerun macs2.
else
  cmd_broad="$cmd_macs2 callpeak -t $TREAT_BAM -c $INPUT_BAM -g $g --outdir $OUT_DIR -n $PREFIX -$stat_broad $stat_val_broad --broad --broad-cutoff $BROADCUTOFF $extra"
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
	TREATMENT=$OUT_DIR/$PREFIX"_treat_pileup.bdg"
	CONTROL=$OUT_DIR/$PREFIX"_control_lambda.bdg"
	METHOD=FE

	echo treat file - $TREATMENT
	echo control file - $CONTROL
	echo method - $METHOD
	CMP_BDG=$OUT_DIR/$(basename $TREATMENT)
	CMP_BDG="${CMP_BDG/_treat_pileup.bdg/}"_"$METHOD".bdg
  CMP_BW=${CMP_BDG/.bdg/.bw}

  CMP_BDG_local=$OUT_DIR_local/$(basename $TREATMENT)
	CMP_BDG_local="${CMP_BDG_local/_treat_pileup.bdg/}"_"$METHOD".bdg
  CMP_BW_local=${CMP_BDG_local/.bdg/.bw}

	if [ -f $CMP_BDG_local ] || [ -f ${CMP_BW_local} ]; then
		echo skipping bdgcmp for "$TREATMENT", file "$CMP_BDG_local" or $CMP_BW_local exists
	else
    run_bdgcmp="$cmd_macs2 bdgcmp -t $TREATMENT -c $CONTROL -m $METHOD -o $CMP_BDG"
    echo run_bdgcmp is:
    echo $run_bdgcmp
		$run_bdgcmp
    run_sortBed="$cmd_sortBed --max-mem 19G $CMP_BDG"
    echo run_sortBed is:
    echo $run_sortBed to $CMP_BDG_local
    $run_sortBed > ${CMP_BDG_local}.tmp
    mv ${CMP_BDG_local}.tmp ${CMP_BDG_local}

	fi

	#required inputs:
	echo input - $CMP_BDG
	echo chrSizes - $CHR_SIZES
	echo output - $CMP_BW
	if [ -f $CMP_BW_local ]; then
		echo file $CMP_BW_local exists so bdg2bw not necessary for $inputBegGraph
		echo nothing done.
	else
    run_bedGraphToBigWig="$cmd_bedGraphToBigWig $CMP_BDG $CHR_SIZES $CMP_BW"
    echo run_bedGraphToBigWig is:
    echo $run_bedGraphToBigWig
    $run_bedGraphToBigWig
	fi
fi
echo FINISHED

