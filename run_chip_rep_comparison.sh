#!/bin/bash
#$ -o idr."$JOB_ID".out
#$ -e idr."$JOB_ID".error

POOLNAME=$1

PEAKS=$2
arr=( ${PEAKS//,/" "} )
arrLen=${#arr[@]}

BAMS=$3

CHR_SIZES=$4

OUTDIR=$5
if [ -z $OUTDIR ]; then
  OUTDIR=$(pwd)
fi
if [ ! -d $OUTDIR ]; then
  mkdir -p $OUTDIR
fi

SCRIPT_DIR=idr_package

i=0; 
while [ $i -lt $arrLen ]; do 
  i2=$(( $i + 1 )); 
  while [ $i2 -lt $arrLen ]; do 
    echo $i $i2
	PEAKS1=${arr[$i]}
	PEAKS2=${arr[$i2]}
        echo PEAKS1 is $PEAKS1
        echo PEAKS2 is $PEAKS2
	if [ -z $PEAKS1 ] || [ ! -f $PEAKS1 ]; then
		echo PEAKS1 $PEAKS1 not found! stop
		exit 1
	fi
	if [ -z $PEAKS2 ] || [ ! -f $PEAKS2 ]; then
		echo PEAKS2 $PEAKS2 not found! stop
		exit 1
	fi
	
	PREFIX=${POOLNAME}.$(basename $PEAKS1)_vs_$(basename $PEAKS2)_IDR
	if [ -f $PREFIX.npeaks-aboveIDR.txt ]; then
		echo file $PREFIX.npeaks-aboveIDR.txt exists so skip IDR.  deleted this file to run.
		exit 0
	fi
	echo PREFIX is $PREFIX
	#num of peaks is capped at 250k, sort by pval, head 250k, resort by chrm position
	#cat output/MCF10A_H3K4AC_R1/MCF10A_H3K4AC_R1_loose_peaks.narrowPeak | sort -n -r -k 8 | head
	#sort -nr -k 8,8 $peaks | head -n $thresh | sort -n -k 2 | sort -k 1,1V -k 2,2n > $trunc
	MAX=250000
	np1=$(wc -l $PEAKS1 | cut -d " " -f 1)
	np2=$(wc -l $PEAKS2 | cut -d " " -f 1)
	IDR_IN1=$PEAKS1
	IDR_IN2=$PEAKS2
	if [ $np1 -gt $MAX ]; then #if > MAX peaks, sort by pval and cut then resort by chrm
		IDR_IN1="$PEAKS1".idr_in
		echo $(basename $PEAKS1) too many peaks. $np1.  filter top $MAX to $(basename $IDR_IN1).
		cat $PEAKS1 | sort -nr -k 8,8 | head -n $MAX | sort -n -k 2 | sort -k 1,1V -k 2,2n > $IDR_IN1
	fi
	if [ $np2 -gt $MAX ]; then #if > MAX peaks, sort by pval and cut then resort by chrm
		IDR_IN2="$PEAKS2".idr_in
		echo $(basename $PEAKS2) too many peaks. $np2.  filter top $MAX to $(basename $IDR_IN2).
		cat $PEAKS2 | sort -nr -k 8,8 | head -n $MAX | sort -n -k 2 | sort -k 1,1V -k 2,2n > $IDR_IN2
	fi
	isBroad=F
	Rscript $SCRIPT_DIR/batch-consistency-analysis.r $SCRIPT_DIR $IDR_IN1 $IDR_IN2 -1 0 $isBroad p.value $CHR_SIZES ${PREFIX}.Rout.txt ${PREFIX}.overlapped-peaks.txt ${PREFIX}.npeaks-aboveIDR.txt ${PREFIX}.em.sav ${PREFIX}.uri.sav
	
	i2=$(( $i2 + 1 ))
  done; 
  i=$(( $i + 1 )); 
done

#exit 0

### bamcorrelate
#!/bin/bash
#$ -o bamcorr."$JOB_ID".out
#$ -e bamcorr."$JOB_ID".error
#$ -pe threads 8

#pools to input semicolon delimited list of bam files into the output file
#arg 1 is semicolon separated list of bam files
#arg 2 is output dir
#arg 3 is short descriptor appended to outputs


#called with qsub_out=$(qsub -v BAMS=$BAMS,OUTDIR=$OUTDIR,DESC=global -wd $OUTDIR -hold_jid $DEP_JIDS job_scripts/run_bamCorrelate.sh)


#if [ -z $REF ]; then
#REF=$3
#fi
echo BAMS is $BAMS
echo OUTDIR is $OUTDIR
echo REF is $REF
bams=${BAMS//";"/" "}
bams_array=( $bams )
for b in ${bams_array[@]}; do
	if [ ! -f $b ]; then echo bam file $b not found! stop.; exit 1; fi
done
if [ ! -d $OUTDIR ]; then echo OUTDIR $OUTDIR not found! stop.; exit 1; fi

#rand=$(date +%N | sed -e 's/000$//' -e 's/^0//')
#bamCorrelate min args are "bamCorrelate bins ---bamfiles file1.bam file2.bam --corMethod spearman -o heatmap.png"
#cmd="bamCorrelate bins --bamfiles $bams --corMethod spearman -o $OUTDIR/bamCorrelate_heatmap_"$rand".pdf  --outFileCorMatrix $OUTDIR/bamCorrelate_values_"$rand".txt --plotFileFormat pdf"
outf=global
if [ ! -z $REF ]; then outf=$(basename $REF | awk 'BEGIN {FS="."} {print $1} '); fi
if [ -f $OUTDIR/bamCorrelate_"$outf"_values.txt ]; then
	echo bamCorrelate output exists already! 
	echo delete $OUTDIR/bamCorrelate_"$outf"_values.txt to rerun
	exit 1
fi
if [ -z $REF ]; then
	echo no REF set, running bamCorrelate globally
	cmd="bamCorrelate bins --bamfiles $bams --corMethod spearman -o $OUTDIR/bamCorrelate_"$outf"_heatmap.pdf  --outFileCorMatrix $OUTDIR/bamCorrelate_"$outf"_values.txt --plotFileFormat pdf"
else
	if [ ! -f $REF ]; then echo REF $REF not found, stop.; exit 1; fi
	cmd="bamCorrelate BED-file --BED $REF --bamfiles $bams --corMethod spearman -o $OUTDIR/bamCorrelate_"$outf"_heatmap.pdf  --outFileCorMatrix $OUTDIR/bamCorrelate_"$outf"_values.txt --plotFileFormat pdf"
fi
echo CMD is "$cmd"
$cmd

for b in ${bams_array[@]}; do
echo $b
done
echo bamCorrelate finished
#echo "$INPUT" > $pooled
