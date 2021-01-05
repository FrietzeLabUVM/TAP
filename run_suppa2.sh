#!/bin/bash
#$ -N suppa2
#py=/slipstream/home/dbgap/anaconda3/bin/python3

echo $0 $@

SUPPA_BIN=~/lab_bin

suppa="$SUPPA_BIN/python $SUPPA_BIN/suppa.py"
salmon_out=$1
gtf=$2
suppa_ref=$3

echo salmon_out is $salmon_out
echo gtf is $gtf
echo suppar_ref is $suppa_ref

if [ -z $salmon_out ]; then echo must have salmon quant result directory as arg1. quit.; exit 1; fi
if [ ! -d $salmon_out ]; then echo could not find salmon quant directory $salmon_out. quit.; exit 1; fi
raw=$salmon_out/quant.sf
tpm=$salmon_out/tpm.txt
if [ ! -f $raw ]; then echo supplied salmon quant directory, $salmon_out, is not valid.  must contain quant.sf. quit.; exit 1; fi

#create a tpm file for suppa2
echo $(basename $salmon_out ) > $tpm
cut -f 1,4 $raw | awk 'NR > 1' >> $tpm

echo running isoforms
$suppa psiPerIsoform -g ${gtf} -e $tpm -o $salmon_out/suppa2


for ioe in $suppa_ref/*ioe; do
  name=$(basename $ioe)
  name=${name/$(basename $gtf)"._"/""}
  name=${name/suppa2./""}
  name=${name/_strict.ioe/""}
  #name=$(dirname $salmon_out)
  echo running $name
  $suppa psiPerEvent -i ${ioe} -e $tpm -o $salmon_out/suppa2.alternate_events.${name}
done

