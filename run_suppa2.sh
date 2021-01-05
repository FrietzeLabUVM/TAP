#!/bin/bash
#$ -N suppa2
#py=/slipstream/home/dbgap/anaconda3/bin/python3

SUPPA_BIN=~/lab_bin

suppa="$SUPPA_BIN/python $SUPPA_BIN/suppa.py"
salmon_out=$1
gtf=$2
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


for ioe in ~/annotation/hg38/suppa2/gencode.v35.annotation.suppa2.*ioe; do
  name=$(basename $ioe)
  name=${name/gencode.v35.annotation.suppa2.ioe_/""}
  name=${name/_strict.ioe/""}
  echo running $name
  $suppa psiPerEvent -i ${ioe} -e $tpm -o $salmon_out/suppa2.alternate_events.${name}
done

