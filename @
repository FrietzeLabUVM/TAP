#!/bin/bash
SUPPA_BIN=~/lab_bin
suppa="$SUPPA_BIN/python $SUPPA_BIN/suppa.py"

gtf=$1
dest=$2
if [ -z $gtf ]; then echo need gtf at arg1! exit; exit 1; fi
if [ ! -f $gtf ]; then echo gtf $gtf not found! exit; exit 1; fi
if [ -z $dest ]; then echo output will be in current directory; dest=$(pwd); fi
o=${dest}/"suppa2."$(basename ${gtf})"."

#use -p for ensembl and UCSC, not gencode
$suppa generateEvents -i $gtf -o $o  -f ioi -p


for e in SE SS MX RI FL; do
  $suppa generateEvents -i $gtf -o $o  -f ioe -e $e -p
done
