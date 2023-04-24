#!/bin/bash
#$ -N suppa2
#py=/slipstream/home/dbgap/anaconda3/bin/python3

echo $0 $@

SALMON_RESULT=$1
GTF=$2
SUPPA_REF=$3
SUPPA_REF_local=$SUPPA_REF
docker=$5

echo SALMON_RESULT is $SALMON_RESULT
echo GTF is $GTF
echo suppa_ref is $SUPPA_REF

if [ -z $SALMON_RESULT ]; then echo must have salmon quant result directory as arg1. quit.; exit 1; fi
if [ ! -d $SALMON_RESULT ]; then echo could not find salmon quant directory $SALMON_RESULT. quit.; exit 1; fi
raw=$SALMON_RESULT/quant.sf
tpm=$SALMON_RESULT/tpm.txt
if [ ! -f $raw ]; then echo supplied salmon quant directory, $SALMON_RESULT, is not valid.  must contain quant.sf. quit.; exit 1; fi

#create a tpm file for suppa2
echo $(basename $SALMON_RESULT ) > $tpm
cut -f 1,4 $raw | awk 'NR > 1' >> $tpm

# docker for suppa2 v1.1
container_type=""
if [ -n "$docker" ]; then container_type="docker"; fi
if [[ "$docker" == *.sif ]]; then container_type="singularity"; fi

echo container_type is $container_type
if [ -n "$container_type" ]; then
  dGTF=/input_gtf/$(basename $GTF)
  dSALMON_RESULT=/output_salmon
  dSUPPA_REF=/suppa_ref

  if [ $container_type = "docker" ]; then
  base_cmd="docker run \
    -u $(id -u):$(id -g) \
    -v $(dirname $GTF):$(dirname $dGTF) \
    -v $SALMON_RESULT:$dSALMON_RESULT \
    -v $SUPPA_REF:$dSUPPA_REF \
    --entrypoint"
  elif [ $container_type = "singularity" ]; then
  base_cmd="singularity exec \
    --bind $(dirname $GTF):$(dirname $dGTF),$SALMON_RESULT:$dSALMON_RESULT,$SUPPA_REF:$dSUPPA_REF \
    $singularity"
  else
    echo "Unrecognized container_type $container_type";
    exit 1;
  fi
    
  cmd_suppa="$base_cmd suppa.py $docker"

  GTF=$dGTF
  SALMON_RESULT=$dSALMON_RESULT
  SUPPA_REF=$dSUPPA_REF
else
  cmd_suppa="suppa.py"
fi

echo running isoforms
cmd_full="$cmd_suppa psiPerIsoform -g ${GTF} -e $SALMON_RESULT/tpm.txt -o $SALMON_RESULT/suppa2"
echo "$cmd_full"
$cmd_full

for ioe in $SUPPA_REF_local/*ioe; do
  name=$(basename "$ioe")
  name=${name/$(basename "$GTF")"._"/""}
  name=${name/suppa2./""}
  name=${name/_strict.ioe/""}
  #name=$(dirname $SALMON_RESULT)
  echo running "$name"
  #local ioe needs to be translated to docker path
  if [ -n "$docker" ]; then
    ioe=${SUPPA_REF}/$(basename "${ioe}")
  fi

  cmd_full="$cmd_suppa psiPerEvent -i ${ioe} -e $SALMON_RESULT/tpm.txt -o $SALMON_RESULT/suppa2.alternate_events.${name}"
  echo "$cmd_full"
  $cmd_full
done

