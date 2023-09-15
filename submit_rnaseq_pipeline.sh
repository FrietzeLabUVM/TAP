#!/bin/bash

SCRIPT_PATH=$(dirname "$(readlink -f "$0")")
#echo $SCRIPT_PATH

#parameters parsing
#1 get config param if specified
#2 if config, parse params
#3 override any config params with remaining params

in_param="$@"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help) cat $SCRIPT_PATH/help_msg.txt; exit 0; shift ;;
    esac
    shift
done

#check if config specified
set -- ${in_param}
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c|--config) cfg="$2"; if [ ! -f $cfg ]; then echo cannnot find config file $cfg. quit!; exit 1; fi; shift ;;
    esac
    shift
done

read_mode=PE
sub_mode=sbatch
docker=""
singularity=""

#parse args specified in config file by lines starting with #CFG
if [ ! -z $cfg ]; then
  echo gathering parameters from config file $cfg
  args=$(cat $cfg | awk -v cfg_prefix="#CFG" -v ORS=" " '{if ($1 == cfg_prefix){$1 = ""; print $0}}')
  args="${args//\~/$HOME}"
  if [ ! -z "$args" ]; then
    echo $args
    set -- $args
    while [[ "$#" -gt 0 ]]; do
      case $1 in
	-c|--config) echo ignoring config file specified in config file.; shift ;;
	-f1s|--f1_suffix) F1_suff="$2"; shift ;;
	-f2s|--f2_suffix) F2_suff="$2"; shift ;;
	-i|--inDir) input="$2"; shift ;;
	-o|--outDir) align_path="$2"; shift ;;
	-ref|--reference) ref="$2"; shift ;;
	-idx|--starIndex) star_index="$2"; shift ;;
	-s|--suppaRef) suppa_ref="$2"; shift ;;
	-g|--gtf) gtf="$2"; shift ;;
	-fa|--fasta) fasta="$2"; shift ;;
  -rDNA|--rDNA_starIndex) rDNA_index="$2"; shift ;;
  -PE|--PE) read_mode=PE ;;
  -SE|--SE) read_mode=SE ;;
  -noSub|--noSub) sub_mode=bash ;;
  -p|--pipeline) pipeline="$2"; shift ;;
  -sl|--scriptLocation) scripts="$2"; shift ;;
  -docker|--docker) docker="$2"; shift ;;
  -singularity|--singularity) singularity="$2"; shift ;;
	*) echo "Unknown parameter passed: $1"; cat $SCRIPT_PATH/help_msg.txt; exit 1 ;;
      esac
      shift
    done
  fi
fi

#parse remaining args from command
set -- ${in_param}
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c|--config) shift ;;
        -f1s|--f1_suffix) F1_suff="$2"; shift ;;
        -f2s|--f2_suffix) F2_suff="$2"; shift ;;
        -i|--inDir) input="$2"; shift ;;
        -o|--outDir) align_path="$2"; shift ;;
        -ref|--reference) ref="$2"; shift ;;
        -idx|--starIndex) star_index="$2"; shift ;;
        -s|--suppaRef) suppa_ref="$2"; shift ;;
        -g|--gtf) gtf="$2"; shift ;;
        -fa|--fasta) fasta="$2"; shift ;;
        -rDNA|--rDNA_starIndex) rDNA_index="$2"; shift ;;
	-PE|--PE) read_mode=PE ;;
        -SE|--SE) read_mode=SE ;;
        -noSub|--noSub) sub_mode=bash ;;
        -p|--pipeline) pipeline="$2"; shift ;;
        -sl|--scriptLocation) scripts="$2"; shift ;;
        -docker|--docker) docker="$2"; shift ;;
        -singularity|--singularity) singularity="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; cat $SCRIPT_PATH/help_msg.txt; exit 1 ;;
    esac
    shift
done

#apply fallback defaults
if [ -z $F1_suff ]; then F1_suff=_R1_001.fastq.gz; fi
if [ -z $F2_suff ]; then F2_suff=_R2_001.fastq.gz; fi
if [ -z $scripts ]; then scripts=${SCRIPT_PATH}; fi
if [ -z $pipeline ]; then pipeline=${scripts}/rnaseq_pipeline.sh; fi

echo pipeline is $pipeline

if [ -z $align_path ]; then echo output directory was never set, creating TAP_output in current directory.; align_path=$(pwd)/TAP_output; fi

#check validity, must have input and (ref or all of idx,s,g,fa)
if [ -z $input ]; then echo input directory to find fastq in was never set, using current directory. use -i \(--inDir\) to specify.; input=$(pwd); fi;
if [ ! -d $input ]; then echo cannot find input directory ${input}. quit!; exit 1; fi;

#only allow 1 container type
if [ -n "$docker" ] && [ -n "$singularity" ]; then
  echo Only 1 of docker or signularity should be set. Quit!
  exit 1
fi

#build final command
cmd=""
if [ ! -z $F1_suff ]; then cmd="$cmd --f1_suffix $F1_suff"; fi
if [ ! -z $F2_suff ] && [ $read_mode != SE ] ; then cmd="$cmd --f2_suffix $F2_suff"; fi
if [ ! -z $root ]; then cmd="$cmd --outPrefix $root"; fi
if [ ! -z $align_path ]; then cmd="$cmd --outDir $align_path"; fi
if [ ! -z $ref ]; then cmd="$cmd --reference $ref"; fi
if [ ! -z $star_index ]; then cmd="$cmd --starIndex $star_index"; fi
if [ ! -z $suppa_ref ]; then cmd="$cmd --suppaRef $suppa_ref"; fi
if [ ! -z $gtf ]; then cmd="$cmd --gtf $gtf"; fi
if [ ! -z $fasta ]; then cmd="$cmd --fasta $fasta"; fi
if [ ! -z $rDNA_index ]; then cmd="$cmd --rDNA_starIndex $rDNA_index"; fi
if [ ! -z $scripts ]; then cmd="$cmd --scriptLocation $scripts"; fi
if [ ! -z $docker ]; then cmd="$cmd --docker $docker"; fi
if [ ! -z $singularity ]; then cmd="$cmd --singularity $singularity"; fi

if [ $read_mode = SE ]; then cmd="$cmd -SE"; fi
if [ $sub_mode = bash ]; then cmd="$cmd -noSub"; fi 
#trim off leading space
cmd=${cmd/" "/""}

if [ ! -z $cfg ]; then 
  todo=$(cat $cfg | awk '/^[^#]/ {print $0 }')
else
  todo=$input/*$F1_suff
fi

for f1 in $todo; do
  if [ ! -f $f1 ]; then 
    root=$(echo $f1 | awk -v FS="," '{print $2}');
    f1=$(echo $f1 | awk -v FS="," '{print $1}');# | awk -v FS=" " -v OFS="&" '$1=$1; {print $0}');
  fi
  #echo $f1
  #f1=${f1//" "/"&"}
  #echo $f1
  f1=${f1//"&"/" "}
  #echo $f1
  #exit 0
  ff1=""
  for f in $f1; do if [ -z "$ff1" ]; then ff1="$input/$(basename $f)"; else ff1="$ff1 $input/$(basename $f)"; fi; done
  #if [ ! -f $f1 ]; then f1=$input/$f1; fi
  f1=$ff1
  for f in $f1; do if [ ! -f $f ]; then echo fastq1 was not found, $f. quit!; exit 1; fi; done
  if [ $read_mode != SE ]; then
    f2=${f1//$F1_suff/$F2_suff}
    for f in $f2; do if [ ! -f $f ]; then echo fastq2 was not found, $f. quit!; exit 1; fi; done
    cmd_full="bash $pipeline -f1 ${f1//" "/&} -f2 ${f2//" "/&} $cmd"
  else
    cmd_full="bash $pipeline -f1 ${f1//" "/&} $cmd"
  fi
  if [ ! -z $root ]; then cmd_full="$cmd_full --outPrefix $root"; fi
  echo $cmd_full
  $cmd_full
done
