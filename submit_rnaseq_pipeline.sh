#!/bin/bash
#parameters parsing
#1 get config param if specified
#2 if config, parse params
#3 override any config params with remaining params

in_param="$@"

#check if config specified
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c|--config) cfg="$2"; if [ ! -f $cfg ]; then echo cannnot find config file $cfg. quit!; exit 1; fi; shift ;;
    esac
    shift
done

#parse args specified in config file by lines starting with #CFG
if [ ! -z $cfg ]; then
  args=$(cat $cfg | awk -v cfg_prefix="#CFG" -v ORS=" " '{if ($1 == cfg_prefix){$1 = ""; print $0}}')
  if [ ! -z $args ]; then
    echo $args
    set -- $args
    while [[ "$#" -gt 0 ]]; do
      case $1 in
	-c|--config) echo ignoring config file specified in config file.; shift ;;
	-f1s|--f1_suffix) F1_suff="$2"; shift ;;
	-f2s|--f2_suffix) F2_suff="$2"; shift ;;
	-i|--inDir) input="$2"; shift ;;
	-p|--outPrefix) root="$2"; shift ;;
	-o|--outDir) align_path="$2"; shift ;;
	-ref|--reference) ref="$2"; shift ;;
	-idx|--starIndex) star_index="$2"; shift ;;
	-s|--suppaRef) suppa_ref="$2"; shift ;;
	-g|--gtf) gtf="$2"; shift ;;
	-fa|--fasta) fasta="$2"; shift ;;
	*) echo "Unknown parameter passed: $1"; exit 1 ;;
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
        -i|--inDir) input="$2"; if [ ! -d $input ]; then echo cannot find input directory ${input}. quit!; exit 1; fi; ;;
        -p|--outPrefix) root="$2"; shift ;;
        -o|--outDir) align_path="$2"; shift ;;
        -ref|--reference) ref="$2"; shift ;;
        -idx|--starIndex) star_index="$2"; shift ;;
        -s|--suppaRef) suppa_ref="$2"; shift ;;
        -g|--gtf) gtf="$2"; shift ;;
        -fa|--fasta) fasta="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

#apply fallback defaults
if [ -z $F1_suff ]; then F1_suff=_R1_001.fastq.gz; fi
if [ -z $F2_suff ]; then F2_suff=_R2_001.fastq.gz; fi

#check validity, must have input and (ref or all of idx,s,g,fa)

if [ -z $input ]; then echo input directory to find fastq in was never set, using current directory. use -i \(--inDir\) to specify.; input=$(pwd); fi;
if [ ! -d $input ]; then echo cannot find input directory ${input}. quit!; exit 1; fi;


#build final command
if [ ! -z $F1_suff ]; then cmd="$cmd --f1_suffix $F1_suff"; fi
if [ ! -z $F2_suff ]; then cmd="$cmd --f2_suffix $F2_suff"; fi


for f1 in $input/*$F1_suff; do
  if [ ! -f $f1 ]; then echo fastq1 was not valid, $f1. quit!; exit 1; fi
  f2=${f1/$F1_suff/""}$F2_suff
  echo $f1 $f2
  #bash rnaseq_pipeline.sh -f1 $f1 -ref $ref
done
