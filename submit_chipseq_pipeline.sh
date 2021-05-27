#!/bin/bash
set -e
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
no_model=

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
        -SE|--SE) read_mode=SE ;;
        -noSub|--noSub) sub_mode=bash ;;
        -p|--pipeline) pipeline="$2"; shift ;;
        -sl|--scriptLocation) scripts="$2"; shift ;;
        -noModel|--noModel) no_model="--noModel" ;;
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
        -SE|--SE) read_mode=SE ;;
        -noSub|--noSub) sub_mode=bash ;;
        -p|--pipeline) pipeline="$2"; shift ;;
        -sl|--scriptLocation) scripts="$2"; shift ;;
        -noModel|--noModel) no_model="--noModel" ;;
        *) echo "Unknown parameter passed: $1"; cat $SCRIPT_PATH/help_msg.txt; exit 1 ;;
    esac
    shift
done

#apply fallback defaults
if [ -z $F1_suff ]; then F1_suff=_R1_001.fastq.gz; fi
if [ -z $F2_suff ]; then F2_suff=_R2_001.fastq.gz; fi
if [ -z $pipeline ]; then pipeline=${SCRIPT_PATH}/chipseq_pipeline.sh; fi
if [ -z $pipeline2 ]; then pipeline2=${SCRIPT_PATH}/chipseq_pipeline.pooled.sh; fi
if [ -z $scripts ]; then scripts=${SCRIPT_PATH}; fi
if [ -z $align_path ]; then align_path=$(pwd); fi

echo pipeline is $pipeline

#check validity, must have input and (ref or all of idx,s,g,fa)
if [ -z $input ]; then echo input directory to find fastq in was never set, using current directory. use -i \(--inDir\) to specify.; input=$(pwd); fi;
if [ ! -d $input ]; then echo cannot find input directory ${input}. quit!; exit 1; fi;


#build final command
cmd=""
if [ ! -z $F1_suff ]; then cmd="$cmd --f1_suffix $F1_suff"; fi
if [ ! -z $F2_suff ] && [ $read_mode != SE ] ; then cmd="$cmd --f2_suffix $F2_suff"; fi
if [ ! -z $align_path ]; then cmd="$cmd --outDir $align_path"; fi
if [ ! -z $ref ]; then cmd="$cmd --reference $ref"; fi
if [ ! -z $star_index ]; then cmd="$cmd --starIndex $star_index"; fi
if [ ! -z $gtf ]; then cmd="$cmd --gtf $gtf"; fi
if [ ! -z $fasta ]; then cmd="$cmd --fasta $fasta"; fi
if [ ! -z $rDNA_index ]; then cmd="$cmd --rDNA_starIndex $rDNA_index"; fi
if [ ! -z $scripts ]; then cmd="$cmd --scriptLocation $scripts"; fi
if [ ! -z $no_model ]; then cmd="$cmd --noModel"; fi

if [ -z $star_index ]; then star_index=$ref/STAR_INDEX; echo guessing star index as $star_index; fi
if [ ! -d $star_index ]; then echo star_index $star_index not found!; exit 1; fi
CHRSIZES=$star_index/chrNameLength.txt

if [ $read_mode = SE ]; then cmd="$cmd -SE"; fi
if [ $sub_mode = bash ]; then cmd="$cmd -noSub"; fi 
#trim off leading space
cmd=${cmd/" "/""}

if [ ! -z $cfg ]; then 
  todo=$(cat $cfg | awk '/^[^#]/ {print $0 }')
else
  todo=$input/*$F1_suff
fi


#need arrays mapping pool_name to rep_bams, rep_jids, and input
declare -Ag input_pool2rep_bams
declare -Ag input_pool2rep_jids
declare -Ag input_pool2pool_jids
declare -Ag chip_rep2input_bams
declare -Ag chip_rep2input_jids
declare -Ag chip_pool2rep_bams
declare -Ag chip_pool2rep_jids
declare -Ag chip_pool2loose_jids
declare -Ag chip_pool2macs2_jids
declare -Ag chip_pool2input_bams
declare -Ag chip_pool2input_jids
declare -Ag chip_pool2pool_jids

#need some functions
parse_jid () { #parses the job id from output of qsub
        #echo $1
        if [ -z "$1" ]; then
          echo parse_jid expects output of qsub as first input but input was empty! stop
          exit 1
        fi
        JOBID=$(echo $1 | egrep -o -e "\b[0-9]+$")
        echo $JOBID;
}
parse_jid_by_name () { #parses the job id from output of qsub
        #echo $1
        if [ -z "$1" ]; then
          echo parse_jid_by_name expects output of pipeline script as first input but input was empty! stop
          exit 1
        fi
        if [ -z "$2" ]; then
          echo parse_jid_by_name expects string to match! stop
          exit 1
        fi
        #JOBID_LINE=$(echo $1 | grep $2)
        #JOBID=$(echo $JOBID_LINE | awk -v FS=" +" '{print $2}')
        line="${1}"
	tm="${2}"
	JOBID=$(echo ${line} | awk -v tomatch=${tm} -v RS="[ \n]+" -v found=0 '{if(found == 1){print $0; found=0} ; if($0 == tomatch){found=1};}')
        echo $JOBID;
}


suf_sort_bam=".Aligned.sortedByCoord.out.bam"
all_rep_names=""

#initial run through to validate format
for f_line in $todo; do
  echo $f_line
  if [ ! -f $f_line ]; then
    f1=$(echo $f_line | awk -v FS="," '{print $1}');
    rep_name=$(echo $f_line | awk -v FS="," '{print $2}');
    pool_name=$(echo $f_line | awk -v FS="," '{print $3}');
    input_name=$(echo $f_line | awk -v FS="," '{print $4}');
  else
    echo only config submission supported for ChIP right now! quit.
    exit 1
  fi
  if [ -z $f1 ] | [ -z $rep_name ] | [ -z $pool_name ] | [ -z $input_name ]; then
    echo check input format, need 4 csv fields: file,rep_name,pool_name,input_name
    echo file       : $f1
    echo rep_name   : $rep_name
    echo pool_name  : $pool_name
    echo input_name : $input_name
    exit 1
  fi
  all_rep_names="$all_rep_names $rep_name"
done
#check for duplicate rep_names
n_dupes=$(echo $all_rep_names | awk -v FS=" " -v OFS="\n" '{$1=$1; print $0}' | sort | uniq -d | wc -l)
if [ ! $n_dupes -eq 0 ]; then echo not all rep_name were unique! quit; exit 1; fi

echo RUN INPUT REPS
#run input reps and gather jids for pool inputs
for f_line in $todo; do
  f1=$(echo $f_line | awk -v FS="," '{print $1}');
  rep_name=$(echo $f_line | awk -v FS="," '{print $2}');
  pool_name=$(echo $f_line | awk -v FS="," '{print $3}');
  input_name=$(echo $f_line | awk -v FS="," '{print $4}');
  if [ $pool_name != $input_name ]; then #this is an input
    continue
  fi
  echo "input   : $rep_name -> $pool_name"

  ff1=""
  for f in $f1; do if [ -z "$ff1" ]; then ff1="$input/$(basename $f)"; else ff1="$ff1 $input/$(basename $f)"; fi; done
  f1=$ff1
  cmd_full="bash $pipeline -f1 ${f1//" "/&} --outPrefix $rep_name $cmd"
  input_rep_pipeout=$($cmd_full)
  align_jid=$(parse_jid_by_name "$input_rep_pipeout" index_jid)
  echo rep align_jid $align_jid
  echo rep bam ${rep_name}${suf_sort_bam}
  if [ -z ${input_pool2rep_bams[${pool_name}]} ]; then
    input_pool2rep_bams[$pool_name]=${align_path}/${rep_name}${suf_sort_bam}
    input_pool2rep_jids[$pool_name]=${align_jid}
  else
    input_pool2rep_bams[$pool_name]="${input_pool2rep_bams[${pool_name}]},${align_path}/${rep_name}${suf_sort_bam}"
    input_pool2rep_jids[$pool_name]="${input_pool2rep_jids[${pool_name}]}:${align_jid}"
  fi
done

#run pool input using rep input and gather jids reps/pool chip
echo POOL INPUT
for samp in "${!input_pool2rep_bams[@]}"; do
  jid="${input_pool2rep_jids["$samp"]}"
  if [ -z ${jid//:/} ]; then
    dep=""
  else
    dep="-d afterok:$jid"
  fi
  bam="${input_pool2rep_bams["$samp"]}"
  echo pooled input $samp
  echo $jid---$bam
  if [ $sub_mode != "bash" ]; then
    log_path=${align_path}/${samp}.logs
    mkdir -p $log_path
    pool_qsub=$(sbatch $dep -J pool_bams -o $log_path/%x.%j.out -e $log_path/%x.%j.error --export=PATH=$PATH ${scripts}/run_pool_bams.sh $bam ${align_path}/${samp}${suf_sort_bam})
    pool_jid=$(parse_jid "$pool_qsub")
    echo pool_jid $pool_jid
    input_pool2pool_jids[$samp]=${pool_jid}
  else
    bash ${scripts}/run_pool_bams.sh $bam ${align_path}/${samp}${suf_sort_bam}
  fi
done

echo RUN REP CHIP VS INPUT
#run rep chip using pool input and gather jids for pool chip
for f_line in $todo; do
  f1=$(echo $f_line | awk -v FS="," '{print $1}');
  rep_name=$(echo $f_line | awk -v FS="," '{print $2}');
  pool_name=$(echo $f_line | awk -v FS="," '{print $3}');
  input_name=$(echo $f_line | awk -v FS="," '{print $4}');
  if [ $pool_name = $input_name ]; then #this is an input
    continue
  fi
  echo "chip    : $rep_name -> $pool_name vs $input_name"

  ff1=""
  for f in $f1; do if [ -z "$ff1" ]; then ff1="$input/$(basename $f)"; else ff1="$ff1 $input/$(basename $f)"; fi; done
  f1=$ff1
  
  echo input_name  : ${input_name}
  echo input_jid   : ${input_pool2pool_jids[${input_name}]}

  cmd_full="bash $pipeline -f1 ${f1//" "/&} --outPrefix $rep_name -input_bam ${align_path}/${input_name}${suf_sort_bam} -input_jid ${input_pool2pool_jids[${input_name}]} $cmd"
  input_rep_pipeout=$($cmd_full)
  align_jid=$(parse_jid_by_name "$input_rep_pipeout" index_jid)
  echo rep align_jid $align_jid
  macs2_jid=$(parse_jid_by_name "$input_rep_pipeout" macs2_jid)
  echo rep macs2_jid $macs2_jid
  echo rep bam ${rep_name}${suf_sort_bam}
  if [ -z ${chip_pool2rep_bams[${pool_name}]} ]; then
    chip_pool2rep_bams[$pool_name]=${align_path}/${rep_name}${suf_sort_bam}
    chip_pool2rep_jids[$pool_name]=${align_jid}
    chip_pool2loose_jids[${pool_name}]=${macs2_jid}
  else
    chip_pool2rep_bams[$pool_name]="${chip_pool2rep_bams[${pool_name}]},${align_path}/${rep_name}${suf_sort_bam}"
    chip_pool2rep_jids[$pool_name]="${chip_pool2rep_jids[${pool_name}]}:${align_jid}"
    chip_pool2loose_jids[$pool_name]="${chip_pool2loose_jids[${pool_name}]}:${macs2_jid}"
  fi

  if [ -z ${chip_pool2pool_input_bam[${pool_name}]} ]; then
    chip_pool2pool_input_bam[$pool_name]=${align_path}/${input_name}${suf_sort_bam}
    chip_pool2pool_input_jid[$pool_name]=${input_pool2pool_jids[${input_name}]}
  else #check that new values match existing
    if [ ${chip_pool2pool_input_bam[$pool_name]} != ${align_path}/${input_name}${suf_sort_bam} ]; then
      echo found input mismatch for pooled bams! quit
      echo ${align_path}/${input_name}${suf_sort_bam} for $rep_name is not ${chip_pool2pool_input_bam[$pool_name]}
      exit 1
    fi
    if [ "${chip_pool2pool_input_jid[$pool_name]}" != "${input_pool2pool_jids[${input_name}]}" ]; then
      echo found input mismatch for pooled jids! quit
      echo ${input_pool2pool_jids[${input_name}]} for $rep_name is not ${chip_pool2pool_input_jid[$pool_name]}
      exit 1
    fi
  fi
done

echo POOL CHIP
#run pool chip using pool input and rep chip
for samp in "${!chip_pool2rep_bams[@]}"; do
  jid="${chip_pool2rep_jids["$samp"]}"
  if [ -z ${jid//:/} ]; then
    dep=""
  else
    dep="-d afterok:$jid"
  fi

  bam="${chip_pool2rep_bams["$samp"]}"
  echo pooled chip $samp
  echo $jid---$bam
  if [ $sub_mode != "bash" ]; then
    log_path=${align_path}/${samp}.logs
    mkdir -p $log_path
    pool_qsub=$(sbatch $dep -J pool_bams -o $log_path/%x.%j.out -e $log_path/%x.%j.error --export=PATH=$PATH ${scripts}/run_pool_bams.sh $bam ${align_path}/${samp}${suf_sort_bam})
    pool_jid=$(parse_jid "$pool_qsub")
    echo pool_jid $pool_jid
    chip_pool2pool_jids[$samp]=${pool_jid}
  else
    bash ${scripts}/run_pool_bams.sh $bam ${align_path}/${samp}${suf_sort_bam}
  fi
done

echo RUN POOL CHIP VS INPUT
for samp in "${!chip_pool2pool_jids[@]}"; do
  chip_jid="${chip_pool2pool_jids["$samp"]}"
  chip_bam=${align_path}/${samp}${suf_sort_bam}
  echo pooled chip $samp
  input_jid="${chip_pool2pool_input_jid["$samp"]}"
  input_bam="${chip_pool2pool_input_bam["$samp"]}"
  
  echo chip    : $chip_jid---$chip_bam
  echo vs
  echo input   : $input_jid---$input_bam
  cmd_pooled="bash $pipeline2 -chip_bam $chip_bam -chip_jid $chip_jid -input_bam $input_bam -input_jid $input_jid $cmd"
  echo pooled_cmd is $cmd_pooled

  chip_pool_pipeout=$($cmd_pooled)
  macs2_jid=$(parse_jid_by_name "$chip_pool_pipeout" macs2_jid)
  echo pool macs2_jid $macs2_jid

  chip_pool2macs2_jids[$samp]=${macs2_jid}
done

declare -Ag diff2rep_bam
declare -Ag diff2peaks
declare -Ag diff2pool_name
declare -Ag diff2jids

echo RUN REP COMPARISON
for samp in "${!chip_pool2rep_bams[@]}"; do
  jid="${chip_pool2loose_jids["$samp"]}"
  if [ -z ${jid//:/} ]; then
    dep=""
  else
    dep="-d afterok:$jid"
  fi

  bam="${chip_pool2rep_bams["$samp"]}"
  loose=${bam//.bam/_macs2_loose_peaks.narrowPeak}
  echo pooled chip $samp
  echo jid $jid
  echo bam $bam
  echo loose $loose
  if [ $sub_mode != "bash" ]; then
    log_path=${align_path}/${samp}.logs
    mkdir -p $log_path
    sbatch $dep -J IDR -o $log_path/%x.%j.out -e $log_path/%x.%j.error --export=PATH=$PATH ${scripts}/run_chip_rep_comparison.sh ${samp} ${loose} ${bam} ${CHRSIZES} ${align_path}
  else
    bash ${scripts}/run_chip_rep_comparison.sh ${samp} ${loose} ${bam} ${CHRSIZES} ${align_path}
  fi
done


echo RUN DIFFERENTIAL
for f_line in $todo; do
  f1=$(echo $f_line | awk -v FS="," '{print $1}');
  rep_name=$(echo $f_line | awk -v FS="," '{print $2}');
  pool_name=$(echo $f_line | awk -v FS="," '{print $3}');
  input_name=$(echo $f_line | awk -v FS="," '{print $4}');
  diff_group=$(echo $f_line | awk -v FS="," '{print $5}');
  if [ $pool_name = $input_name ]; then continue; fi #this is an input
  if [ -z $diff_group ]; then continue; fi
  if [ $diff_group = 0 ]; then continue; fi

  diff_jids="${chip_pool2loose_jids["$pool_name"]}":${chip_pool2macs2_jids["$pool_name"]}
  diff_jids=$(echo $diff_jids | sed 's/:\+/:/g')
  if [ -z ${diff_jids//:/} ]; then
    dep=""
  else
    dep="-d afterok:$diff_jids"
  fi

  rep_bam=${align_path}/${rep_name}${suf_sort_bam}
  pool_peak=${align_path}/${pool_name}_macs2_peaks.narrowPeak

  if [ -z ${diff2rep_bam[${diff_group}]} ]; then
    diff2rep_bam[$diff_group]=${rep_bam}
    diff2pool_name[$diff_group]=${pool_name}
    diff2jids[$diff_group]=${diff_jids}
    diff2peaks[$diff_group]=${pool_peak}
  else
    diff2rep_bam[$diff_group]="${diff2rep_bam[${diff_group}]},${rep_bam}"
    diff2pool_name[$diff_group]="${diff2pool_name[${diff_group}]},${pool_name}"
    diff2jids[$diff_group]="${diff2jids[${diff_group}]}:${diff_jids}"
    diff2peaks[$diff_group]="${diff2peaks[${diff_group}]},${pool_peak}"
  fi
done

if [ ! -z "${!diff2rep_bam[@]}" ]; then #check if any diff groups set
  for diff_group in ${!diff2rep_bam[@]}; do
    rep_bams=${diff2rep_bam[${diff_group}]}
    pool_names=${diff2pool_name[${diff_group}]}
    jids=${diff2jids[${diff_group}]}

    jids=$(echo $jids | sed 's/:\+/:/g')
    if [ -z ${jids//:/} ]; then
      dep=""
    else
      dep="-d afterok:$jids"
    fi


    peaks=${diff2peaks[${diff_group}]}
    echo pool_names ${pool_names}
    #pool_names=$(echo $pool_names | awk -v FS="," -v OFS="\n" '{$1=$1; print $0}' | sort | uniq | awk -v ORS="," '{$1=$1; print $0}') | sed 's/,$//'
    echo rep_bams $rep_bams
    #echo unique pool_names $pool_names
    log_path=${align_path}/diff_${diff_group}.logs
    mkdir -p $log_path
    sbatch $dep -J diff -o $log_path/%x.%j.out -e $log_path/%x.%j.error --export=PATH=$PATH ${scripts}/run_chip_differential.sh $pool_names $rep_bams $peaks
#    bash ${scripts}/run_chip_differential.sh
  done
fi
