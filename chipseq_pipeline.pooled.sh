#!/bin/bash
#SLURM pipeline for RNAseq

export PATH=/gpfs2/pi-sfrietze/bin:$PATH

SCRIPTS=$(dirname "$(readlink -f "$0")")

mode=SE
sub_mode=sbatch
# umask 077 # rw permission for user only
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -chip_bam|--chip_bam) sort_bam="$2"; shift ;;
        -chip_jid|--chip_jid) sort_jid="$2"; shift ;;
        -input_bam|--input_bam) input_bam="$2"; shift ;;
        -input_jid|--input_jid) input_jid="$2"; shift ;;
        -p|--outPrefix) root="$2"; echo root is $root; shift ;;
        -o|--outDir) align_path="$2"; shift ;;
        -ref|--reference) ref=$2; shift ;;
        -idx|--starIndex) star_index="$2"; shift ;;
        -s|--suppaRef) suppa_ref="$2"; shift ;; 
        -g|--gtf) gtf="$2"; shift ;;
        -fa|--fasta) fasta="$2"; shift ;;
        -rDNA|--rDNA_starIndex) rDNA_index="$2"; shift ;;
        -PE|--PE) mode=PE ;;
        -noSub|--noSub) sub_mode=bash ;;
        -sl|--scriptLocation) SCRIPTS="$2"; shift ;;
        -h|--help) cat $SCRIPTS/help_msg.txt; exit 0; shift ;;
        -f1s|--f1_suffix) shift ;; #these get dropped by pooled pipeline
        -f2s|--f2_suffix) shift ;;
        -i|--inDir) shift ;;
        *) echo "Unknown parameter passed: $1"; cat $SCRIPTS/help_msg.txt; exit 1 ;;
    esac
    shift
done

if [ ! -d $SCRIPTS ]; then echo could not find script directory $SCRIPTS, quit!; exit 1; fi
if [ -z $sort_bam ]; then echo need chip_bam as -chip_bam! quit; exit 1; fi
if [ -z $sort_jid ]; then echo need chip_jid as -chip_jid! quit; exit 1; fi
if [ -z $input_bam ]; then echo need input_bam as -input_bam! quit; exit 1; fi
if [ -z $input_jid ]; then echo need input_jid as -input_jid! quit; exit 1; fi

if [ -z $star_index ]; then star_index=$ref/STAR_INDEX; echo guessing star index as $star_index; fi
if [ ! -d $star_index ]; then echo star_index $star_index not found!; exit 1; fi
if [ -z $gtf ]; then gtf=$(readlink -m -f $ref/GTF/current.gtf); echo guessing gtf as $gtf; fi
if [ ! -f $gtf ]; then echo gtf $gtf not found! exit; exit 1; fi
if [ -z $fasta ]; then fasta=$(readlink -m -f $ref/FASTA/genome.fa); echo guessing fasta as $fasta; fi
if [ ! -f $fasta ]; then echo fasta $fasta not found! exit; exit 1; fi
if [ -z $rDNA_index ]; then
  rDNA_index=$(readlink -m -f ${ref}.rDNA/STAR_INDEX); echo guess rDNA star index as $rDNA_index;
else
  #rDNA star index must exist if specified instead of guessed
  if [ ! -d $rDNA_index ]; then echo Specified rDNA star index $rDNA_index was not found! quit; fi
fi


echo chip_bam is $sort_bam
echo chip_jid is $sort_jid
echo input_bam is $input_bam
echo input_jid is $input_jid
echo outPrefix is $root
echo outDir is $align_path
echo ref is $ref
echo star_index is $star_index

suf_sort_bam=".Aligned.sortedByCoord.out.bam"
suf_sort_bai=".Aligned.sortedByCoord.out.bam.bai"

if [ -z $gtf ]; then gtf=$(readlink -m -f $ref/GTF/current.gtf); echo guessing gtf as $gtf; fi
if [ ! -f $gtf ]; then echo gtf $gtf not found! exit; exit 1; fi

#output locations
if [ -z $root ]; then root=$(basename ${sort_bam/$suf_sort_bam/""}); echo guessing root prefix is $root. provide with -p if incorrect;  fi
log_path=${align_path}/${root}.logs

#all must exist
mkdir -p $log_path

parse_jid () { #parses the job id from output of qsub
        #echo $1
        if [ -z "$1" ]; then
        echo parse_jid expects output of qsub as first input but input was empty! stop
        exit 1
        fi
        #JOBID=$(echo $1 | awk '{print $3}') #returns JOBID
        JOBID=$(echo $1 | egrep -o -e "\b[0-9]+$")
        echo $JOBID;
}

if [ $sub_mode != "sbatch" ] && [ $sub_mode != "bash" ]; then echo sub_mode was $sub_mode, not one of sbatch or bash. quit!; exit 1; fi
if [ $sub_mode = "sbatch" ]; then 
  qsub_cmd="sbatch -o $log_path/%x.%j.out -e $log_path/%x.%j.error --export=PATH=$PATH"
elif [ $sub_mode = "bash" ]; then
  qsub_cmd="bash"
fi

if [ -f ${align_path}/${root}.complete ]; then 
  echo completion file found!
  echo delete ${align_path}/${root}.complete to rerun. quitting.
  exit 0
else
  echo no completion file found, starting run for ${align_path}/${root}
fi

date > ${align_path}/${root}.start

$qsub_cmd $SCRIPTS/echo_submission.sh $0 $#

#bigwigs
bw_sub_args="-d afterok:$sort_jid -J make_bigwigs"
if [ $sub_mode = "bash" ]; then bw_sub_args=""; fi
if [ $mode = SE ]; then
  bw_qsub=$($qsub_cmd $bw_sub_args $SCRIPTS/run_bam_to_bigwig.chip.sh -b $sort_bam -s $star_index/chrNameLength.txt -o ${sort_bam/.bam/""}.bigwigs)
else
  bw_qsub=$($qsub_cmd $bw_sub_args $SCRIPTS/run_bam_to_bigwig.chip.sh -b $sort_bam -s $star_index/chrNameLength.txt -o ${sort_bam/.bam/""}.bigwigs -pe)
fi
bw_jid=$(parse_jid "$bw_qsub")
echo bw_jid $bw_jid

#SNPs
exactSNP_sub_args="-d afterok:$sort_jid -J exactSNP"
if [ $sub_mode = "bash" ]; then exactSNP_sub_args=""; fi
exact_jid=$(parse_jid "$($qsub_cmd $exactSNP_sub_args $SCRIPTS/run_exactSNP.all.sh $sort_bam $fasta)")

echo exactSNP_jid $exact_jid

if [ ! -z $input_bam ]; then #treat as chip sample and call peaks
  macs2_cmd="$SCRIPTS/run_chip_vs_input.sh -t $sort_bam -i $input_bam  -o $align_path -p $(basename $sort_bam .bam)_macs2 -g $(basename $ref) -s $star_index/chrNameLength.txt"
  if [ -z "$sort_jid$input_jid" ]; then
    macs2_sub_args="-J macs2"
  elif [ -z $input_jid ]; then #no input job dependency
    macs2_sub_args="-d afterok:$sort_jid -J macs2"
  elif [ -z $sort_jid ]; then 
    macs2_sub_args="-d afterok:$input_jid -J macs2"
  else #has input job dependency
    macs2_sub_args="-d afterok:$sort_jid:$input_jid -J macs2"
  fi
  if [ $sub_mode = "bash" ]; then macs2_sub_args=""; fi
  macs2_qsub=$($qsub_cmd $macs2_sub_args $macs2_cmd)
  macs2_jid=$(parse_jid "$macs2_qsub")
  echo macs2_jid $macs2_jid

  #loose peak
  #broad peak
  #run_bdgcmp.sh
  #run_bdg2bw.sh
  #run_np2bb.sh
fi


if [ -z $macs2_jid ]; then
  finish_sub_args="-d afterany:$exact_jid:$bw_jid -J finish"
  completion_sub_args="-d afterok:$exact_jid:$bw_jid -J completion"
else
  finish_sub_args="-d afterany:$exact_jid:$bw_jid:$macs2_jid -J finish"
  completion_sub_args="-d afterok:$exact_jid:$bw_jid:$macs2_jid -J completion"
fi
if [ $sub_mode = "bash" ]; then completion_sub_args=""; finish_sub_args=""; fi
complete_qsub=$($qsub_cmd $completion_sub_args $SCRIPTS/write_completion.sh ${align_path}/${root})
$qsub_cmd $finish_sub_args $SCRIPTS/write_finish.sh ${align_path}/${root}

complete_jid=$(parse_jid "$complete_qsub")
echo complete_jid $complete_jid
