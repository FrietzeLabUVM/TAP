#!/bin/bash
#SLURM pipeline for RNAseq

mode=PE
sub_mode=sbatch
# umask 077 # rw permission for user only
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f1|--fastq1) F1="$2"; shift ;;
        -f2|--fastq2) F2="$2"; shift ;;
        -f1s|--f1_suffix) F1_suff="$2"; shift ;;
        -f2s|--f2_suffix) F2_suff="$2"; shift ;;
        -i|--inDir) in_path="$2"; shift ;;
        -p|--outPrefix) root="$2"; echo root is $root; shift ;;
        -o|--outDir) align_path="$2"; shift ;;
        -ref|--reference) ref=$2; shift ;;
        -idx|--starIndex) star_index="$2"; shift ;;
        -s|--suppaRef) suppa_ref="$2"; shift ;; 
        -g|--gtf) gtf="$2"; shift ;;
        -fa|--fasta) fasta="$2"; shift ;;
        -SE|--SE) mode=SE ;;
        -noSub|--noSub) sub_mode=bash ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

#if [ -z $ref ]; then echo need star index location for -idx; exit 1; fi
if [ -z $F1 ]; then echo need fastq1 as -f1! quit; exit 1; fi
if [ ! -z $in_path ]; then F1=${in_path}/${F1}; fi
F1=${F1//"&"/" "}
for f in $F1; do if [ ! -f $f ]; then echo fastq1 $f could not be found! quit; exit 1; fi; done
#if [ -z $gtf ]; then echo need gtf as -g; exit 1; fi
#if [ -z $fasta]; then echo need reference fasta as -fa; exit 1; fi
#if [ -z $align_path ]; then echo need alignment output path as -o; exit 1; fi

SCRIPTS=/users/j/r/jrboyd/dbgap_scripts/vacc_scripts

#relevant suffixes
if [ -z $F1_suff ]; then F1_suff="_R1_001.fastq.gz"; fi
if [ -z $F2_suff ]; then F2_suff="_R2_001.fastq.gz"; fi
suf_gz1="$F1_suff"
suf_gz2="$F2_suff"

if [ $mode = PE ]; then
 F2=${F1//$suf_gz1/$suf_gz2}
 #check that all F2 exist
 for f in $F2; do if [ ! -f $f ]; then echo fastq2 $f could not be found! quit; exit 1; fi; done
 F1a=($F1)
 F2a=($F2)
 #check same number of F1 and F2
 if [ ${#F1a[@]} != ${#F2a[@]} ]; then
   echo Differing numbers of R1 and R2 fastqs supplied! quit
   echo $F1
   echo $F2
   exit 1
 fi
 #check that F1 and F2 are different
 i=0
 while [ $i -lt ${#F1a[@]} ]; do
   if [ ${F1a[$i]} = ${F2a[$i]} ]; then echo problem with fastq pairing. quit; exit 1; fi
    i=$(( $i + 1 ))
 done
fi

suf_out_bam=".Aligned.out.bam"
suf_tx_bam=".Aligned.toTranscriptome.out.bam"
suf_sort_bam=".Aligned.sortedByCoord.out.bam"
suf_sort_bai=".Aligned.sortedByCoord.out.bam.bai"
suf_salmon_quant=".salmon_quant"
suf_log=".Log.out"
suf_align_stats=".Log.final.out"
suf_featr_count=".Aligned.sortedByCoord.out.featureCounts.txt"

#star_index=$ref/STAR_INDEX
#suppa_ref=$ref/SUPPA2
if [ -z $star_index ]; then star_index=$ref/STAR_INDEX; echo guessing star index as $star_index; fi
if [ ! -d $star_index ]; then echo star_index $star_index not found!; exit 1; fi
if [ -z $suppa_ref ]; then suppa_ref=$ref/SUPPA2; echo guessing suppa_ref as $suppa_ref; fi
if [ ! -d $suppa_ref ]; then echo suppa_ref $suppa_ref not found!; exit 1; fi
if [ -z $gtf ]; then gtf=$(readlink -m -f $ref/GTF/current.gtf); echo guessing gtf as $gtf; fi
if [ ! -f $gtf ]; then echo gtf $gtf not found! exit; exit 1; fi
if [ -z $fasta ]; then fasta=$(readlink -m -f $ref/FASTA/genome.fa); echo guessing fasta as $fasta; fi
if [ ! -f $fasta ]; then echo fasta $fasta not found! exit; exit 1; fi
if [ -z $tx ]; then tx=$(readlink -m -f $ref/FASTA/transcriptome.fa); echo guessing transcriptome fasta as $tx; fi
if [ ! -f $tx ]; then echo transcriptome fasta $tx not found! exit; exit 1; fi

for d in $star_index $suppa_ref; do
  if [ ! -d $d ]; then
    echo reference location $d must exists! exit
    exit 1
  fi
done

#output locations
if [ -z $align_path ]; then align_path=~/scratch/alignment_RNA-seq; echo using default alignment output path of $align_path; fi
if [ -z $root ]; then root=$(basename $(echo $F1 | awk -v FS=" +" '{print $1}') $suf_gz1); echo guessing root prefix is $root. provide with -p if incorrect;  fi
log_path=${align_path}/${root}.logs

#all must exist
mkdir -p $align_path
mkdir -p $log_path

#important files
out_bam=${align_path}/${root}${suf_out_bam}
tx_bam=${align_path}/${root}${suf_tx_bam}
sort_bam=${align_path}/${root}${suf_sort_bam}
salmon_out=${align_path}/${root}${suf_salmon_quant}
featr_out=${align_path}/${root}${suf_featr_count}

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
  qsub_cmd="sbatch -o $log_path/%x.%j.out -e $log_path/%x.%j.error"
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

$qsub_cmd echo_submission.sh $0 $#

#align script
F1=${F1//" "/"&"}
se_mode=""
if [ $mode = SE ]; then se_mode="-SE"; fi
align_qsub=$($qsub_cmd $SCRIPTS/run_STAR.noSort.sh -f1 $F1 -wd $align_path -idx $star_index -o $root -f1s $F1_suff -f2s $F2_suff $se_mode)
echo align_qsub $align_qsub
align_jid=$(parse_jid "$align_qsub")
echo align_jid $align_jid

#tx quant
salmon_sub_args="-d afterok:$align_jid -J salmon_quant"
if [ $sub_mode = "bash" ]; then salmon_sub_args=""; fi
salmon_jid=$(parse_jid "$($qsub_cmd $salmon_sub_args $SCRIPTS/run_salmon_quant.sh $tx_bam $tx)")
echo salmon_jid $salmon_jid
suppa2_sub_args="-d afterok:$salmon_jid -J suppa2"
if [ $sub_mode = "bash" ]; then suppa2_sub_args=""; fi
suppa2_jid=$(parse_jid "$($qsub_cmd $suppa2_sub_args $SCRIPTS/run_suppa2.sh $salmon_out $gtf $suppa_ref)")
echo suppa2_jid $suppa2_jid
#sort and index
index_sub_args="-d afterok:$align_jid -J bsortindex"
if [ $sub_mode = "bash" ]; then index_sub_args=""; fi
index_jid=$(parse_jid "$($qsub_cmd $index_sub_args $SCRIPTS/run_bam_sort_index.sh $out_bam)")
echo index_jid $index_jid

#counting
featureCounts_sub_args="-d afterok:$index_jid -J featureCounts"
if [ $sub_mode = "bash" ]; then featureCounts_sub_args=""; fi
featr_jid=$(parse_jid "$($qsub_cmd $featureCounts_sub_args $SCRIPTS/run_featureCounts.sh $sort_bam $gtf)")
exactSNP_sub_args="-d afterok:$index_jid -J exactSNP"
if [ $sub_mode = "bash" ]; then exactSNP_sub_args=""; fi
exact_jid=$(parse_jid "$($qsub_cmd $exactSNP_sub_args $SCRIPTS/run_exactSNP.all.sh $sort_bam $fasta)")

#subsetting bams to a region of interest is useful if only specific genes are relevant for tracks etc.
#subst_jid=$(parse_jid "$($qsub_cmd -d afterok:$index_jid -J subset_bam $SCRIPTS/run_ikaros_subset_bam.sh $sort_bam)")

echo featureCount_jid $featr_jid
echo exactSNP_jid $exact_jid


#cleanup #currently omitted, run_cmd.sh not updated for SLURM and only critical for bulk processing
#$qsub_cmd -hold_jid $salmon_jid -J cleanup_tx_bam $SCRIPTS/run_cmd.sh "if [ -d $salmon_out ]; then rm $tx_bam; fi"
#$qsub_cmd -hold_jid $featr_jid,$exact_jid,$subst_jid -J cleanup_bams $SCRIPTS/run_cmd.sh "if [ -s $featr_out ]; then rm $sort_bam; rm ${sort_bam}.bai; fi;"
#$qsub_cmd -hold_jid  $salmon_jid,$featr_jid,$exact_jid,$subst_jid -J finalize $SCRIPTS/run_cmd.sh "date > ${align_path}/${root}.complete"

completion_sub_args="-d afterok:$salmon_jid:$suppa2_jid:$featr_jid:$exact_jid -J completion"
if [ $sub_mode = "bash" ]; then completion_sub_args=""; fi
$qsub_cmd $completion_sub_args $SCRIPTS/write_completion.sh ${align_path}/${root}

