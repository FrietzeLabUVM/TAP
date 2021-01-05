#!/bin/bash
#SLURM pipeline for RNAseq

# umask 077 # rw permission for user only
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f1|--fastq1) F1="$2"; shift ;;
        -f2|--fastq2) F2="$2"; shift ;;
        -p|--prefix) root="$2"; shift ;;
        -wd|--workdir) align_path="$2"; shift ;;
        -idx|--starindex) star_idx="$2"; shift ;;
        -g|--gtf) gtf="$2"; shift ;;
        -fa|--fasta) fasta="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ -z $star_idx ]; then echo need star index location for -idx; exit 1; fi
if [ -z $F1 ]; then echo need fastq1 as -f1; exit 1; fi
if [ -z $gtf ]; then echo need gtf as -g; exit 1; fi
if [ -z $fasta]; then echo need reference fasta as -fa; exit 1; fi
#if [ -z $align_path ]; then echo need alignment output path as -wd; exit 1; fi

SCRIPTS=/users/j/r/jrboyd/dbgap_scripts/vacc_scripts

#relevant suffixes
suf_gz1="_R1_001.fastq.gz"
suf_gz2="_R2_001.fastq.gz"
suf_out_bam=".Aligned.out.bam"
suf_tx_bam=".Aligned.toTranscriptome.out.bam"
suf_sort_bam=".Aligned.sortedByCoord.out.bam"
suf_sort_bai=".Aligned.sortedByCoord.out.bam.bai"
suf_salmon_quant=".salmon_quant"
suf_log=".Log.out"
suf_align_stats=".Log.final.out"
suf_featr_count=".Aligned.sortedByCoord.out.featureCounts.txt"

#output locations
if [ -z $align_path ]; then align_path=~/scratch/alignment_RNA-seq; echo using default alignment output path of $align_path; fi
if [ -z $root ]; then root=$(basename $F1 $suf_gz1); echo guessing root prefix is $root. provide with -p if incorrect;  fi
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

qsub_cmd="sbatch -o $log_path/%x.%j.out -e $log_path/%x.%j.error"

if [ -f ${align_path}/${root}.complete ]; then 
  echo completion file found!
  echo delete ${align_path}/${root}.complete to rerun. quitting.
  exit 0
else
  echo no completion file found, starting run for ${align_path}/${root}
fi

date > ${align_path}/${root}.start

#align script
align_qsub=$($qsub_cmd $SCRIPTS/run_STAR.noSort.sh -f1 $F1 -wd $align_path -idx $star_idx)
echo align_qsub $align_qsub
align_jid=$(parse_jid "$align_qsub")
echo align_jid $align_jid

#tx quant
salmon_jid=$(parse_jid "$($qsub_cmd -d afterok:$align_jid -J salmon_quant $SCRIPTS/run_salmon_quant.sh $tx_bam $gtf)")
echo $salmon_jid
suppa2_jid=$(parse_jid "$($qsub_cmd -d afterok:$salmon_jid -J suppa2 $SCRIPTS/run_suppa2.sh $salmon_out)")
echo $suppa2_jid
#sort and index
index_jid=$(parse_jid "$($qsub_cmd -d afterok:$align_jid -J bsortindex $SCRIPTS/run_bam_sort_index.sh $out_bam)")
echo $index_jid
#counting
featr_jid=$(parse_jid "$($qsub_cmd -d afterok:$index_jid -J featureCounts $SCRIPTS/run_featureCounts.sh $sort_bam $gtf)")
exact_jid=$(parse_jid "$($qsub_cmd -d afterok:$index_jid -J exactSNP $SCRIPTS/run_exactSNP.all.sh $sort_bam $fasta)")

#subsetting bams to a region of interest is useful if only specific genes are relevant for tracks etc.
#subst_jid=$(parse_jid "$($qsub_cmd -d afterok:$index_jid -J subset_bam $SCRIPTS/run_ikaros_subset_bam.sh $sort_bam)")

echo $featr_jid
echo $exact_jid


#cleanup #currently omitted, run_cmd.sh not updated for SLURM and only critical for bulk processing
#$qsub_cmd -hold_jid $salmon_jid -J cleanup_tx_bam $SCRIPTS/run_cmd.sh "if [ -d $salmon_out ]; then rm $tx_bam; fi"
#$qsub_cmd -hold_jid $featr_jid,$exact_jid,$subst_jid -J cleanup_bams $SCRIPTS/run_cmd.sh "if [ -s $featr_out ]; then rm $sort_bam; rm ${sort_bam}.bai; fi;"
#$qsub_cmd -hold_jid  $salmon_jid,$featr_jid,$exact_jid,$subst_jid -J finalize $SCRIPTS/run_cmd.sh "date > ${align_path}/${root}.complete"



