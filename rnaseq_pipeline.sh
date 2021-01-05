#!/bin/bash
#SLURM pipeline for RNAseq

# umask 077 # rw permission for user only
srr=$1
root=$2

#output locations
data_path=/slipstream/home/dbgap/data
align_path=${data_path}/alignment_RNA-Seq
log_path=${align_path}/${root}.logs

#all must exist
mkdir -p $data_path
mkdir -p $align_path
mkdir -p $log_path

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
        #JOBID=$(awk -v RS=[0-9]+ '{print RT+0;exit}' <<< "$1") #returns JOBID  #original method was dependent on specific awk version apparently for RT variable
#        JOBID=$(echo $1 | awk 'BEGIN {RS=" "} {if ($1 ~ "[0-9]+"){print $1; exit}}') #returns JOBID
        JOBID=$(echo $1 | awk '{print $3}') #returns JOBID
#        echo $JOBID >> $ALL_JIDS
        echo $JOBID;
}
#export -f parse_jid

qsub_cmd="qsub -V -cwd -p -500 -o $log_path -e $log_path"

if [ -f ${align_path}/${root}.complete ]; then 
  echo completion file found!
  echo delete ${align_path}/${root}.complete to rerun. quitting.
  exit 0
else
  echo no completiong file, starting run for ${align_path}/${root}
fi
#exit 0
date > ${align_path}/${root}.start

#align script
align_qsub=$($qsub_cmd ~/scripts/run_STAR.noSort.sh $fq_gz1 $align_path)
echo align_qsub $align_qsub
align_jid=$(parse_jid "$align_qsub")
echo align_jid $align_jid

#tx quant
salmon_jid=$(parse_jid "$($qsub_cmd -hold_jid $align_jid -N salmon_quant ~/scripts/run_salmon_quant.sh $tx_bam)")
echo $salmon_jid
#$qsub_cmd -hold_jid $salmon_jid -N cleanup_tx_bam ~/scripts/run_cmd.sh "if [ -d $salmon_out ]; then rm $tx_bam; fi"
#qsub -hold_jid $salmon_jid -V -cwd -N cleanup_tx_bam ~/scripts/run_cmd.sh "echo rm $tx_bam;"
suppa2_jid=$(parse_jid "$($qsub_cmd -hold_jid $salmon_jid -N suppa2 ~/scripts/run_suppa2.sh $salmon_out)")
echo $suppa2_jid
#sort and index
index_jid=$(parse_jid "$($qsub_cmd -hold_jid $align_jid ~/scripts/run_bam_sort_index.sh $out_bam)")
echo $index_jid
#counting
featr_jid=$(parse_jid "$($qsub_cmd -hold_jid $index_jid ~/scripts/run_featureCounts.sh $sort_bam)")
exact_jid=$(parse_jid "$($qsub_cmd -hold_jid $index_jid ~/scripts/run_exactSNP.all.sh $sort_bam)")
subst_jid=$(parse_jid "$($qsub_cmd -hold_jid $index_jid ~/scripts/run_ikaros_subset_bam.sh $sort_bam)")

echo $featr_jid
echo $exact_jid


#cleanup #currently omitted, run_cmd.sh not updated for SLURM and only critical for bulk processing
#$qsub_cmd -hold_jid $salmon_jid -N cleanup_tx_bam ~/scripts/run_cmd.sh "if [ -d $salmon_out ]; then rm $tx_bam; fi"
#$qsub_cmd -hold_jid $featr_jid,$exact_jid,$subst_jid -N cleanup_bams ~/scripts/run_cmd.sh "if [ -s $featr_out ]; then rm $sort_bam; rm ${sort_bam}.bai; fi;"
#$qsub_cmd -hold_jid  $salmon_jid,$featr_jid,$exact_jid,$subst_jid -N finalize ~/scripts/run_cmd.sh "date > ${align_path}/${root}.complete"



