#!/bin/bash
wd=~/lab_shared/data/Frietze_SK_iLabs_15333_Cahan05_12092020
idx=~/dbgap_scripts/vacc_scripts/DM6_STAR_INDEX
out=~/pipelines/dm6_v2

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


for f in $wd/*/*_R1_001.fastq.gz; do
  echo $f;
  bam=$(basename $f)
  bam=${bam/$suf_gz1/""}${suf_out_bam}
  star_cmd="run_STAR.noSort.sh -f1 $f -idx $idx -wd $out"
  echo $star_cmd
  star_sb="sbatch $star_cmd"
  star_out=$($star_sb)
  id_star=$(parse_jid "$star_out")
  echo STAR id is $id_star
  sbatch -d afterok:$id_star run_bam_sort_index.sh ${out}/${bam}
  
done
