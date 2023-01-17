rm -r /slipstream_old/home/joeboyd/lab_shared/scripts/vacc_rnaseq_pipeline/testing/test_alignment.pool2
bash run_STAR.noSort.sh \
  -f1 "/slipstream_old/home/joeboyd/lab_shared/scripts/vacc_rnaseq_pipeline/testing/test_data/test_dm6_1_R1_001.fastq.gz&/slipstream_old/home/joeboyd/lab_shared/scripts/vacc_rnaseq_pipeline/testing/test_data/test_dm6_2_R1_001.fastq.gz" \
  -wd /slipstream_old/home/joeboyd/lab_shared/scripts/vacc_rnaseq_pipeline/testing/test_alignment.pool2 \
  -idx /slipstream/home/joeboyd/indexes/DM6/STAR_INDEX \
  -o pool_1and2.withDocker \
  -f1s _R1_001.fastq.gz \
  -f2s _R2_001.fastq.gz \
  --docker tap 

ls -lha /slipstream_old/home/joeboyd/lab_shared/scripts/vacc_rnaseq_pipeline/testing/test_alignment.pool2/pool_1and2.withDocker*bam

bash run_STAR.noSort.sh \
  -f1 "/slipstream_old/home/joeboyd/lab_shared/scripts/vacc_rnaseq_pipeline/testing/test_data/test_dm6_1_R1_001.fastq.gz&/slipstream_old/home/joeboyd/lab_shared/scripts/vacc_rnaseq_pipeline/testing/test_data/test_dm6_2_R1_001.fastq.gz" \
  -wd /slipstream_old/home/joeboyd/lab_shared/scripts/vacc_rnaseq_pipeline/testing/test_alignment.pool2 \
  -idx /slipstream/home/joeboyd/indexes/DM6/STAR_INDEX \
  -o pool_1and2.noDocker \
  -f1s _R1_001.fastq.gz \
  -f2s _R2_001.fastq.gz

ls -lha /slipstream_old/home/joeboyd/lab_shared/scripts/vacc_rnaseq_pipeline/testing/test_alignment.pool2/pool_1and2.noDocker*bam