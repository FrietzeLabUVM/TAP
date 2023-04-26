dock_test_dir=$(dirname "$0")
dock_test_dir=$(readlink -f $dock_test_dir)
scripts=$(dirname $dock_test_dir)
index_dir=~/lab_shared/indexes/DM6
#dock_test_dir=$scripts/testing_docker
#cd $dock_test_dir
test_dir=${scripts}/testing

rm -r ${test_dir}/test_alignment.pool2

for container_type in docker singularity native; do

if [ $container_type = docker ]; then
#docker environment
  test_cmd=docker
  out=pool_1and2.withDocker
  container_cmd="--docker jrboyd/tap"
elif [ $container_type = singularity ]; then
#singularity environment
  test_cmd=singularity
  out=pool_1and2.withSingularity
  singularity pull docker://jrboyd/tap 
  container_cmd="--singularity tap_latest.sif"
elif [ $container_type = native ]; then
#native environment
  test_cmd=STAR
  out=pool_1and2.withNative
  container_cmd=""
else
  echo container_type $container_type not recognized
fi

if ! command -v $test_cmd &> /dev/null
then
echo $container_type is not available
else
bash run_STAR.noSort.sh \
  -f1 "${test_dir}/test_data/test_dm6_1_R1_001.fastq.gz&${test_dir}/test_data/test_dm6_2_R1_001.fastq.gz" \
  -wd ${test_dir}/test_alignment.pool2 \
  -idx ${index_dir}/STAR_INDEX \
  -o $out \
  -f1s _R1_001.fastq.gz \
  -f2s _R2_001.fastq.gz \
  $container_cmd 

out_prefix=${test_dir}/test_alignment.pool2/${out}

bash run_salmon_quant.sh \
  ${out_prefix}.Aligned.toTranscriptome.out.bam \
  ${index_dir}/FASTA/transcriptome.fa \
  $container_cmd

bash run_suppa2.sh ${out_prefix}.salmon_quant ${index_dir}/GTF/current.gtf ${index_dir}/SUPPA2 $container_cmd

bash run_bam_sort_index.sh ${out_prefix}".Aligned.out.bam" $container_cmd

ls -lha ${out_prefix}*bam
fi
done

exit 0
#docker environment
if ! command -v docker &> /dev/null
then
echo docker is not available
else
#rm -r ${test_dir}/test_alignment.pool2
bash run_STAR.noSort.sh \
  -f1 "${test_dir}/test_data/test_dm6_1_R1_001.fastq.gz&${test_dir}/test_data/test_dm6_2_R1_001.fastq.gz" \
  -wd ${test_dir}/test_alignment.pool2 \
  -idx ~/lab_shared/indexes/DM6/STAR_INDEX \
  -o pool_1and2.withDocker \
  -f1s _R1_001.fastq.gz \
  -f2s _R2_001.fastq.gz \
  --docker jrboyd/tap 

ls -lha ${test_dir}/test_alignment.pool2/pool_1and2.withDocker*bam
fi

#singularity environment
if ! command -v singularity &> /dev/null
then
echo singularity is not available
else
#rm -r ${test_dir}/test_alignment.pool2
singularity pull docker://jrboyd/tap 
bash run_STAR.noSort.sh \
  -f1 "${test_dir}/test_data/test_dm6_1_R1_001.fastq.gz&${test_dir}/test_data/test_dm6_2_R1_001.fastq.gz" \
  -wd ${test_dir}/test_alignment.pool2 \
  -idx ~/lab_shared/indexes/DM6/STAR_INDEX \
  -o pool_1and2.withSingularity \
  -f1s _R1_001.fastq.gz \
  -f2s _R2_001.fastq.gz \
  --singularity tap_latest.sif

ls -lha ${test_dir}/test_alignment.pool2/pool_1and2.withSingularity*bam
fi

#native environment
if ! command -v STAR &> /dev/null
then
echo STAR is not available
else
bash run_STAR.noSort.sh \
  -f1 "${test_dir}/test_data/test_dm6_1_R1_001.fastq.gz&${test_dir}/test_data/test_dm6_2_R1_001.fastq.gz" \
  -wd ${test_dir}/test_alignment.pool2 \
  -idx ~/lab_shared/indexes/DM6/STAR_INDEX \
  -o pool_1and2.noDocker \
  -f1s _R1_001.fastq.gz \
  -f2s _R2_001.fastq.gz

ls -lha ${test_dir}/test_alignment.pool2/pool_1and2.noDocker*bam
fi

