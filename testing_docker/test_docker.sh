dock_test_dir=$(dirname "$0")
dock_test_dir=$(readlink -f $dock_test_dir)
scripts=$(dirname $dock_test_dir)
#dock_test_dir=$scripts/testing_docker
#cd $dock_test_dir
test_dir=${scripts}/testing

#docker environment
if ! command -v docker &> /dev/null
then
echo docker is not available
else
rm -r ${test_dir}/test_alignment.pool2
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
rm -r ${test_dir}/test_alignment.pool2
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

