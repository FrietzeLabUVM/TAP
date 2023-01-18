#!/bin/bash
#$ -cwd

fq=$1
if [ -z $fq ]; then
  echo arg1 must be fastq.gz file.
  exit 1
fi
if [ ! -f $fq ]; then
  echo cannot find $fq. stop.
  exit 1
fi
newLength=$2
mkdir -p trimmed_fastq

trim_fq=trimmed_fastq/$(basename $fq .fastq.gz)".${newLength}.fastq"

if [ -z $newLength ]; then echo "2nd arg must me target read length"; exit 1; fi
if [ $fq == $trim_fq ]; then echo "could not derive new filename for $fq"; exit 1; fi

echo f is $fq
echo trim_fq is $trim_fq
echo newLength is $newLength
gunzip -c $fq | awk -v left_file=$trim_fq -v newLength=$newLength '{
        sub("length=[0-9]+", "length="newLength); 
        if(NR % 4 == 2 || NR % 4 == 0){
                left_read=substr($0, 1, newLength);
                print left_read > left_file; 
        }else{
                print $0 > left_file}
        }'

gzip $trim_fq
