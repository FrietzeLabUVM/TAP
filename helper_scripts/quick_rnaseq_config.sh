#!/bin/bash
wd=/netfiles02/lockwood_lab/Kylie_Finnegan_transcriptomics/data/raw/seq
idx_dir=/gpfs2/pi-sfrietze/indexes/DM6
config_f=~/tmp_data/KF_fly_config.csv
#location to look for previous inputs
#reuse_dir="/users/c/g/cgao1/BRCA_progression/MCF10_Core/Output_dir_2023Jan"
reuse_dir=""

f1_suff=_1.fq.gz
f2_suff=_2.fq.gz


if [ ! -d $idx_dir ]; then
  echo "index directory not found: ${idx_dir}" 
fi
#data_path="/gpfs2/pi-sfrietze/data"
data_path=~/tmp_data
project_dir=$data_path/$(basename $wd)
mkdir -p "$project_dir"
fastq_dir="$project_dir"/"fastqs"
mkdir -p "$fastq_dir"
output_dir="$project_dir"/"pipeline_output"
mkdir -p "$output_dir"

#to systematically replace strings, create comma delimmited list of target,replacement
str_fixes="
AT1,MCF10AT1
CA1a,MCF10CA1a
DCIS,MCF10DCIS
r1,rep1
r2,rep2
"
str_fixes=""

str_targets=($(echo "$str_fixes" | awk -v RS="\n" -v FS="," '{print $1}'))
str_replacements=($(echo "$str_fixes" | awk -v RS="\n" -v FS="," '{print $2}'))

#1 based indexes
rep_pos=3
mark_pos=2
max_pos=3

fq_files=$(ls $wd/*$f1_suff)
echo "# config generated by: $(readlink -f $0)" > $config_f
echo "# $(date)" >> $config_f
echo "#CFG -i ${fastq_dir} -o ${output_dir} -ref $idx_dir" >> $config_f
for f in $fq_files; do
  #chmod a-w "$f"
  fq_file=$(basename "$f")
  if [ ! -f "$fastq_dir"/"$fq_file" ]; then
    ln -s "$f" "$fastq_dir"/"$fq_file"
  fi
  rep_name=${fq_file/$f1_suff/""}
  i=0
  while [ $i -lt ${#str_targets[@]} ]; do
    rep_name=${rep_name/${str_targets[$i]}/${str_replacements[$i]}}
    i=$(( i + 1 ))
  done
  rep_name=$(echo "$rep_name" | cut -d "_" -f 1-$max_pos)

  echo "$fq_file,$rep_name" >> $config_f

  #pool_name=$(echo "$rep_name" | awk -v FS="_" -v OFS="_" -v pos=$rep_pos '{$pos="pooled"; print $0}')
  #input_name=$(echo "$pool_name" | awk -v FS="_" -v OFS="_" -v pos=$mark_pos '{$pos="input"; print $0}')
  #echo "$fq_file,$rep_name,$pool_name,$input_name" >> $config_f
done

input_pre=$(cat $config_f | grep -v "#" | cut -d "," -f 4 | sort | uniq)
if [ ! -z $reuse_dir ]; then
for ip in $input_pre; do
  
  ip2=$(echo $ip | cut -d "_" -f 1-$(( max_pos - 1 )))
  #echo "$ip2"
  #link existing outputs
  for f in "$reuse_dir"/"$ip2"*; do
    outf=$output_dir/$(basename "$f")
    if [ "${outf/".logs"/""}" = "$outf" ]; then
        if [ ! -e "$outf" ]; then
            ln -s "$f" "$outf"
        fi
    fi
  done
  
  #create dummy fastq so no error from pipeline and append to config
  for f in "$reuse_dir"/"$ip2"*.Aligned.sortedByCoord.out.bam; do
    input_prefix=$(echo $(basename ${f/.Aligned.sortedByCoord.out.bam/""}) | grep -v "pooled")
    if [ -z $input_prefix ]; then continue; fi
    input_fastq=$fastq_dir/$input_prefix"_R1_001.fastq.gz"
    if [ ! -e "$input_fastq" ]; then
      echo "this is a dummy fastq" > "$input_fastq"
    fi

    fq_file=$(basename "$input_fastq")
    #append to config
    rep_name=$(echo "$input_prefix" | cut -d "_" -f 1-$max_pos)
    pool_name=$(echo "$rep_name" | awk -v FS="_" -v OFS="_" -v pos=$rep_pos '{$pos="pooled"; print $0}')
    input_name=$(echo "$pool_name" | awk -v FS="_" -v OFS="_" -v pos=$mark_pos '{$pos="input"; print $0}')
    echo "$fq_file,$rep_name,$pool_name,$input_name" >> $config_f
  done 
done
fi

echo "validate and use config file: $config_f"
#cat $config_f
