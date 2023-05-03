#!/bin/bash
#SBATCH --nodes=1                               # Request one core
#SBATCH --ntasks-per-node=1                               # Request one node (if you request more than one core with -n, also using
#SBATCH --cpus-per-task=8                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 1-06:00                         # Runtime in D-HH:MM format
#SBATCH -p bluemoon                           # Partition to run in
#SBATCH --mem=31000                        # Memory total in MB (for all cores)
#SBATCH -o fasterq_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e fasterq_%j.err                 # File to which STDERR will be written, including job ID

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -s|--ssr) srr="$2"; shift ;;
        -p|--prefix) root="$2"; shift ;;
        -o|--outdir) fetch_path="$2"; shift ;;
        -n|--ncores) ncores="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done


DUMP_PATH=/slipstream/home/joeboyd/bin/sratoolkit
if [ -z "$srr" ]; then echo "(-s|--ssr) parameter is required!"; exit 1; fi
if [ -z "$root" ]; then echo "(-p|--prefix) parameter is required!"; exit 1; fi

if [ -z "$fetch_path" ]; then
  fetch_path=$(pwd)
fi
if [ -z "$ncores" ]; then
  ncores=8
fi

echo srr is "$srr"
echo root is "$root"
echo fetch_path is "$fetch_path"

mkdir -p "$fetch_path"
cd "$fetch_path" || exit

dl_fq1=${fetch_path}/${srr}".sra_1.fastq"
dl_fq2=${fetch_path}/${srr}".sra_2.fastq"
dl_fqse=${fetch_path}/${srr}".sra.fastq"

fq1=${fetch_path}/${root}"_R1_001.fastq"
fq2=${fetch_path}/${root}"_R2_001.fastq"
fqse=${fetch_path}/${root}"_R1_001.fastq"

if [ -f "${fq1}".gz ] && [ -f "${fq2}".gz ]; then
 echo final fastqs are present, will not rerun.;
 exit 0;
elif [ -f "${fqse}".gz ]; then
 echo final fastq is present, will not rerun.;
 exit 0;
fi

if [ -d "${fetch_path}"/"${srr}" ]; then
  echo found prefetch for "$srr", will not rerun.;
else
  $DUMP_PATH/prefetch --max-size 500G -O "$fetch_path"  "$srr"
fi

if [ -f "${dl_fq1}" ] && [ -f "${dl_fq2}" ]; then
  echo found dumped fastqs, will not rerun.
else
  if [ -f "${dl_fqse}" ]; then
    echo found dumped fastq, will not rerun.
  else
    $DUMP_PATH/fasterq-dump -O "$fetch_path" --split-3 "${fetch_path}"/"${srr}".sra --threads $ncores
  fi
fi

if [ -f "${dl_fq1}" ]; then
  mv "${dl_fq1}" "${fq1}"
  pigz -p $ncores "${fq1}"
fi

if [ -f "${dl_fq2}" ]; then
  mv "${dl_fq2}" "${fq2}"
  pigz -p $ncores "${fq2}"
fi

if [ -f "${dl_fqse}" ]; then
  mv "${dl_fqse}" "${fqse}"
  pigz -p $ncores "${fqse}"
fi

