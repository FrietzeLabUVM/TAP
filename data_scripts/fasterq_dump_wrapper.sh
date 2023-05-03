#!/bin/bash
#SBATCH --nodes=1                               # Request one core
#SBATCH --ntasks-per-node=1                               # Request one node (if you request more than one core with -n, also using
#SBATCH --cpus-per-task=8                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 1-06:00                         # Runtime in D-HH:MM format
#SBATCH -p bluemoon                           # Partition to run in
#SBATCH --mem=31000                        # Memory total in MB (for all cores)
#SBATCH -o fasterq_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e fasterq_%j.err                 # File to which STDERR will be written, including job ID

usage="$(basename "$0") [-h] -s SSR1234 -p fastq_file_prefix -o output -n 8\n
\n
Runs fasterq-dump and pigz (parellel gzip) with 8 cores by default.\n
Final fastq.gz file will have specified prefix.\n
\n
where:\n
    -s|--ssr     An SSR accession for sra  \n
    -p|--prefix  Prefix for fastq output. Will have _R1_001.fastq.gz appended, and _R2_001.fastq.gz if data is PE.\n
    -o|--outdir  Output directory, will default to current directory if not specified.\n
    -n|--ncores  The number of cores to use for fasterq-dump and pigz, defaults to 8. Remember to add --cpus-per-task=? if using sbatch.\n
    -h|--help    Print this help and exit.\n
 "

if [ -z $1 ]; then echo -e $usage; exit 1; fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -s|--ssr) srr="$2"; shift ;;
        -p|--prefix) root="$2"; shift ;;
        -o|--outdir) fetch_path="$2"; shift ;;
        -n|--ncores) ncores="$2"; shift ;;
        -h|--help) echo -e $usage; exit 0 ;;
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
echo prefix is "$root"
echo fetch_path is "$fetch_path"

mkdir -p "$fetch_path"
cd "$fetch_path" || exit

#dl_fq1=${fetch_path}/${srr}".sra_1.fastq"
#dl_fq2=${fetch_path}/${srr}".sra_2.fastq"
#dl_fqse=${fetch_path}/${srr}".sra.fastq"
dl_fq1=${fetch_path}/${srr}".1.fastq"
dl_fq2=${fetch_path}/${srr}".2.fastq"
dl_fqse=${fetch_path}/${srr}".fastq"


fq1=${fetch_path}/${root}"_R1_001.fastq"
fq2=${fetch_path}/${root}"_R2_001.fastq"
fqse=${fetch_path}/${root}"_R1_001.fastq"

if [ -f "${fq1}".gz ] && [ -f "${fq2}".gz ]; then
 echo final fastqs are present, will not rerun anything.;
 exit 0;
elif [ -f "${fqse}".gz ]; then
 echo final fastq is present, will not rerun anything.;
 exit 0;
fi

if [ -d "${fetch_path}"/"${srr}" ]; then
  echo found prefetch for "$srr", will not rerun prefetch.;
else
  prefetch --max-size 500G -O "$fetch_path"  "$srr"
fi

if [ -f "${dl_fq1}" ] && [ -f "${dl_fq2}" ]; then
  echo found dumped fastqs, will not rerun fasterq-dump.
else
  if [ -f "${dl_fqse}" ]; then
    echo found dumped fastq, will not rerun fasterq-dump.
  else
    fasterq-dump -O "$fetch_path" --split-3 "${fetch_path}"/"${srr}"/"${srr}".sra --threads $ncores
  fi
fi

if [ -f "${dl_fq1}" ]; then
  mv "${dl_fq1}" "${fq1}"
  pigz -p $ncores "${fq1}"
  ${fq1}.gz finised!
fi

if [ -f "${dl_fq2}" ]; then
  mv "${dl_fq2}" "${fq2}"
  pigz -p $ncores "${fq2}"
  ${fq2}.gz finised!
fi

if [ -f "${dl_fqse}" ]; then
  mv "${dl_fqse}" "${fqse}"
  pigz -p $ncores "${fqse}"
  ${fqse}.gz finised!
fi

