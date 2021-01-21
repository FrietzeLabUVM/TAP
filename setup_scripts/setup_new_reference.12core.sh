#!/bin/bash
#SBATCH --nodes=1                               # Request one core
#SBATCH --ntasks-per-node=1                               # Request one node (if you request more than one core with -n, also using
#SBATCH --cpus-per-task=12                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p bluemoon                           # Partition to run in
#SBATCH --mem=32000                        # Memory total in MB (for all cores)
#SBATCH -o ref_build_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e ref_build_%j.err                 # File to which STDERR will be written, including job ID


while [[ "$#" -gt 0 ]]; do
  case $1 in
    -o|--output) out="$2"; shift ;;
    -gens|--gtf_ensemble) gtf="$2"; suppa_p="-p"; shift ;;
    -gucsc|--gtf_ucsc) gtf="$2"; suppa_p="-p"; shift ;;
    -genc|--gtf_gencode) gtf="$2"; suppa_p=""; shift ;;
    -f|--fasta) fa="$2"; shift ;;
    --genomeSAindexNbases) saIndex="$2"; shift ;;
    *) echo "Unknown parameter passed: $1"; exit 1 ;;
  esac
  shift
done

if [ -z $out ] || [ -z $gtf ] || [ -z $fa ]; then echo output \(-o\|--output\) gtf \(--gtf_ensemble\|--gtf_ucsc\|--gtf_gencode\) and fasta \(-f\|--fasta\) are all required. quit.; exit 1; fi

echo out will be $out
echo gtf will be $gtf
echo fasta will be $fa
echo suppa_p will be $suppa_p \(ok if empty\)

#if [ -d $out ]; then echo "output location $out exists! delete to or change location. quit."; exit 1; fi
if [ ! -f $gtf ]; then echo gtf file $gtf not found! quit.; exit 1; fi
if [ ! -f $fa ]; then echo fasta file $fa not found! quit.; exit 1; fi
if [ -z $saIndex ]; then saIndex=14; fi

echo saIndex will be $saIndex \(set below 14 for smaller genomes\)


gtf=$(readlink -f $gtf)
fa=$(readlink -f $fa)

mkdir -p $out
out=$(readlink -f $out)
cd $out

bfa=$(basename $fa)
bgtf=$(basename $gtf .gz)

echo staging gtf and fasta...

#gtf
mkdir GTF
cd GTF
if [ $gtf = ${gtf/%.gz/""} ]; then #not a gzipped file
  cp $gtf $bgtf
else #is gzipped
  gunzip -c $gtf > $bgtf
fi
ln -sf $bgtf current.gtf
cd ..

#fasta
mkdir FASTA
cd FASTA
cp $fa $bfa
ln -sf $bfa genome.fa

tfa=$bfa
tfa=${tfa/%.fasta/""}
tfa=${tfa/%.fa/""}
if [ $tfa = $bfa ]; then echo $bfa is named weird, should end in .fa or .fasta, quit; exit 1; fi
tfa=${tfa}.transcriptome.fa

echo generate transcriptome fasta...

gffread -F -w $tfa -g $bfa ../GTF/current.gtf
ln -sf $tfa transcriptome.fa
cd ..

echo making SUPPA2 reference...

#suppa2
SUPPA_BIN=~/lab_bin
suppa="$SUPPA_BIN/python $SUPPA_BIN/suppa.py"

dest=SUPPA2
mkdir $dest

o=${dest}/"suppa2."$(basename ${gtf})"."

#use -p for ensembl and UCSC, not gencode
$suppa generateEvents -i GTF/current.gtf -o $o  -f ioi $suppa_p


for e in SE SS MX RI FL; do
  $suppa generateEvents -i GTF/current.gtf -o $o  -f ioe -e $e $suppa_p
done

echo making STAR index...

#star index
STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir $out/STAR_INDEX \
--genomeFastaFiles FASTA/genome.fa \
--sjdbGTFfile GTF/current.gtf \
--sjdbOverhang 99 \
--genomeSAindexNbases $saIndex
