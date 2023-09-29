#!/bin/bash
#SBATCH --nodes=1                               # Request one core
#SBATCH --ntasks-per-node=1                               # Request one node (if you request more than one core with -n, also using
#SBATCH --cpus-per-task=12                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p bluemoon                           # Partition to run in
#SBATCH --mem=32000                        # Memory total in MB (for all cores)
#SBATCH -o ref_build_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e ref_build_%j.err                 # File to which STDERR will be written, including job ID

CPUS=12

while [[ "$#" -gt 0 ]]; do
  case $1 in
    -o|--output) out="$2"; shift ;;
    -gens|--gtf_ensemble) gtf="$2"; suppa_p="-p"; shift ;;
    -gucsc|--gtf_ucsc) gtf="$2"; suppa_p="-p"; shift ;;
    -genc|--gtf_gencode) gtf="$2"; suppa_p=""; shift ;;
    -f|--fasta) fa="$2"; shift ;;
    --genomeSAindexNbases) saIndex="$2"; shift ;;
    -docker|--docker) docker="$2"; shift ;;
    -singularity|--singularity) singularity="$2"; shift ;;
    -p|--runThreadN) CPUS="$2"; shift ;;
    *) echo "Unknown parameter passed: $1"; exit 1 ;;
  esac
  shift
done

if [ -z $out ] || [ -z $gtf ] || [ -z $fa ]; then echo output \(-o\|--output\) gtf \(--gtf_ensemble\|--gtf_ucsc\|--gtf_gencode\) and fasta \(-f\|--fasta\) are all required. quit.; exit 1; fi
#only allow 1 container type
if [ -n "$docker" ] && [ -n "$singularity" ]; then
  echo Only 1 of docker or signularity should be set. Quit!
  exit 1
fi
container_type=""
if [ -n "$docker" ]; then container_type="docker"; fi
if [ -n "$singularity" ]; then container_type="singularity"; singularity=$(readlink -f $singularity); fi


echo out will be $out
echo gtf will be $gtf
echo fasta will be $fa
echo suppa_p will be $suppa_p \(ok if empty\)
echo container_type is $container_type

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
mkdir -p GTF
cd GTF
if [ $gtf = ${gtf/%.gz/""} ]; then #not a gzipped file
  cp $gtf $bgtf
else #is gzipped
  gunzip -c $gtf > $bgtf
fi
ln -sf $bgtf current.gtf

awk '$3 != "gene" ' current.gtf > current.no_genes.gtf

cd ..

#fasta
mkdir -p  FASTA
cd FASTA
cp $fa $bfa
ln -sf $bfa genome.fa

tfa=$bfa
tfa=${tfa/%.fasta/""}
tfa=${tfa/%.fa/""}
if [ $tfa = $bfa ]; then echo $bfa is named weird, should end in .fa or .fasta, quit; exit 1; fi
tfa=${tfa}.transcriptome.fa

echo generate transcriptome fasta...

gffread -F -w $tfa -g $bfa ../GTF/current.no_genes.gtf
ln -sf $tfa transcriptome.fa
cd ..

echo making SUPPA2 reference...

#suppa2
GTF=GTF/current.gtf
SUPPA_REF_DIR=SUPPA2
mkdir -p  $SUPPA_REF_DIR

# container for suppa2 v1.1
if [ -n "$container_type" ]; then
  
  dGTF=/input_gtf/$(basename $GTF)
  dSUPPA_REF_DIR=/suppa_ref

  if [ $container_type = "docker" ]; then
    cmd_suppa="docker run \
      -u $(id -u):$(id -g) \
      -v $(dirname $GTF):$(dirname $dGTF) \
      -v $SUPPA_REF_DIR:$dSUPPA_REF_DIR \
      --entrypoint \
      suppa.py $docker"
  elif [ $container_type = "singularity" ]; then
    cmd_suppa="singularity exec \
      --bind $(dirname $GTF):$(dirname $dGTF),$SUPPA_REF_DIR:$dSUPPA_REF_DIR \
      $singularity suppa.py"
  else
    echo "Unrecognized container_type $container_type";
    exit 1;
  fi
  GTF=$dGTF
  SUPPA_REF_DIR=$dSUPPA_REF_DIR
else
  cmd_suppa="suppa.py"
fi

SUPPA_REF=${SUPPA_REF_DIR}/"suppa2."$(basename ${gtf})"."
#use -p for ensembl and UCSC, not gencode
$cmd_suppa generateEvents -i $GTF -o $SUPPA_REF  -f ioi $suppa_p

for e in SE SS MX RI FL; do
  $cmd_suppa generateEvents -i $GTF -o $SUPPA_REF  -f ioe -e $e $suppa_p
done

STAR_IDX_OUT_DIR=$out/STAR_INDEX
mkdir -p $STAR_IDX_OUT_DIR
STAR_FASTA_DIR=$out/FASTA
STAR_GTF_DIR=$out/GTF

echo making STAR index...
if [ -n "$container_type" ]; then
  # Derive mount points for input files inside docker
  #--genomeDir $out/STAR_INDEX \
  #--genomeFastaFiles FASTA/genome.fa \
  #--sjdbGTFfile GTF/current.gtf \
  dSTAR_IDX_OUT_DIR=STAR_INDEX
  dSTAR_FASTA_DIR=FASTA
  dSTAR_GTF_DIR=GTF

   
  echo $container_type STAR_IDX_OUT_DIR is $dSTAR_IDX_OUT_DIR
  echo $container_type STAR_FASTA_DIR is $dSTAR_FASTA_DIR
  echo $container_type STAR_GTF_DIR is $dSTAR_GTF_DIR

  dir_B=$(dirname "$B")
  dir_dB=$(dirname "$dB")

  if [ $container_type = "docker" ]; then
  cmd_star="docker run \
    -u $(id -u):$(id -g) \
    -v $STAR_IDX_OUT_DIR:$dSTAR_IDX_OUT_DIR \
    -v $STAR_FASTA_DIR:$dSTAR_FASTA_DIR \
    -v $STAR_GTF_DIR:$dSTAR_GTF_DIR \
    --entrypoint STAR\
    $docker \
    "
  elif [ $container_type = "singularity" ]; then
  cmd_star="singularity exec \
    --bind $STAR_IDX_OUT_DIR:$dSTAR_IDX_OUT_DIR,$STAR_FASTA_DIR:$dSTAR_FASTA_DIR,$STAR_GTF_DIR:$dSTAR_GTF_DIR \
    $singularity \
    STAR \
    "
  else
    echo "Unrecognized container_type $container_type";
    exit 1;
  fi
  #update STAR command inputs to use docker paths
  STAR_IDX_OUT_DIR=$dSTAR_IDX_OUT_DIR
  STAR_FASTA_DIR=$dSTAR_FASTA_DIR,
  STAR_GTF_DIR=$dSTAR_GTF_DIR
else 
  cmd_star=STAR
fi


echo creating STAR_INDEX
#star index
$cmd_star --runThreadN $CPUS \
--runMode genomeGenerate \
--genomeDir $STAR_IDX_OUT_DIR \
--genomeFastaFiles $STAR_FASTA_DIR/genome.fa \
--sjdbGTFfile $STAR_GTF_DIR/current.gtf \
--sjdbOverhang 99 \
--genomeSAindexNbases $saIndex

chmod a+rx $out/STAR_INDEX
