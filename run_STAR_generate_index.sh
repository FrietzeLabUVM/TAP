#!/bin/bash

#SBATCH -p short 		# partition name
#SBATCH -t 0-2:00 		# hours:minutes runlimit after which job will be killed
#SBATCH -c 6 		# number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem 16G
#SBATCH --job-name STAR_index 		# Job name
#SBATCH -o star_index.%j.out			# File to which standard out will be written
#SBATCH -e star_index.%j.err 		# File to which standard err will be written

STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir DM6_STAR_INDEX \
--genomeFastaFiles /users/j/r/jrboyd/lab_shared/indexes/DM6/FASTA/dm6.fa \
--sjdbGTFfile /users/j/r/jrboyd/lab_shared/indexes/DM6/GTF/dm6.ensGene.gtf \
--sjdbOverhang 99 \
--genomeSAindexNbases 12
