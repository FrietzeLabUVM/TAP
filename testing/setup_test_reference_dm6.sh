#!/bin/bash
if [ -f dm6.fa ]; then
 wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
 gunzip dm6.fa.gz
fi
if [ -f dm6.ensGene.gtf.gz ]; then
 wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/genes/dm6.ensGene.gtf.gz
fi
bash setup_scripts/setup_new_reference.sh -o references/dm6 --gtf_ensemble dm6.ensGene.gtf.gz -f dm6.fa --genomeSAindexNbases 12
