    For reference setup
    
    Required:
    -o|--output The output directory for the new index structure.
    There are 3 "styles" of GTF file I've come across. Choose the appropriate one for where you downloaded your gene reference from.  My preference is GENCODE.
      -gens|--gtf_ensemble [GTF from ENSEMBL](https://useast.ensembl.org/)
      -gucsc|--gtf_ucsc [GTF from UCSC](https://genome-euro.ucsc.edu/cgi-bin/hgTables) 
      -genc|--gtf_gencode [GTF from GENCODE](https://www.gencodegenes.org/)
    -f|--fasta The single fasta file containing all genomic sequences
    Optional:
    --genomeSAindexNbases Important for STAR indexing.  Default of 14 is appropriate for mouse/human size genomes. Must be lowered to avoid error for smaller genomes.
    
