# Create input and output directories

Inputs must be in a single directory.  If your fastq files are not, do not mv move them. Instead, create symbolic links for archived fastq.gz files.

mkdir input_fastqs
ln -s fullpath/to/\*fastq.gz input_fastqs/

# Create a configuration file

The ChIPseq pipeline assumes single-end sequencing so that is the default (RNAseq assumes paired-end). If you have paired-end ChIPseq (or cut&run or ATAC), set the -PE flag in your config.

    #CFG -PE

In the body of your conifg, every line should a biological replicate. If you have sequencing or technical replicates, combine them using & on a single line.  

    fq1.gz&fq2.gz&fq3.gz,rep_name,pool_name,input_pool_name; 

If you prefer, you could combine sequencing reps ahead of time like so:
    
    cat file1.gz file2.gz file3.gz > pooled.gz

A sample config csv is here: /users/c/g/cgao1/BRCA_progression/BRCA_H3K9me3/H3K9me3_config.params.csv

General config file info is found here: https://github.com/FrietzeLabUVM/vacc_rnaseq_pipeline

# Submit for ChIP-Seq pipeline
bash /gpfs2/pi-sfrietze/scripts/vacc_chipseq_pipeline/submit_chipseq_pipeline.sh -c YOUR_CONFIG.csv

The pipeline submits several SLURM jobs for every sample specified.

Files will appear in the output location as jobs finish. When all jobs for a sample finish successfully, a \*.complete file gets written. This \*.complete file will prevent the pipeline from running again for that same sample. Therefore you can add files to the same configuration file later and resubmit without wasting time reprocessing the same data again. Or resubmit in case some samples processing OK and other do not. If you do want to replace a sample, for example if you have done more sequencing or made an error, you will have to manually delete this \*.complete file. *Do not resubmit the same config without first killing jobs in progress with scancel. Crazy things will happen.*

# Manage jobs

    squeue -u username #see your jobs
    scancel jobID #cancel a single job
    scancel -u username #cancel all of your jobs

# Transfer files from VACC to Galaxy server 

Since we don't have Rstudio via VACC, most analysis is done on Galaxy.  Only copy what you need.  Maybe you just need count files?

Use rsync while logged into the Galaxy via ssh  (https://www.tecmint.com/rsync-local-remote-file-synchronization-commands/)

While logged into Galaxy, copy a Directory from VACC to Galaxy

    rsync -trP YOUR_VACC_USER@vacc-user1.uvm.edu:my/dir/or/file/\*pattern location/on/galaxy
    
While logged into Galaxy, copy a Directory from Galxy to VACC

    rsync -trP my/dir/or/file/\*pattern YOUR_VACC_USER@vacc-user1.uvm.edu:location/on/VACC

You can accomplish the same thing from VACC supplying your Galaxy login instead. YOUR_GALAXY_USER@galaxy.med.uvm.edu or YOUR_GALAXY_USER@galaxy2.med.uvm.edu


