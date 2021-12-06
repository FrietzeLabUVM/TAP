# Create input and output directories
## In input directory, create symbolic links for archived fastq.gz files (multiple at once)
ln -s fullpath/to/fastq.gz .

# Create a configuration file
## Assume single-end sequencing in default
## Firstly input samples, then chip target samples, each line is for one replicate sample
## sample config csv is here: /users/c/g/cgao1/BRCA_progression/BRCA_H3K9me3/H3K9me3_config.params.csv
## Explanation for configuration could be found here: https://github.com/FrietzeLabUVM/vacc_rnaseq_pipeline
## If you have multiple sequencing fastqs for a single sample, include all the fastqs on a single line separated by &. For example a line pooling 3 samples would be fq1.gz&fq2.gz&fq3.gz,rep_name,pool_name,input_pool_name; 
## or you can merge them in advance with the command line: cat file1.gz file2.gz file3.gz > allfiles.gz

# Submit for ChIP-Seq pipeline
bash /gpfs2/pi-sfrietze/scripts/vacc_chipseq_pipeline/submit_chipseq_pipeline.sh -c YOUR_CONFIG.csv
## The pipeline submits several SLURM jobs for every sample specified.
## Files will appear in the output location as jobs finish. When all jobs for a sample finish successfully, a *.complete file gets written. This *.complete file will prevent the pipeline from running again for that same sample. Therefore you can add files to the same configuration file later and resubmit without wasting time reprocessing the same data again. Or resubmit in case some samples processing OK and other do not. If you do want to replace a sample, for example if you have done more sequencing or made an error, you will have to manually delete this *.complete file.

# Manage jobs
squeue -u username
scancel jobID
scancel -u username

# Transfer files from VACC to galaxy server to do downsteam analysis (https://www.tecmint.com/rsync-local-remote-file-synchronization-commands/)
## Copy a Directory from Local Server to a Remote Server
rsync -avzh /root/rpmpkgs root@192.168.0.141:/root/
## Copy a File from a Remote Server to a Local Server with SSH
rsync -avzhe ssh root@192.168.0.141:/root/anaconda-ks.cfg /tmp
## Select (find) specific files using regex, then transfer (eg. broadPeak files for target H3K9me3)
find . -type f -regex ".+_H3K9me3_.+broadPeak$" | rsync -avzh --files-from=- ./ root@XXX:/path/to/destination/



