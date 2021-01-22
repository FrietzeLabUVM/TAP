# Basic Usage

These instructions are for people who are running the pipeline "as is" and have no need to modify it.

## Setup your PATH

PATH controls where linux searches for programs when you try to run them.  The lab maintains a directory with all executables needed to run this pipeline.

The "set it and forget it" way to do this is to add this line to your ~/.bashrc file at the end.

`export PATH=/gpfs2/pi-sfrietze/bin:$PATH`

Then `source .bashrc` and your PATH is setup.  (.bashrc gets sourced everytime you login to the server)

Alternatively you could run `PATH=/gpfs2/pi-sfrietze/bin:$PATH` everytime before running the pipeline.

## Create a configuration file

The configuration file controls all parameters used to run the pipeline.  An example config file is at `/gpfs2/pi-sfrietze/scripts/vacc_rnaseq_pipeline/configs/AI_RNAseq_config.csv`.

Configuration files are comma delimited and can contain comment lines starting with `#` and special comments starting with `#CFG` that supply parameters to the pipeline.  It is possible and recommended that you run the pipeline solely using parameters supplied by configuration files and preserve the configs as a record.

Here's an example config:

\#CFG -i /gpfs2/pi-sfrietze/data/AI_RNAseq -o /gpfs2/pi-sfrietze/data/AI_RNAseq_processed -ref /gpfs2/pi-sfrietze/references/MM10 -SE
31-wt-t0_S93_L006_R1_001.fastq.gz,WT_U0_R1
sfrietze_spring16_20160413_mouse_76_R1.fastq.gz,WT_U0_R2
PR8_S24_L007_R1_001.fastq.gz,WT_U0_R3

The #CFG supplies the input location (-i) output destination (-o) and genome reference (-ref) as well as specifying the data is single end (-SE).

-i, -o, and -ref are always necessary to run the pipeline.  There are additional optional parameters, like -SE in this example that are documented in detail later.

The other 3 lines specify R1 fastq files to use as input and the final output prefix.  The prefix is actually optional but highly recommended to generate "pretty" file names that are easier to use downstream. Easy to use filenames use underscores to separate descriptive variables - with all files specifying the same number of variables in the same order.

If you have multiple sequencing fastqs for a single sample, include all the fastqs on a single line separated by `&`.  For example a line pooling 3 samples would be `fq1.gz&fq2.gz&fq3.gz,sample1`.

If this had been paired-end data the R2 files would be guessed.  See details for --f1_suffix and --f2_suffix for details and limitations.

## Run the pipeline

With our PATH setup and configuration file created, you are ready to run.

The pipeline submission script is `/gpfs2/pi-sfrietze/scripts/vacc_rnaseq_pipeline/submit_rnaseq_pipeline.sh`

Running is simple now:

`bash /gpfs2/pi-sfrietze/scripts/vacc_rnaseq_pipeline/submit_rnaseq_pipeline.sh -c YOUR_CONFIG.csv`

The pipeline submits several SLURM jobs for every sample specified.

It's a good idea to check on your submission periodically.  This command will show your current jobs `squeue -u $(basename ~)`.  If you have jobs showing up in the NODELIST(REASON) column with (DependencyNeverSatisfied), something has gone wrong.  You'll need to `scancel` those jobs and go through the logs in the output folder to figure out what has gone wrong.
