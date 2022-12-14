# TAP 

Transcriptomic Analysis Pipeline

A variation of this pipeline is [available for ChIPseq](https://github.com/FrietzeLabUVM/TAP/tree/chipseq)

# Basic Usage

These instructions are for people who are running the RNA-seq and ChIP-seq pipelines "as is" and have no need to modify them. The ChIP-seq pipeline is an elaboration on the RNA-seq pipeline and intructions here generally apply to both.  Please read the ChIP-seq specific section for what is different and additional instructions.

It may also be helpful to consult Sophie and Alyssa's tips document [here](https://docs.google.com/document/d/1t8kc3DhJnxHebRJB4b6FNDKfz_HYsP8HyzlB_1A_5I0/view "Useful tips").

Cong has similarly added notes on running the ChIPseq pipeline [here](https://github.com/FrietzeLabUVM/vacc_rnaseq_pipeline/blob/main/VACC_ChIP-Seq.md)

## Setup your PATH

PATH controls where linux searches for programs when you try to run them.  The lab maintains a directory with all executables needed to run this pipeline.

The "set it and forget it" way to do this is to add this line to your ~/.bashrc file at the end.

`export PATH=/gpfs2/pi-sfrietze/bin:$PATH`

Or you could just run this line once to modify your .bashrc:

`echo 'export PATH=/gpfs2/pi-sfrietze/bin:$PATH' >> ~/.bashrc`

Then `source .bashrc` and your PATH is setup.  (.bashrc gets sourced everytime you login to the server)

Alternatively you could run `PATH=/gpfs2/pi-sfrietze/bin:$PATH` everytime before running the pipeline.

## Create a configuration file

The configuration file controls all parameters used to run the pipeline.  An example config file is at `/gpfs2/pi-sfrietze/scripts/vacc_rnaseq_pipeline/configs/AI_RNAseq_config.csv`.

Configuration files are comma delimited and can contain comment lines starting with `#` and special comments starting with `#CFG` that supply parameters to the pipeline.  It is possible and recommended that you run the pipeline solely using parameters supplied by configuration files and preserve the configs as a record.

Here's an example config:

    #CFG -i /gpfs2/pi-sfrietze/data/AI_RNAseq -o /gpfs2/pi-sfrietze/data/AI_RNAseq_processed -ref /gpfs2/pi-sfrietze/references/MM10 -SE
    #comments like this aren't interpeted as configuration
    31-wt-t0_S93_L006_R1_001.fastq.gz,WT_U0_R1
    sfrietze_spring16_20160413_mouse_76_R1.fastq.gz,WT_U0_R2
    PR8_S24_L007_R1_001.fastq.gz,WT_U0_R3

The #CFG supplies the input location (-i) output destination (-o) and genome reference (-ref) as well as specifying the data is single end (-SE).

-i, -o, and -ref are always necessary to run the pipeline.  There are additional optional parameters, like -SE in this example that are documented in detail later.

The other 3 lines specify R1 fastq files to use as input and the final output prefix.  The prefix is actually optional but highly recommended to generate "pretty" file names that are easier to use downstream. Easy to use filenames use underscores to separate descriptive variables - with all files specifying the same number of variables in the same order.

If you have multiple sequencing fastqs for a single sample, include all the fastqs on a single line separated by `&`.  For example a line pooling 3 samples would be `fq1.gz&fq2.gz&fq3.gz,sample1`.

If this had been paired-end data the R2 files would be guessed.  See details for --f1_suffix and --f2_suffix for details and limitations.

### ChIPseq configuration files

Processing ChIPseq data is more complicated RNAseq files as it requires knowing how to pool replicates and how each samples pairs up to a control or input sample.

Unlike RNAseq, the second column of prefixes is required and an additional 2 columns specifc to ChIPseq are also required: the pooled ChIP prefixes and pooled input prefixes.

The body of a ChIPseq confiugration file could look like this:

    #1) fastq file, 2) unique sample prefix, 3) pooled sample prefix, 4) pooled input prefix to compare to
    fq1a.gz&fq1b.gz,MCF7_ATAD2B_rep1,MCF7_ATAD2B_pooled,MCF7_input_pooled
    fq2.gz,MCF7_ATAD2B_rep2,MCF7_ATAD2B_pooled,MCF7_input_pooled
    fq3.gz,MCF7_input_rep1,MCF7_input_pooled,MCF7_input_pooled
    fq4.gz,MCF7_input_rep2,MCF7_input_pooled,MCF7_input_pooled

Note how the input samples have identical entries for the pooled sample prefix and input prefix

This file would result in 3 sets of peak calls: ATAD2B rep1 vs pooled input, ATAD2B rep2 vs pooled input, and ATAD2B pooled vs pooled input

## Run the pipeline

With our PATH setup and configuration file created, you are ready to run.

The pipeline submission script is `/gpfs2/pi-sfrietze/scripts/vacc_rnaseq_pipeline/submit_rnaseq_pipeline.sh`

The ChIPseq submission script is `/gpfs2/pi-sfrietze/scripts/vacc_chipseq_pipeline/submit_chipseq_pipeline.sh`

Running is simple now:

`bash /gpfs2/pi-sfrietze/scripts/vacc_rnaseq_pipeline/submit_rnaseq_pipeline.sh -c YOUR_CONFIG.csv`

The pipeline submits several SLURM jobs for every sample specified.

It's a good idea to check on your submission periodically.  This command will show your current jobs `squeue -u $(basename ~)`.  If you have jobs showing up in the NODELIST(REASON) column with (DependencyNeverSatisfied), something has gone wrong.  You'll need to `scancel` those jobs and go through the logs in the output folder to figure out what has gone wrong.

Files will appear in the output location as jobs finish.  When all jobs for a sample finish successfully, a \*.complete file gets written.  This \*.complete file will prevent the pipeline from running again for that same sample.  Therefore you can add files to the same configuration file later and resubmit without wasting time reprocessing the same data again.  Or resubmit in case some samples processing OK and other do not.  If you do want to replace a sample, for example if you have done more sequencing or made an error, you will have to manually delete this \*.complete file.

# Parameters

## Special (may replace parameters that are otherwise Required)
-c, --config				
: A valid configuration file, each line of which specifies R1 fastq files, optionally comma delimited with a second entry for final file name prefix. A # commented header is also allowed with lines that start with #CMD being parsed for command line parameters.  The config file can specify all required and optional parameters in this way.

-ref, --reference
: Path to parent directory for all reference components. I have several references setup at `/gpfs2/pi-sfrietze/indexes_jrb`.  They are for HG38, MM10, and DM6.  The "canon" variations don't include haplotypes or patches, just the basic somatic, sex, and mitochondrial chromosomes.  I generally use canon and these are the default if you don't specify "full".

### Building your own reference

To build your own reference, you will need a FASTA file and a correpsponding gene annotation GTF file.  Use the included setup script `bash setup_scripts/setup_new_reference.12core.sh -o OUTPUT_LOCATION -f YOUR_FASTA --gtf_gencode GENCODE GTF`.  ENSEMBLE and UCSC GTF files are also supported via --gtf_ensemble or --gtf_ucsc. Similarly `setup_new_rDNA.sh` can be used to install an rDNA  reference.
  
## Required
-o, --outDir
: Relative or absolute path where pipeline results shuold be written.

-idx, --starIndex
: Path to STAR index, will be derived from -ref if not supplied.

-s, --suppaRef
: Path to SUPPA2 ioi and ioe references, will be derived from -ref if not supplied.

-g, --gtf
: Path to GTF reference file, will be derived from -ref if not supplied.

-fa, --fasta
: Path to genome fasta file, will be derived from -ref if not supplied.

## Optional
-f1s, --f1_suffix
: <_R1_001.fastq.gz> The suffix used for R1 fastq files.  This suffix will be replaced with the R2 suffix to guess R2 fastq files.  If no final file prefix is supplied, removal of the R1 suffix generates the final file prefix.

-f2s, --f2_suffix
: <_R2_001.fastq.gz> The suffix used for R2 fastq files.  Will replace the R1 suffix when guessing R2 files.

-i, --inDir
: <current directory> The directory in which all fastq files are located.

-rDNA, --rDNA_starIndex
: Path to STAR index for organism's rDNA, will be derived from -ref if not supplied.  Without a rDNA STAR index, the rDNA alignment step is skipped.

-SE, --SE
: If activated, alignment will be in single-end mode instead of the RNA-seq default of paired-end.  For ChIP-seq, -SE is the default.

-noSub, --noSub
: If activated, bash will be used to run all pipeline steps in serial instead of sbatch to run in parallel via the job scheduler.  For debugging only or if SLURM's sbatch is not available.
