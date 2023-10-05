## TAP 

Welcome to the github page for TAP (Transcriptomic Analysis Pipeline).

## What is TAP?

TAP is a pipeline for processing RNA-seq datasets, starting from gzipped fastq files (fastq.gz) and resulting in several useful output files. Such as:

1) sorted aligned reads (.bam)
2) gene count files (.tab)
3) TPM (transcripts per million) (.txt)
4) SNV quantification (.vcf)
5) bigWig (.bw) files for visualization in genome browser tracks.

In short, these are the outputs you would need for almost any analysis; diffential gene expression, differential splicing, sequence variants, etc.

![schematic of TAP pipeline steps.](https://github.com/FrietzeLabUVM/TAP/blob/main/TAP_Figure_1.png?raw=true)

## Why TAP?

There are many RNA-seq pipelines and workflows available at this point. If you're happy with your current setup for processing RNA-seq data, frankly, you probably don't need TAP. 

TAP aims to make setup of genomic references and software dependencies as simple as possible. I've setup TAP on several systems at this point, from my daily-use Windows machine to my university's HPCC (high performance computing cluster). 

Dependencies are intended to be handled via either docker or singularity, although there are instructions for manual installation if you insist.

All references are generated by a single script that requires 2 commonly available inputs. A genomic sequence fasta (.fasta/.fa) and a gene reference file (.gtf/.gff).

# Basic Usage

1) Clone the TAP repository
   * `git clone git@github.com:FrietzeLabUVM/TAP.git; cd TAP`
3) Deploy job scripts.
   * ``
   * some extra steps are required if you plan to [use SLURM](### SLURM)
5) Obtain dependencies.
6) Install references.
7) Create a configuration file.

## Deployment

Deployment is alawys required before running TAP for the first time. Deployment utitlizes a configuration file and template job scripts to tailor TAP to the hardware of your system.

There are additional critical details if you'll be using SLURM.

### SLURM

With SLURM, you also need to consider the job header files.  The content of these files gets prepended to job templates and instructs SLURM how to handle each job. The most important things to look at are below:

* --cpus-per-task : the number of threads to use, should match the "processors" column of your deployment config.
* -t : max time the job is allowed to run. this depends of the partition, use `sinfo` to check.
* -p : partition name, again use `sinfo` to see options.
* --mem : memory in MB the job is allowed to use.

## Dependencies

The easiest way to handle all dependencies is via either [Docker](https://www.docker.com/) or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html). With either of these approaches you can even run this pipeline on Windows or Mac.

That said, you may already have all or most of the dependencies installed on your system, or simply prefer to install things the old fashioned way. You can certainly do so, the only requirement is that all executables be on your PATH at the time the pipeline is run.

### Docker

Satisfying all dependencies using Docker is as simple as:

```bash
docker pull jrboyd/tap
```

That's it!

Now you need to let TAP know you're using Docker by appending ```--docker jrboyd/tap``` to all submission calls or including ```--docker jrboyd/tap``` in the configuration header.

### Singularity

Like Docker, Singularity make it simple to satisfy dependencies:

```bash
singularity pull docker://jrboyd/tap
```

And also similar to Docker, you need to let TAP know you're using Singularity by appending ```--singularity tap_latest.sif``` to all submission calls or including ```--singularity tap_latest.sif``` in the configuration header. Where tap_latest.sif is the path to the .sif file generated by ```singularity pull```.

### Conventional Installation

See the included script for an outline of installing all dependencies: `setup_scripts/manual_install_dependencies.sh`. You will need to elevate to sudo privleges for some of these installation steps as written. It's certainly possible, without sudo, to install all dependencies on any system where you're already able to run other bioinformatics software. I've done it. You'll need to make adjustments as necessary.

The PATH variable controls where Linux searches for programs when you try to run them. Check it with ```echo $PATH```.  

If you happen to be a user of the VACC HPC at UVM, the Frietze lab maintains a directory with all executables needed to run this pipeline: ```/gpfs2/pi-sfrietze/bin```.

The "set it and forget it" way to do this is to add this line to your ~/.bashrc file at the end.

`export PATH=/gpfs2/pi-sfrietze/bin:$PATH`

Then `source .bashrc` and your PATH is setup.  (.bashrc gets sourced everytime you login to the server)

Alternatively you could run `PATH=/gpfs2/pi-sfrietze/bin:$PATH` every time before running the pipeline.

## Create a configuration file

The configuration file can potentially control all parameters used to run the pipeline.  An example config file is [here](https://github.com/FrietzeLabUVM/TAP/blob/main/testing/test_configs/test_dm6_config.params.csv)`. You may also supply all pipeline parameters when calling the submit script. If you have parameters defined in the config and in the script call, parameters in the script call will overwrite duplicated parameters in the config file.

Configuration files are comma delimited and can contain comment lines starting with `#` and special comments starting with `#CFG` that supply parameters to the pipeline.  It is possible and recommended that you run the pipeline solely using parameters supplied by configuration files and preserve the configs as a record.

Here's an example config:

    #CFG -i AI_RNAseq_fastqs -o AI_RNAseq_processed -ref ~/references/MM10 -SE
    #comments like this aren't interpeted as configuration
    31-wt-t0_S93_L006_R1_001.fastq.gz,WT_U0_R1
    sfrietze_spring16_20160413_mouse_76_R1.fastq.gz,WT_U0_R2
    PR8_S24_L007_R1_001.fastq.gz,WT_U0_R3

The #CFG supplies the input location (-i) output destination (-o) and genome reference (-ref) as well as specifying the data is single end (-SE).

-i, -o, and -ref are always necessary to run the pipeline.  There are additional optional parameters, like -SE in this example that are [documented in detail later](#parameters).

The other 3 lines specify R1 fastq files to use as input and the final output prefix.  The prefix is actually optional but highly recommended to generate "pretty" file names that are easier to use downstream. Easy to use filenames use underscores to separate descriptive variables - with all files specifying the same number of variables in the same order.

If you have multiple sequencing fastqs for a single sample, include all the fastqs on a single line separated by `&`.  For example a line pooling 3 samples would be `fq1.gz&fq2.gz&fq3.gz,sample1`.

Only R1 files should be specfied. If this had been paired-end data the R2 files would be derived be substituting the R2 suffix with the R1 suffix.  See parameters [--f1_suffix and --f2_suffix](#optional) for details and limitations.

## Run the pipeline

With your dependecies available, references generated, and configuration file created, you are ready to run.

The pipeline submission script is `submit_rnaseq_pipeline.sh`

Running is simple now:

`bash submit_rnaseq_pipeline.sh -c YOUR_CONFIG.csv`

The pipeline submits several SLURM jobs for every sample specified.

It's a good idea to check on your submission periodically.  This command will show your current jobs `squeue -u $(basename ~)`.  If you have jobs showing up in the NODELIST(REASON) column with (DependencyNeverSatisfied), something has gone wrong.  You'll need to `scancel` those jobs and go through the logs in the output folder to figure out what has gone wrong.

Files will appear in the output location as jobs finish.  When all jobs for a sample finish successfully, a \*.complete file gets written.  This \*.complete file will prevent the pipeline from running again for that same sample.  Therefore you can add fastq files to the same configuration file later and resubmit without wasting time reprocessing the same data again.  Or resubmit in case some samples processing OK and other do not.  If you do want to replace a sample, for example if you have done more sequencing or made an error, you will have to manually delete this \*.complete file and all other impacted outputs.

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
: <_R1_001.fastq.gz> The suffix used for R1 fastq files.  This suffix will be replaced with the R2 suffix to guess R2 fastq files.  If no final file prefix is supplied in column 2, removal of the R1 suffix generates the final file prefix.

-f2s, --f2_suffix
: <_R2_001.fastq.gz> The suffix used for R2 fastq files.  Will replace the R1 suffix when guessing R2 files.

-i, --inDir
: <current directory> The directory in which all fastq files are located.

-rDNA, --rDNA_starIndex
: Path to STAR index for organism's rDNA, will be derived from -ref if not supplied.  Without a rDNA STAR index, the rDNA alignment step is skipped.

-SE, --SE
: If activated, alignment will be in single-end mode instead of the RNA-seq default of paired-end.

-PE, --PE
: If activated, alignment will be in paired-end mode, already the default for RNA-seq.

-noSub, --noSub
: If activated, bash will be used to run all pipeline steps in serial instead of sbatch to run in parallel via the job scheduler.  For debugging only or if SLURM's sbatch is not available.

-docker, --docker
: Name of an installed docker image. Likely jrboyd/tap. Docker will provide all dependencies. `docker` itself still needs to be on your PATH.

-singularity, --singularity
: Path to a singularity .sif file. Likely named tap_latest.sif. Singularity will provide all dependencies. `singularity` itself still needs to be on your PATH.
