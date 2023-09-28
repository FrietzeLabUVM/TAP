#!/bin/bash
#SBATCH -p short                           # Partition to run in
#SBATCH -J echo_sub
#SBATCH -o submission_params_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e submission_params_%j.err 

date
echo $0 $@

