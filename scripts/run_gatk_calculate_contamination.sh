#!/bin/bash

#SBATCH --job-name=sm_gatk_contamination_main_job
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=1200M
#SBATCH --output=slurm_logs/%x-%j.log

# Handle temporary files and directories as recommended by the HPC documentation
# https://hpc-docs.cubi.bihealth.org/best-practice/temp-files/#tmpdir-and-the-scheduler
# https://bihealth.github.io/bih-cluster/slurm/snakemake/#custom-logging-directory

# First, point TMPDIR to the scratch in your home as mktemp will use this
export TMPDIR=$HOME/scratch/tmp
# Second, create another unique temporary directory within this directory
export TMPDIR=$(mktemp -d)
# Finally, setup the cleanup trap to remove the temporary directory on exit
trap "rm -rf $TMPDIR" EXIT

# Create a directory for logs if it doesn't exist
mkdir -p slurm_logs
export SBATCH_DEFAULTS=" --output=slurm_logs/%x-%j.log"

# Print the current date and time
date
# Execute the Snakemake script with the desired options
srun snakemake -s snakemake/gatk_calculate_contamination.smk --use-conda --profile=cubi-v1 -j100
# Print the current date and time after completion
date
