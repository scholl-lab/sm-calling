#!/bin/bash
#
#SBATCH --job-name=sm_mutect2_main_job
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --mem=2000M
#SBATCH --output=slurm_logs/%x-%j.log

################################################################################
# Usage:
#   sbatch submit_mutect2_pipeline.sh [SNAKEMAKE_FILE] [CONFIG_FILE] [MAX_JOBS]
#
#   - SNAKEMAKE_FILE:   Path to the Snakemake workflow file
#                       (default: "mutect2_germline_pipeline.smk")
#   - CONFIG_FILE:      Path to a Snakemake config file (default: "config.yaml")
#   - MAX_JOBS:         Number of Snakemake jobs (in parallel) to use (default: 20)
#
################################################################################

# ------------------------------------------------------------------------------
# 1) Parse command-line arguments with defaults
# ------------------------------------------------------------------------------
SNAKEMAKE_FILE=${1:-"mutect2_germline_pipeline.smk"}
CONFIG_FILE=${2:-"config.yaml"}
MAX_JOBS=${3:-20}

# ------------------------------------------------------------------------------
# 2) HPC environment setup
# ------------------------------------------------------------------------------
export TMPDIR="$HOME/scratch/tmp"
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Create the slurm_logs directory if it doesn't exist
mkdir -p slurm_logs

# Export default SBATCH outputs (optional; some clusters auto-handle this)
export SBATCH_DEFAULTS="--output=slurm_logs/%x-%j.log"

echo "DEBUG: Running mutect2 pipeline with:"
echo "       Snakefile:    $SNAKEMAKE_FILE"
echo "       Config file:  $CONFIG_FILE"
echo "       Max jobs:     $MAX_JOBS"
echo "       TMPDIR:       $TMPDIR"

date

# ------------------------------------------------------------------------------
# 3) Run Snakemake (e.g. via srun), using conda environment & cluster profile
# ------------------------------------------------------------------------------
srun snakemake \
    -s "$SNAKEMAKE_FILE" \
    --use-conda \
    --profile=cubi-v1 \
    -j "$MAX_JOBS" \
    --configfile "$CONFIG_FILE"

date
