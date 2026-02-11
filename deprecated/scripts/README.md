# Scripts Folder README

This folder contains shell scripts to run the Snakemake workflows for the MuTect2 Variant Calling Pipeline. Each script is designed to launch a specific Snakemake workflow and is configured to run on a SLURM-based cluster.

## Scripts Overview

### `run_mutect2_calling.sh`
- **Description**: Runs the `mutect2_calling.smk` Snakemake workflow.
- **Function**: Calls variants using MuTect2.
- **SLURM Configuration**:
  - Job name: `sm_mutect2_calling_main_job`
  - Time: 168 hours
  - Memory: 1200M per CPU
  - Output logs: `slurm_logs/%x-%j.log`
- **Usage**:
  ```bash
  sbatch run_mutect2_calling.sh
  ```

### `run_merge_mutect2_calls.sh`
- **Description**: Runs the `merge_mutect2_calls.smk` Snakemake workflow.
- **Function**: Merges scattered VCF files, f1r2 files, and stats files.
- **SLURM Configuration**:
  - Job name: `sm_merge_mutect2_calls_main_job`
  - Time: 168 hours
  - Memory: 1200M per CPU
  - Output logs: `slurm_logs/%x-%j.log`
- **Usage**:
  ```bash
  sbatch run_merge_mutect2_calls.sh
  ```

### `run_gatk_calculate_contamination.sh`
- **Description**: Runs the `gatk_calculate_contamination.smk` Snakemake workflow.
- **Function**: Calculates contamination using GATK's CalculateContamination.
- **SLURM Configuration**:
  - Job name: `sm_gatk_calculate_contamination_main_job`
  - Time: 168 hours
  - Memory: 1200M per CPU
  - Output logs: `slurm_logs/%x-%j.log`
- **Usage**:
  ```bash
  sbatch run_gatk_calculate_contamination.sh
  ```

### `run_gatk_filter_mutect_calls.sh`
- **Description**: Runs the `gatk_filter_mutect_calls.smk` Snakemake workflow.
- **Function**: Filters the merged raw MuTect2 variant calls.
- **SLURM Configuration**:
  - Job name: `sm_gatk_filter_mutect_calls_main_job`
  - Time: 168 hours
  - Memory: 1200M per CPU
  - Output logs: `slurm_logs/%x-%j.log`
- **Usage**:
  ```bash
  sbatch run_gatk_filter_mutect_calls.sh
  ```

### `run_bcftools_concat.sh`
- **Description**: Concatenates VCF files using bcftools.
- **Function**: Combines VCF files from different intervals into a single VCF file.
- **SLURM Configuration**:
  - Job name: `sm_bcftools_concat_main_job`
  - Time: 168 hours
  - Memory: 1200M per CPU
  - Output logs: `slurm_logs/%x-%j.log`
- **Usage**:
  ```bash
  sbatch run_bcftools_concat.sh
  ```

## Shell Script Example

Below is an example structure of a typical run script (`run_mutect2_calling.sh`):

```bash
#!/bin/bash
#
#SBATCH --job-name=sm_mutect2_calling_main_job
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=1200M
#SBATCH --output=slurm_logs/%x-%j.log

# Setup temporary directory
export TMPDIR=$HOME/scratch/tmp
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

mkdir -p slurm_logs
export SBATCH_DEFAULTS=" --output=slurm_logs/%x-%j.log"

date
srun snakemake -s snakemake/mutect2_calling.smk --use-conda --profile=cubi-v1 -j100
date
```
