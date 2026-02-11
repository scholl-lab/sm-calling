# Snakemake Scripts Folder

This folder contains Snakemake workflows for the MuTect2 Variant Calling Pipeline. Each script is responsible for a specific part of the variant calling process, from calling variants to filtering and merging results. Below is a detailed description of each script and its functionality.

## Scripts Overview

### `mutect2_calling.smk`
- **Description**: This script calls variants using MuTect2.
- **Function**: It processes BAM files to generate VCF files with variant calls for each individual and analysis type.
- **Main Rules**:
  - `all`: Specifies the final output files.
  - `call_variants`: Runs MuTect2 to call variants for each chromosome.

### `merge_mutect2_calls.smk`
- **Description**: This script merges scattered VCF files, stats files, and f1r2 files.
- **Function**: It combines the results of variant calling into a single VCF file per individual and analysis type.
- **Main Rules**:
  - `all`: Specifies the final merged output files.
  - `merge_vcfs`: Merges VCF files across chromosomes.
  - `merge_stats`: Merges stats files corresponding to VCFs.
  - `merge_f1r2`: Merges f1r2 orientation files across chromosomes.

### `gatk_calculate_contamination.smk`
- **Description**: This script calculates contamination using GATK's CalculateContamination.
- **Function**: It calculates the fraction of reads coming from cross-sample contamination.
- **Main Rules**:
  - `all`: Specifies the final contamination table output files.
  - `get_pileup_summaries`: Generates pileup summaries needed for contamination calculation.
  - `calculate_contamination`: Runs the CalculateContamination tool to estimate contamination levels.

### `gatk_filter_mutect_calls.smk`
- **Description**: This script filters the merged raw MuTect2 variant calls.
- **Function**: It processes the merged VCF files to produce filtered VCF files with high-confidence variant calls.
- **Main Rules**:
  - `all`: Specifies the final filtered VCF output files.
  - `filter_mutect_calls`: Runs the FilterMutectCalls tool to filter variant calls based on various criteria.

### `bcftools_concat.smk`
- **Description**: Concatenates VCF files using bcftools.
- **Function**: Combines VCF files from different intervals into a single VCF file.
- **Main Rules**:
  - `all`: Specifies the final concatenated VCF output files.
  - `concatenate_vcfs`: Uses bcftools to concatenate VCF files and index the result.

## Usage

To run any of the above scripts, use the corresponding `run_*.sh` shell script, which sets up the environment and submits the Snakemake workflow to a SLURM-based cluster. For example, to run the `mutect2_calling.smk` script, use:

```bash
sbatch run_mutect2_calling.sh
```
