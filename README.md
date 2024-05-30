# MuTect2 Variant Calling Pipeline

This pipeline is designed to perform variant calling using GATK's MuTect2. It is implemented using Snakemake and allows for a high degree of customization through a configuration file.

## Requirements

- Snakemake
- Conda
- GATK (installed through Conda in the pipeline)

## Installation

1. Clone this repository to your local system.
2. Ensure you have all the required software installed.

## Configuration

Before running the pipeline, you need to configure the following files:

### `config.yaml`

This file contains the settings for the pipeline. Here are the details of the settings that you can configure:

- `final_bam_folder`: The folder containing the final BAM files.
- `final_bam_file_extension`: (Optional) The extension of the BAM files (default: ".bam").
- `output_folder`: The folder where the output files will be stored.
- `reference_unpacked`: The path to the reference genome file.
- `panel_of_normals`: The path to the Panel of Normals file.
- `af_only_gnomad`: The path to the allele frequency only gnomAD file.
- `mutect_scatter_by_chromosome`: (Optional) Set to `True` to enable scattering by chromosome, `False` otherwise (default: `False`).

### `calling_metadata.tsv`

This file contains the metadata for the analyses to be run. It should contain the following columns:

- `sample1`: The name of the first sample (tumor sample).
- `sample2`: (Optional) The name of the second sample (normal sample). Leave empty for tumor-only analyses.
- `bam1_file_basename`: The basename of the BAM file for the first sample.
- `bam2_file_basename`: (Optional) The basename of the BAM file for the second sample. Leave empty for tumor-only analyses.
- `individual1`: The identifier for the first individual.
- `individual2`: (Optional) The identifier for the second individual. Leave empty for tumor-only analyses.
- `analysis`: The type of analysis to be performed (e.g., "To" for tumor-only).

## Running the Pipeline

To run the pipeline, use the following command:

```sh
sbatch run_mutect2_calling.sh
```

The run_mutect2_calling.sh shell script contains the Snakemake command to run the workflow with the appropriate settings and resource allocations.
You may need to edit this script to specify the number of cores and other resources based on your system's configuration.

## Output
The pipeline produces the following outputs in the output_folder specified in the config.yaml:

- `variant_calls`: A folder containing the VCF files with the variant calls.
- `logs`: A folder containing the log files for the MuTect2 runs.

Each VCF file is named with the format `<individual1>_<analysis>_<chromosome>.vcf.gz`.

# TODO
- [x] script to merge VCF files, f1r2 files, and stats files
    - use GatherVcfs (Picard) to merge VCF files
    - use gatk MergeMutectStats to merge stats files
    - the f1r2 files are not merged but are all used as input for LearnReadOrientationModel
- [x] script for LearnReadOrientationModel
    - this should be part of the merge script
- [ ] script for FilterMutectCalls
- [x] script for CalculateContamination (plus GetPileupSummaries)

--> see: https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2