# sm-calling

Snakemake 8+ variant calling pipeline for somatic (Mutect2) and germline
(FreeBayes) analysis on HPC clusters.

Part of the [scholl-lab](https://github.com/scholl-lab) bioinformatics suite.

## Quick Start

```bash
# 1. Install Snakemake 8+ and the SLURM executor plugin
conda create -n snakemake -c conda-forge -c bioconda 'snakemake>=8.0'
pip install snakemake-executor-plugin-slurm

# 2. Edit configuration
cp config/config.yaml config/my_config.yaml   # customize paths
vim config/samples.tsv                         # add your samples

# 3. Run (auto-detects BIH / Charite / local)
sbatch scripts/run_snakemake.sh config/my_config.yaml
```

## Supported Callers

| Caller | Analysis Types | Output |
|--------|---------------|--------|
| **Mutect2** | tumor_only, tumor_normal | Filtered somatic VCFs per sample |
| **FreeBayes** | germline | Merged + normalized germline VCF |

Set `caller` in `config/config.yaml` to `"mutect2"`, `"freebayes"`, or `"all"`.

## Configuration

### `config/config.yaml`

Unified hierarchical configuration. Key sections:

| Section | Description |
|---------|-------------|
| `caller` | Which caller(s) to run |
| `ref` | Reference genome path and build |
| `gatk_resources` | Panel of normals, gnomAD files (Mutect2) |
| `paths` | Input BAM folder, output folder, samples TSV |
| `bam` | BAM file extension |
| `scatter` | Scatter mode (`chromosome` / `interval` / `none`) |
| `params` | Tool-specific extra arguments |

### `config/samples.tsv`

Tab-separated sample sheet:

```
sample          tumor_bam          normal_bam        analysis_type
IND001_To       IND001.tumor       .                 tumor_only
IND002_TN       IND002.tumor       IND002.normal     tumor_normal
IND003_G        IND003             .                 germline
```

- `sample`: Unique identifier (used in output filenames)
- `tumor_bam`: BAM basename (without extension)
- `normal_bam`: Matched normal basename, or `.` if none
- `analysis_type`: One of `tumor_only`, `tumor_normal`, `germline`

## Running on HPC

The launcher auto-detects the cluster environment:

```bash
# Submit via SLURM (auto-detects BIH or Charite):
sbatch scripts/run_snakemake.sh

# Custom config:
sbatch scripts/run_snakemake.sh config/my_config.yaml

# Dry-run (local):
bash scripts/run_snakemake.sh config/config.yaml -n

# Pass extra Snakemake flags:
sbatch scripts/run_snakemake.sh config/config.yaml --forcerun mutect2_call
```

Or invoke Snakemake directly:

```bash
snakemake \
    -s workflow/Snakefile \
    --configfile config/config.yaml \
    --workflow-profile profiles/default \
    --profile profiles/bih
```

## Pipeline Architecture

```
workflow/Snakefile
    |
    +-- rules/common.smk      Config shortcuts, metadata, helpers
    +-- rules/scatter.smk      Interval scattering (shared)
    +-- rules/mutect2.smk      Mutect2: call -> gather -> contamination -> filter
    +-- rules/freebayes.smk    FreeBayes: call -> merge + normalize
```

### Mutect2 DAG

```
BAM files
    +-> mutect2_call (per sample x scatter_unit)        [temp]
    +-> get_pileup_summaries (per unique BAM)
    +-> gather_mutect2_vcfs (per sample)
    +-> merge_mutect2_stats (per sample)
    +-> learn_read_orientation (per sample)
    +-> calculate_contamination (per sample)
    +-> filter_mutect_calls (per sample)                [protected]
```

### FreeBayes DAG

```
BAM files
    +-> freebayes_call (per scatter_unit)               [temp]
    +-> merge_freebayes_vcfs                            [protected]
```

## Development

```bash
# Lint
make lint

# Format
make format

# Run unit tests
make test-unit

# Run all tests (requires snakemake + conda)
make test
```

## Project Structure

```
sm-calling/
+-- config/              Configuration and sample sheet
+-- workflow/
|   +-- Snakefile        Entry point (Snakemake 8+)
|   +-- rules/           Rule files + helpers.py
|   +-- envs/            Conda environment definitions
|   +-- schemas/          JSON schemas for validation
+-- profiles/            Snakemake profiles (default, bih, charite, local)
+-- scripts/             Universal launcher
+-- tests/               pytest unit and integration tests
+-- deprecated/          Old standalone workflows (preserved for reference)
```

## License

MIT License. See [LICENSE](LICENSE) for details.
