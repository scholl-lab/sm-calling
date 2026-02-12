# sm-calling

Snakemake 8+ variant calling pipeline for somatic (Mutect2) and germline (FreeBayes) analysis on SLURM HPC clusters.

**Analysis-ready BAMs → scatter → call → gather → filter → VCF**

Supports both **BIH HPC** and **Charite HPC** (auto-detected), plus local execution.

Part of the [scholl-lab](https://github.com/scholl-lab) bioinformatics suite — designed to run after [sm-alignment](https://github.com/scholl-lab/sm-alignment).

---

## Setup

### 1. Clone into your project directory

The pipeline should live inside each project's directory on the cluster, not in your home directory (quota is too small for results).

```bash
# Charite HPC
cd /sc-projects/<your-project>
git clone https://github.com/scholl-lab/sm-calling.git
cd sm-calling

# BIH HPC
cd /data/cephfs-1/work/projects/<your-project>
git clone https://github.com/scholl-lab/sm-calling.git
cd sm-calling
```

### 2. Install Snakemake 8+ and the SLURM executor plugin

```bash
conda create -n snakemake -c conda-forge -c bioconda 'snakemake>=8.0'
conda activate snakemake
pip install snakemake-executor-plugin-slurm
```

### 3. Set up reference data

Follow the [lab handbook: Reference Data Setup](https://github.com/scholl-lab/lab-handbook/blob/main/docs/reference-data-setup.md).

**BIH HPC** — shared data already exists:

```
/data/cephfs-1/work/groups/scholl/shared/ref/GRCh38/
/data/cephfs-1/work/projects/apa-sequencing/analysis/GATK_resource_bundle/
```

**Charite HPC** — set up per project:

```
/sc-projects/<your-project>/resources/ref/GRCh38/
/sc-projects/<your-project>/resources/gatk_bundle/
```

### Recommended project layout

```
/sc-projects/<your-project>/
├── sm-alignment/                    # alignment pipeline (upstream)
├── sm-calling/                      # this pipeline (git clone)
│   ├── workflow/
│   ├── config/
│   │   ├── config.yaml              # ← edit or generate for your project
│   │   └── samples.tsv              # ← generate from BAM filenames
│   ├── profiles/
│   └── scripts/
├── resources/                       # reference data (or symlinks to shared)
│   ├── ref/GRCh38/
│   └── gatk_bundle/
└── results/
    └── exomes/
        ├── bqsr/                    # BAMs from sm-alignment (input)
        └── variant_calls/           # outputs from sm-calling
```

---

## Generate Config Files

The config generator scans your BAM directory, classifies tumor/normal roles by filename, pairs tumors with normals, and writes `config/samples.tsv` and `config/config.yaml`.

### Interactive wizard (recommended for first-time setup)

Run without arguments for a guided 3-step workflow:

```bash
python scripts/generate_config.py
```

```
============================================================
  sm-calling -- Config Generator
============================================================

BAM folder [results/exomes/bqsr]: /data/bams
BAM extension [.merged.dedup.bqsr.bam]:

Found 3 BAMs:

  #  Basename              Size      Role     Reason
  1  SAMPLE01-T           45.3 GB   tumor    suffix '-T'
  2  SAMPLE01-N           19.2 GB   normal   suffix '-N'
  3  SAMPLE02_tumor        8.1 GB   tumor    suffix '_tumor'
```

The wizard has three phases:

1. **Propose** — scan BAMs, auto-detect roles (tumor/normal) from filename suffixes (`-T`, `_N1`, `_tumor`, etc.), prompt only for unrecognized BAMs
2. **Review** — show the complete proposed table, accept or edit individual rows, switch analysis mode
3. **Write** — confirm output paths, optionally scan for reference genome and GATK resources, write files

Happy path = **3 interactions**: BAM folder → Accept table → Write confirm.

### From sm-alignment output (flags mode)

When you have analysis-ready BAMs from [sm-alignment](https://github.com/scholl-lab/sm-alignment):

```bash
# Generate samples.tsv from BAM directory
python scripts/generate_config.py \
    --bam-folder results/exomes/bqsr/

# Output:
#   Found 4 BAMs
#   Proposed samples.tsv: 2 tumor_only, 2 tumor_normal (4 rows)
#   Written: config/samples.tsv (4 rows)
```

### Options

```bash
# Preview without writing (dry-run)
python scripts/generate_config.py --bam-folder /path/to/bams --dry-run

# Generate samples.tsv AND config.yaml
python scripts/generate_config.py --bam-folder /path/to/bams --config-template

# Specify caller and scatter mode
python scripts/generate_config.py --bam-folder /path/to/bams --config-template \
    --caller mutect2 --scatter-mode chromosome

# Tumor-only analysis only (no tumor-normal pairs)
python scripts/generate_config.py --bam-folder /path/to/bams --analysis-mode tumor_only

# Germline-only analysis
python scripts/generate_config.py --bam-folder /path/to/bams --analysis-mode germline

# Point to specific reference and GATK resource directories
python scripts/generate_config.py --bam-folder /path/to/bams --config-template \
    --ref-dir /resources/ref/GRCh38 --gatk-resource-dir /resources/gatk_bundle

# Overwrite existing files
python scripts/generate_config.py --bam-folder /path/to/bams --force
```

**Full CLI reference:**

| Flag | Default | Description |
|------|---------|-------------|
| `--bam-folder PATH` | *(interactive)* | BAM directory (triggers non-interactive mode) |
| `--bam-extension EXT` | `.merged.dedup.bqsr.bam` | BAM file extension |
| `--analysis-mode MODE` | `both` | `auto`, `tumor_only`, `tumor_normal`, `both`, `germline` |
| `--caller CALLER` | `mutect2` | `mutect2`, `freebayes`, `all` |
| `--output PATH` | `config/samples.tsv` | Output samples.tsv path |
| `--config-template` | off | Also generate config.yaml |
| `--config-output PATH` | `config/config.yaml` | Config.yaml output path |
| `--output-dir PATH` | *(inferred)* | Pipeline output directory (default: BAM folder parent + `variant_calls`) |
| `--ref-dir PATH` | *(auto-scan)* | Reference genome directory |
| `--gatk-resource-dir PATH` | *(auto-scan)* | GATK resource VCF directory |
| `--scatter-mode MODE` | `chromosome` | `chromosome`, `interval`, `none` |
| `--dry-run` | off | Preview without writing |
| `--force` | off | Overwrite without asking |

When `--config-template` is used, the script scans for reference data:

1. **Explicit** `--ref-dir` and `--gatk-resource-dir` (if provided)
2. **Relative** paths near the project (`resources/ref/GRCh38/`, `analysis/GATK_resource_bundle/`, etc.)
3. **Shared BIH HPC** locations (`/data/cephfs-1/work/groups/scholl/shared/ref/GRCh38/`)

Discovered paths are written into `config/config.yaml`. Paths not found are marked with `EDIT_ME:` placeholders.

### Role detection

The script recognizes tumor and normal BAMs by filename suffix:

| Pattern | Detected role | Examples |
|---------|--------------|----------|
| `-T`, `-T1`, `_TS` | tumor | `SAMPLE01-T`, `SAMPLE-T1` |
| `_tumor`, `.FFPE`, `-met`, `_primary` | tumor | `S1_tumor`, `S1.FFPE` |
| `-L`, `-L1`, `_L2` | tumor (lesion) | `SAMPLE_P1_L1`, `SAMPLE-L` |
| `-N`, `-N1` | normal | `SAMPLE01-N`, `SAMPLE_P1_N1` |
| `_normal`, `-blood`, `_germline`, `.pbmc` | normal | `S1_normal`, `S1-blood` |
| *(no match)* | unknown (prompted) | `SAMPLE_DNA_01_STREAM` |

Pairing uses patient ID extraction (strip role suffix) first, then falls back to string similarity for complex names.

---

## Configuration Reference

### `config/config.yaml`

Review and adjust after generating:

```yaml
caller: "mutect2"                    # "mutect2" | "freebayes" | "all"

ref:
  genome: "/path/to/GRCh38.fna"     # reference FASTA (must have .fai)
  build: "GRCh38"                    # "GRCh37" or "GRCh38"

gatk_resources:                      # required for Mutect2
  panel_of_normals: "/path/to/1000g_pon.hg38.vcf.gz"
  af_only_gnomad: "/path/to/af-only-gnomad.hg38.vcf.gz"
  common_biallelic_gnomad: "/path/to/af-only-gnomad.hg38.common_biallelic.vcf.gz"

paths:
  samples: "config/samples.tsv"      # path to sample sheet
  bam_folder: "results/exomes/bqsr"  # input BAMs
  output_folder: "results/exomes/variant_calls"

bam:
  file_extension: ".merged.dedup.bqsr.bam"

scatter:
  mode: "chromosome"                 # "chromosome" | "interval" | "none"
  count: 400                         # intervals (for interval mode)
  chromosomes: ["chr1", ..., "chrY", "chrM"]

params:                              # extra CLI args passed through to tools
  mutect2:
    extra: "--genotype-germline-sites true --genotype-pon-sites true"
  freebayes:
    extra: "--min-coverage 20 --limit-coverage 500"
  bcftools_norm:
    extra: "-m-any --force -a --atom-overlaps ."
```

| Section | Key settings |
|---------|-------------|
| `caller` | Which caller(s) to run |
| `ref` | Reference genome path and build |
| `gatk_resources` | Panel of normals, gnomAD files (required for Mutect2) |
| `paths` | Input BAM folder, output folder, samples TSV path |
| `bam` | BAM file extension (must match sm-alignment output) |
| `scatter` | Parallelization mode and chromosome list |
| `params` | Tool-specific extra arguments (passthrough strings) |

### `config/samples.tsv`

One row per analysis. Generated by `scripts/generate_config.py` or written manually:

```
sample          tumor_bam          normal_bam        analysis_type
IND001_To       IND001.tumor       .                 tumor_only
IND002_TN       IND002.tumor       IND002.normal     tumor_normal
IND003_G        IND003             .                 germline
```

| Column | Description |
|--------|-------------|
| `sample` | Unique identifier (used in output filenames) |
| `tumor_bam` | BAM basename without extension (not the full path) |
| `normal_bam` | Matched normal BAM basename, or `.` if none |
| `analysis_type` | `tumor_only`, `tumor_normal`, or `germline` |

A single patient can have **multiple rows** (e.g., both `tumor_only` and `tumor_normal`). The config generator's default `both` mode creates this automatically.

### Profiles

Resource and executor settings are separated into layered profiles:

| Profile | Purpose | Applied via |
|---------|---------|-------------|
| `profiles/default/` | Per-rule threads, memory, walltime | `--workflow-profile` |
| `profiles/bih/` | BIH SLURM settings (partition, account) | `--profile` |
| `profiles/charite/` | Charite SLURM settings (partition overrides) | `--profile` |
| `profiles/local/` | Local execution (4 cores, conda) | `--profile` |

To adjust resources without touching workflow code, edit `profiles/default/config.v8+.yaml`.

---

## Supported Callers

| Caller | Analysis types | Output |
|--------|---------------|--------|
| **Mutect2** | `tumor_only`, `tumor_normal` | Filtered somatic VCF per sample |
| **FreeBayes** | `germline` | Merged + normalized germline VCF |

Set `caller: "all"` in config to run both.

---

## Run the Pipeline

### Dry-run first

Always verify the execution plan before submitting:

```bash
conda activate snakemake

snakemake -s workflow/Snakefile --configfile config/config.yaml -n
```

### Submit to SLURM

```bash
# Create log directory
mkdir -p slurm_logs

# Submit — cluster is auto-detected
sbatch scripts/run_snakemake.sh
```

The launcher auto-detects whether you're on BIH HPC or Charite HPC:

| Cluster | Detection | Behavior |
|---------|-----------|----------|
| **BIH HPC** | `/etc/xdg/snakemake/cubi-v1` exists or hostname matches `cubi`/`bihealth` | `--profile profiles/bih` |
| **Charite HPC** | `/etc/profile.d/conda.sh` exists or hostname matches `charite`/`.sc-` | `--profile profiles/charite`, partition overrides for long jobs |
| **Other/local** | Fallback if neither BIH nor Charite is detected | `--profile profiles/local` |

### Usage examples

```bash
# Custom config file
sbatch scripts/run_snakemake.sh config/my_config.yaml

# Custom job name (visible in squeue)
sbatch --job-name=sm_exomes scripts/run_snakemake.sh

# Pass extra Snakemake flags
sbatch scripts/run_snakemake.sh config/config.yaml --forceall

# Override resources for a single run
sbatch scripts/run_snakemake.sh config/config.yaml \
    --set-threads mutect2_call=8 --set-resources mutect2_call:mem_mb=24000

# Run locally without cluster submission
bash scripts/run_snakemake.sh config/config.yaml
```

### Monitor progress

```bash
# Check running jobs
squeue -u $USER

# Follow coordinator log
tail -f slurm_logs/sm_calling-*.log

# Follow individual rule logs
tail -f slurm_logs/slurm-*.log
```

### Invoke Snakemake directly

```bash
snakemake \
    -s workflow/Snakefile \
    --configfile config/config.yaml \
    --workflow-profile profiles/default \
    --profile profiles/bih
```

---

## Pipeline Architecture

```
workflow/Snakefile
    │
    ├── rules/common.smk      Config shortcuts, samples_df, input helpers
    ├── rules/helpers.py       Pure Python functions (unit-testable)
    ├── rules/scatter.smk      Interval scattering (shared by both callers)
    ├── rules/mutect2.smk      Mutect2: call → gather → contamination → filter
    └── rules/freebayes.smk    FreeBayes: call → merge + normalize
```

### Mutect2 DAG

```
BAM files
    ├── mutect2_call (per sample x scatter_unit)        [temp]
    ├── get_pileup_summaries (per unique BAM)
    ├── gather_mutect2_vcfs (per sample)
    ├── merge_mutect2_stats (per sample)
    ├── learn_read_orientation (per sample)
    ├── calculate_contamination (per sample)
    └── filter_mutect_calls (per sample)                [protected]
```

### FreeBayes DAG

```
BAM files
    ├── freebayes_call (per scatter_unit)               [temp]
    └── merge_freebayes_vcfs                            [protected]
```

---

## Tools & Conda Environments

| Conda env | Tools | Used by |
|-----------|-------|---------|
| `workflow/envs/gatk.yaml` | gatk4 4.6.1.0, samtools 1.21 | Mutect2 pipeline |
| `workflow/envs/freebayes.yaml` | freebayes 1.3.8, htslib 1.21 | FreeBayes calling |
| `workflow/envs/bcftools.yaml` | bcftools 1.21, htslib 1.21 | VCF merge + normalize |

Conda environments are created automatically by Snakemake on first run (`software-deployment-method: conda` in the workflow profile).

To skip per-rule conda and use pre-installed tools instead:

```bash
# Install tools into the snakemake environment
mamba install -n snakemake -c bioconda -c conda-forge \
    gatk4=4.6.1.0 samtools=1.21 freebayes=1.3.8 bcftools=1.21 htslib=1.21

# Run with --sdm none
sbatch scripts/run_snakemake.sh config/config.yaml --sdm none
```

---

## Development

```bash
# Install dev tools
pip install ruff mypy snakefmt shellcheck-py pytest

# Run all checks
make lint

# Auto-format
make format

# Run unit tests
make test-unit

# Run all tests (requires snakemake + conda)
make test

# Type check
make typecheck
```

---

## Project Structure

```
sm-calling/
├── config/
│   ├── config.yaml              Main configuration
│   ├── samples.tsv              Sample sheet
│   └── README.md                Config field reference
├── workflow/
│   ├── Snakefile                Entry point (Snakemake 8+)
│   ├── rules/                   Rule files + helpers.py
│   ├── envs/                    Conda environment YAMLs
│   └── schemas/                 JSON schemas for validation
├── profiles/
│   ├── default/                 Per-rule resources (threads, memory, walltime)
│   ├── bih/                     BIH SLURM executor settings
│   ├── charite/                 Charite SLURM executor settings
│   └── local/                   Local execution (no cluster)
├── scripts/
│   ├── generate_config.py       Config generator (interactive + CLI)
│   └── run_snakemake.sh         SLURM launcher (auto-detects cluster)
├── tests/                       pytest unit and integration tests
└── deprecated/                  Old standalone workflows (for reference)
```

## License

MIT License. See [LICENSE](LICENSE) for details.
