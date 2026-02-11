# sm-calling Refactoring Plan

## Executive Summary

Refactor the `sm-calling` variant calling Snakemake workflow to match the modern
architecture established in `sm-alignment`, adopting Snakemake 8+ best practices,
the Snakemake Workflow Catalog standard layout, and production-grade
infrastructure (schema validation, profiles, testing, CI/CD).

**Current state:** 7 standalone `.smk` files in `scripts/snakemake/`, 7 SLURM
launcher scripts, 3 config files, no tests, no schema validation, no conda
environment definitions, extensive code duplication.

**Target state:** Single unified `workflow/Snakefile` targeting **Snakemake 8+**
with modular rule includes, validated configuration, reproducible conda
environments, SLURM executor plugin-based profiles, unit + integration tests,
and a single universal launcher. Adopts modern Snakemake 8 idioms throughout:
`min_version("8.0")`, `localrule: True`, `ensure()`, `lookup()`/`collect()`,
`software-deployment-method`, and the `snakemake-executor-plugin-slurm`.

---

## Table of Contents

1. [Current Architecture Analysis](#1-current-architecture-analysis)
2. [Snakemake 8+ Modernization](#2-snakemake-8-modernization)
3. [Target Architecture](#3-target-architecture)
4. [Phase 1: Directory Structure & Scaffold](#phase-1-directory-structure--scaffold)
5. [Phase 2: Common Infrastructure](#phase-2-common-infrastructure)
6. [Phase 3: Mutect2 Pipeline Rules](#phase-3-mutect2-pipeline-rules)
7. [Phase 4: FreeBayes Pipeline Rules](#phase-4-freebayes-pipeline-rules)
8. [Phase 5: Profiles & Resource Management (Snakemake 8 Executor Plugin)](#phase-5-profiles--resource-management-snakemake-8-executor-plugin)
9. [Phase 6: Launcher & HPC Integration](#phase-6-launcher--hpc-integration)
10. [Phase 7: Testing & CI/CD](#phase-7-testing--cicd)
11. [Phase 8: Documentation & Cleanup](#phase-8-documentation--cleanup)
12. [Migration Guide](#migration-guide)
13. [Risk Assessment](#risk-assessment)

---

## 1. Current Architecture Analysis

### 1.1 Problems Identified

#### Structural

| Issue | Impact | Severity |
|-------|--------|----------|
| Non-standard layout (`scripts/snakemake/`) | Cannot join Workflow Catalog; confuses new users | High |
| No main `Snakefile` entry point | Each pipeline is a standalone invocation | High |
| 7 SLURM launcher scripts | Maintenance burden; inconsistent parameters | Medium |
| Config split between root `config.yaml` and `config/` | Confusing; some files override others | Medium |

#### Code Quality

| Issue | Impact | Severity |
|-------|--------|----------|
| Metadata loading duplicated in every `.smk` file | DRY violation; divergence risk | High |
| `get_mem_from_threads()` duplicated 5 times | Inconsistency risk | Medium |
| Named conda envs (`conda: "gatk"`) instead of YAML defs | Not reproducible across machines | Critical |
| No schema validation for config or metadata | Silent misconfiguration | High |
| No wildcard constraints | Ambiguous rule matching risk | Medium |
| Inline lambdas for input functions | Hard to test, read, and reuse | Medium |
| Manual `os.makedirs()` at module load | Unnecessary; Snakemake creates output dirs | Low |
| No `temp()` on Mutect2 scattered outputs | Wastes disk (25 VCFs per sample persist) | Medium |
| No `retries` on Mutect2 rules | Jobs fail permanently on transient errors | Medium |
| No `benchmark` directives | Cannot profile resource usage | Low |

#### Missing Infrastructure

| Item | Status in sm-alignment | Status in sm-calling |
|------|------------------------|----------------------|
| Schema validation | Yes (config + samples) | None |
| Conda env YAML files | Yes (3 pinned files) | None (named envs) |
| Snakemake profiles | Yes (default + charite) | None (relies on external) |
| Unit tests | Yes (99 tests) | None |
| CI/CD | Makefile targets | None |
| `helpers.py` (testable) | Yes | None |
| `common.smk` | Yes | None |
| Universal launcher | Yes (cluster auto-detect) | 7 separate scripts |
| `generate_config.py` | Yes (interactive wizard) | None |

### 1.2 Redundant Workflows

The repository contains two parallel approaches for Mutect2:

- **Modular sequential:** `mutect2_calling.smk` -> `merge_mutect2_calls.smk` ->
  `gatk_calculate_contamination.smk` -> `gatk_filter_mutect_calls.smk`
  (4 separate SLURM submissions, manual orchestration)
- **End-to-end:** `mutect2_pipeline.smk` (single workflow with all steps)

The end-to-end approach is superior and should be the only supported path.
The modular sequential approach should be deprecated.

Similarly, `bcftools_concat.smk` overlaps with `merge_mutect2_calls.smk` and
the merge rules in `mutect2_pipeline.smk`.

---

## 2. Snakemake 8+ Modernization

This refactoring targets **Snakemake 8+** (current stable: 9.x). Snakemake 8
introduced a plugin-based architecture that fundamentally changes how workflows
interact with clusters, software environments, and storage. Every aspect of
this plan is written for the modern API.

### 2.1 Breaking Changes from Snakemake 7 → 8

| Old (Snakemake 7) | New (Snakemake 8+) | Impact |
|--------------------|---------------------|--------|
| `--use-conda` | `--software-deployment-method conda` (or `--sdm conda`) | Profile + launcher |
| `--cluster "sbatch ..."` | `executor: slurm` via `snakemake-executor-plugin-slurm` | Profile overhaul |
| `subworkflow:` directive | `module:` directive | Not used here |
| `dynamic()` | `checkpoint` | Not used here |
| `localrules: all, ...` (global) | `localrule: True` (per-rule, preferred) | Rule syntax |
| No output validation | `ensure(path, non_empty=True)` | Output safety |
| Lambda input functions | `lookup()` / `collect()` (declarative) | Cleaner rules |
| No version enforcement | `min_version("8.0")` | Snakefile header |

### 2.2 Key Snakemake 8+ Features to Adopt

#### `min_version("8.0")`
Enforce minimum Snakemake version at the top of the Snakefile:
```python
from snakemake.utils import min_version, validate
min_version("8.0")
```

#### SLURM Executor Plugin
Replace the old `cluster:` string with the dedicated executor plugin:
```bash
pip install snakemake-executor-plugin-slurm
```

Profile uses `executor: slurm` and SLURM-specific resources:
```yaml
executor: slurm
default-resources:
  slurm_partition: "default"
  slurm_account: "your_account"
  mem_mb_per_cpu: 2000
  runtime: 60
```

Available SLURM resources: `slurm_partition`, `slurm_account`, `slurm_extra`,
`mem_mb`, `mem_mb_per_cpu`, `runtime`, `cpus_per_task`, `tasks`, `nodes`,
`constraint`.

#### `software-deployment-method`
In profiles (replaces `--use-conda`):
```yaml
software-deployment-method:
  - conda
```

#### `localrule: True` (per-rule)
Modern replacement for global `localrules:` directive. Keeps the local
designation co-located with the rule definition:
```python
rule all:
    localrule: True
    input: get_final_outputs()
```

#### `lookup()` and `collect()`
Declarative input resolution replacing lambda-based metadata lookups:
```python
# Old: lambda-based (hard to test, verbose)
rule mutect2_call:
    input:
        tumor_bam=lambda wc: os.path.join(
            BAM_FOLDER,
            metadata_dict[wc.analysis_key]["bam1_file_basename"] + BAM_EXT
        ),

# New: lookup() (declarative, auto-substitutes wildcards)
rule mutect2_call:
    input:
        tumor_bam=lookup(
            query="sample == '{sample}'",
            cols="tumor_bam",
            within=samples_df,
        ),
```

And `collect()` for target aggregation in `rule all`:
```python
rule all:
    localrule: True
    input:
        collect(
            os.path.join(MUTECT2_FILTERED_DIR, "{item.sample}.filtered.vcf.gz"),
            item=lookup(
                query="analysis_type in ['tumor_only', 'tumor_normal']",
                within=samples_df,
            ),
        ),
```

**When to keep lambdas:** Complex conditional logic (e.g., returning `[]` when
no matched normal exists) still requires lambdas or named functions. Use
`lookup()` for simple "get column value by wildcard" patterns.

#### `ensure(non_empty=True)`
Validate output file properties -- catches truncated/empty outputs:
```python
rule filter_mutect_calls:
    output:
        vcf=ensure(
            protected(os.path.join(FILTERED_DIR, "{sample}.filtered.vcf.gz")),
            non_empty=True,
        ),
```

#### `--shared-fs-usage`
Granular control over shared filesystem assumptions on clusters:
```yaml
shared-fs-usage:
  - persistence
  - software-deployment
  - sources
  - source-cache
  - input-output
```

#### `--rerun-triggers`
Control what triggers job re-execution:
```yaml
# Production: all triggers (default -- full reproducibility)
rerun-triggers:
  - mtime
  - code
  - input
  - params
  - software-env

# Development: skip code trigger to avoid reruns on shell tweaks
rerun-triggers:
  - mtime
  - input
  - params
  - software-env
```

#### Attempt-based Resource Scaling with `retries:`
The `attempt` parameter in resource callables (starts at 1, increments per
retry) enables automatic resource escalation:
```python
rule mutect2_call:
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 17600 * attempt,
        runtime=lambda wildcards, attempt: 4320 * attempt,
```

### 2.3 sm-alignment Snakemake 8 Gap Analysis

The recently refactored `sm-alignment` adopted some modern patterns but still
uses **Snakemake 7-era** constructs in several places. This refactoring of
`sm-calling` should go further:

| Feature | sm-alignment Status | sm-calling Target |
|---------|--------------------|--------------------|
| `min_version()` | Not used | `min_version("8.0")` |
| Executor plugin | Uses `cluster:` string in Charite profile | `executor: slurm` everywhere |
| `--sdm` | Uses `software-deployment-method` in profile | Same |
| `localrule: True` | Uses global `localrules:` | Per-rule `localrule: True` |
| `lookup()` / `collect()` | Not used (lambda functions) | Use where appropriate |
| `ensure()` | Not used | On critical VCF outputs |
| `--shared-fs-usage` | Not configured | Configured in profiles |
| `--rerun-triggers` | Not configured | Documented in profiles |
| Attempt-based resources | Not used | On variant calling rules |

---

## 3. Target Architecture

### 3.1 Directory Layout (Snakemake Workflow Catalog Standard)

```
sm-calling/
├── .github/
│   └── workflows/
│       └── main.yaml                  # CI: lint, format, dry-run, test
├── .snakemake-workflow-catalog.yml    # Catalog metadata
├── config/
│   ├── config.yaml                    # Unified hierarchical configuration
│   ├── samples.tsv                    # Sample metadata (replaces calling_metadata.tsv)
│   └── README.md                      # Configuration reference
├── workflow/
│   ├── Snakefile                      # Entry point: validate, include, rule all
│   ├── rules/
│   │   ├── common.smk                 # Config shortcuts, metadata, helper wrappers
│   │   ├── scatter.smk                # Interval scattering (shared by all callers)
│   │   ├── mutect2.smk                # Mutect2 calling + merge + contamination + filter
│   │   ├── freebayes.smk              # FreeBayes calling + merge + normalize
│   │   └── helpers.py                 # Pure Python functions (testable)
│   ├── envs/
│   │   ├── gatk.yaml                  # gatk4, samtools, tabix (pinned)
│   │   ├── bcftools.yaml              # bcftools (pinned)
│   │   └── freebayes.yaml             # freebayes, htslib (pinned)
│   └── schemas/
│       ├── config.schema.yaml         # JSON Schema for config.yaml
│       └── samples.schema.yaml        # JSON Schema for samples.tsv
├── profiles/
│   ├── default/
│   │   └── config.yaml                # Default threads/mem/runtime per rule
│   └── charite/
│       └── config.yaml                # Charite HPC SLURM settings
├── scripts/
│   ├── run_snakemake.sh               # Universal SLURM launcher (from sm-alignment)
│   └── generate_config.py             # Interactive config + metadata generator
├── tests/
│   ├── conftest.py                    # pytest fixtures
│   ├── test_helpers.py                # Unit tests for helpers.py
│   ├── test_schema_validation.py      # Config/samples schema tests
│   └── test_dryrun.py                 # Integration dry-run tests
├── deprecated/                        # Old standalone workflows (removal target)
│   ├── scripts/                       # Old launchers + .smk files
│   └── README.md                      # Migration guide
├── Makefile                           # lint, format, test targets
├── pyproject.toml                     # ruff, mypy, snakefmt, pytest config
├── CLAUDE.md                          # AI assistant guidance
├── README.md                          # User documentation
└── LICENSE
```

### 3.2 Unified Configuration Design

Replace the current 3 config files with a single hierarchical `config/config.yaml`:

```yaml
# --- Caller selection ---
caller: "mutect2"  # "mutect2" | "freebayes" | "all"

# --- Reference genome ---
ref:
  genome: "path/to/GRCh38.fna"
  build: "GRCh38"

# --- GATK resources (for Mutect2) ---
gatk_resources:
  panel_of_normals: "path/to/1000g_pon.hg38.vcf.gz"
  af_only_gnomad: "path/to/af-only-gnomad.hg38.vcf.gz"
  common_biallelic_gnomad: "path/to/af-only-gnomad.hg38.common_biallelic.vcf.gz"

# --- Paths ---
paths:
  samples: "config/samples.tsv"
  bam_folder: "results/exomes/bqsr"
  output_folder: "results/exomes/variant_calls"
  log_subdir: "logs"

# --- BAM file settings ---
bam:
  file_extension: ".merged.dedup.bqsr.bam"

# --- Scatter settings ---
scatter:
  mode: "chromosome"       # "chromosome" | "interval" | "none"
  count: 400               # For interval mode
  chromosomes:             # For chromosome mode
    - "chr1"
    - "chr2"
    # ... chr3-chr22
    - "chrX"
    - "chrY"
    - "chrM"

# --- Tool-specific parameters (passthrough) ---
params:
  mutect2:
    extra: "--genotype-germline-sites true --genotype-pon-sites true"
  freebayes:
    extra: "--min-coverage 20 --limit-coverage 500 --use-best-n-alleles 4"
  bcftools_norm:
    extra: "-m-any --force -a --atom-overlaps ."
```

### 3.3 Unified Sample Sheet Design

Replace `calling_metadata.tsv` with a cleaner `config/samples.tsv`:

```
sample          tumor_bam               normal_bam              analysis_type
IND001_To       IND001.tumor            .                       tumor_only
IND002_TN       IND002.tumor            IND002.normal           tumor_normal
IND003_G        IND003                  .                       germline
```

Key changes:
- `sample` column replaces `individual1_analysis` composite key
- `tumor_bam` / `normal_bam` replace `bam1_file_basename` / `bam2_file_basename`
- `analysis_type` enum replaces implicit detection from empty fields
- `.` or empty for optional normal (explicit, not implicit)
- FreeBayes germline samples use a BAM list derived from `analysis_type == "germline"`

### 3.4 Rule Architecture

```
workflow/Snakefile
    │
    ├── min_version("8.0")
    ├── configfile + validate()
    ├── samples = pd.read_table() + validate()
    │
    ├── include: rules/common.smk      # Config shortcuts, samples_df, helpers
    ├── include: rules/scatter.smk     # scatter_intervals_by_ns, split_intervals
    ├── include: rules/mutect2.smk     # Full Mutect2 pipeline (if caller includes mutect2)
    └── include: rules/freebayes.smk   # Full FreeBayes pipeline (if caller includes freebayes)

    rule all:
        localrule: True
        input: get_final_outputs()     # Dispatches based on config["caller"]
```

**DAG (Mutect2 path):**
```
BAM files (from bam_folder)
    ├─→ [scatter_intervals if mode=interval]
    ├─→ mutect2_call (per sample x scatter_unit)        → temp()
    ├─→ get_pileup_summaries (per unique BAM)
    ├─→ gather_vcfs (per sample)
    ├─→ merge_stats (per sample)
    ├─→ learn_read_orientation (per sample)
    ├─→ calculate_contamination (per sample)
    └─→ filter_mutect_calls (per sample)                 → protected()
```

**DAG (FreeBayes path):**
```
BAM files (from bam_folder or bam_list)
    ├─→ [scatter_intervals if mode=interval]
    ├─→ freebayes_call (per scatter_unit)                → temp()
    └─→ merge_normalize_vcfs                             → protected()
```

---

## Phase 1: Directory Structure & Scaffold

**Goal:** Create the standard layout, move old files to `deprecated/`.

### Tasks

1.1. Create standard directory tree:
```
workflow/
workflow/rules/
workflow/envs/
workflow/schemas/
profiles/default/
profiles/charite/
scripts/
tests/
deprecated/
```

1.2. Move existing files to `deprecated/`:
```
scripts/snakemake/*.smk           → deprecated/scripts/snakemake/
scripts/run_*.sh                  → deprecated/scripts/launchers/
scripts/submit_*.sh               → deprecated/scripts/launchers/
config.yaml (root)                → deprecated/configs/config.yaml
config/config_mutect2.yaml        → deprecated/configs/config_mutect2.yaml
config/config_freebayes_germline.yaml → deprecated/configs/config_freebayes_germline.yaml
```

1.3. Create empty scaffold files:
- `workflow/Snakefile`
- `workflow/rules/common.smk`
- `workflow/rules/scatter.smk`
- `workflow/rules/mutect2.smk`
- `workflow/rules/freebayes.smk`
- `workflow/rules/helpers.py`
- `config/config.yaml` (new unified)
- `config/samples.tsv` (new unified)

1.4. Create `deprecated/README.md` with migration notes.

1.5. Update `.gitignore` to include:
```
results/
.snakemake/
slurm_logs/
__pycache__/
*.pyc
```

### Acceptance Criteria
- All old files preserved in `deprecated/`
- Empty scaffold compiles (`snakemake -n` with stub rule all)
- Git history preserved (use `git mv` where possible)

---

## Phase 2: Common Infrastructure

**Goal:** Build the shared foundation: config, schemas, helpers, conda envs.

### Tasks

2.1. **Write `config/config.yaml`** (unified hierarchical config as defined in
section 2.2 above). Merge settings from all 3 existing config files. Use
the Mutect2 pipeline config as the primary source, add FreeBayes params.

2.2. **Write `config/samples.tsv`** with the new column schema (section 2.3).
Provide example rows covering tumor-only, tumor-normal, and germline.

2.3. **Write `workflow/schemas/config.schema.yaml`**:
```yaml
$schema: "https://json-schema.org/draft/2020-12/schema"
description: Configuration for sm-calling variant calling pipeline
properties:
  caller:
    type: string
    enum: ["mutect2", "freebayes", "all"]
  ref:
    type: object
    properties:
      genome: { type: string }
      build: { type: string, enum: ["GRCh37", "GRCh38"] }
    required: [genome, build]
  gatk_resources:
    type: object
    properties:
      panel_of_normals: { type: string }
      af_only_gnomad: { type: string }
      common_biallelic_gnomad: { type: string }
    required: [panel_of_normals, af_only_gnomad, common_biallelic_gnomad]
  paths:
    type: object
    properties:
      samples: { type: string }
      bam_folder: { type: string }
      output_folder: { type: string }
      log_subdir: { type: string, default: "logs" }
    required: [samples, bam_folder, output_folder]
  bam:
    type: object
    properties:
      file_extension: { type: string, default: ".merged.dedup.bqsr.bam" }
  scatter:
    type: object
    properties:
      mode: { type: string, enum: ["chromosome", "interval", "none"] }
      count: { type: integer, minimum: 1, default: 400 }
      chromosomes: { type: array, items: { type: string } }
    required: [mode]
  params:
    type: object
    properties:
      mutect2: { type: object }
      freebayes: { type: object }
      bcftools_norm: { type: object }
required: [caller, ref, paths, bam, scatter]
```

2.4. **Write `workflow/schemas/samples.schema.yaml`**:
```yaml
$schema: "https://json-schema.org/draft/2020-12/schema"
description: An entry in the calling sample sheet
properties:
  sample:
    type: string
    description: Unique sample/analysis identifier
  tumor_bam:
    type: string
    description: Basename of tumor/primary BAM file
  normal_bam:
    type: string
    description: Basename of matched normal BAM (empty or '.' if none)
  analysis_type:
    type: string
    enum: ["tumor_only", "tumor_normal", "germline"]
required: [sample, tumor_bam, analysis_type]
```

2.5. **Write `workflow/envs/gatk.yaml`**:
```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - gatk4=4.6.1.0
  - samtools=1.21
```

2.6. **Write `workflow/envs/bcftools.yaml`**:
```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - bcftools=1.21
  - htslib=1.21
```

2.7. **Write `workflow/envs/freebayes.yaml`**:
```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - freebayes=1.3.8
  - htslib=1.21
```

2.8. **Write `workflow/rules/helpers.py`** -- pure Python, no Snakemake imports:
```python
"""Pure helper functions for the sm-calling workflow.

All functions are free of Snakemake imports so they can be unit-tested.
"""
import os
import pandas as pd


def get_java_opts(mem_mb: int, tmpdir: str) -> str:
    """Derive GATK --java-options from allocated resources."""
    xmx = int(mem_mb * 0.8)
    xms = int(mem_mb * 0.2)
    return f"-Xms{xms}m -Xmx{xmx}m -Djava.io.tmpdir={tmpdir}"


def get_samples(samples_df: pd.DataFrame) -> list[str]:
    """Return sorted list of unique sample identifiers."""
    return sorted(samples_df["sample"].unique().tolist())


def get_mutect2_samples(samples_df: pd.DataFrame) -> list[str]:
    """Return samples requiring Mutect2 calling (tumor_only or tumor_normal)."""
    mask = samples_df["analysis_type"].isin(["tumor_only", "tumor_normal"])
    return sorted(samples_df.loc[mask, "sample"].unique().tolist())


def get_germline_samples(samples_df: pd.DataFrame) -> list[str]:
    """Return samples requiring germline calling."""
    mask = samples_df["analysis_type"] == "germline"
    return sorted(samples_df.loc[mask, "sample"].unique().tolist())


def get_unique_bam_basenames(samples_df: pd.DataFrame) -> list[str]:
    """Return all unique BAM basenames across tumor and normal columns."""
    bams = set()
    for _, row in samples_df.iterrows():
        if row["tumor_bam"] and row["tumor_bam"] != ".":
            bams.add(row["tumor_bam"])
        if pd.notna(row.get("normal_bam")) and row["normal_bam"] not in ("", "."):
            bams.add(row["normal_bam"])
    return sorted(bams)


def has_matched_normal(samples_df: pd.DataFrame, sample: str) -> bool:
    """Check whether a sample has a matched normal BAM."""
    row = samples_df.loc[samples_df["sample"] == sample].iloc[0]
    normal = row.get("normal_bam", "")
    return bool(normal) and normal != "."


def get_scatter_units(
    mode: str,
    chromosomes: list[str] | None = None,
    scatter_count: int = 400,
) -> list[str]:
    """Return the list of scatter units based on mode."""
    if mode == "chromosome":
        return chromosomes or []
    elif mode == "interval":
        return [f"{i:04d}-scattered" for i in range(scatter_count)]
    else:
        return ["all"]


def build_freebayes_params(params_dict: dict) -> str:
    """Convert a dict of FreeBayes params to a CLI string."""
    parts = []
    for key, value in params_dict.items():
        flag = f"--{key}" if not key.startswith("-") else key
        if value is True or value == "true":
            parts.append(flag)
        elif value is not False and value != "false":
            parts.append(f"{flag} {value}")
    return " ".join(parts)
```

2.9. **Write `workflow/Snakefile`** (Snakemake 8+ entry point):
```python
"""sm-calling: Somatic and germline variant calling pipeline.

Requires Snakemake 8.0+. Uses the SLURM executor plugin for cluster execution.
"""
from snakemake.utils import min_version, validate

min_version("8.0")

import pandas as pd

# ── Configuration ─────────────────────────────────────────────────────
configfile: "config/config.yaml"

validate(config, "schemas/config.schema.yaml")

samples_df = pd.read_table(config["paths"]["samples"]).set_index(
    "sample", drop=False
)
validate(samples_df, "schemas/samples.schema.yaml")

# ── Includes ──────────────────────────────────────────────────────────
include: "rules/common.smk"
include: "rules/scatter.smk"
include: "rules/mutect2.smk"
include: "rules/freebayes.smk"

# ── Target rule ───────────────────────────────────────────────────────
rule all:
    localrule: True
    input:
        get_final_outputs(),
```

2.10. **Write `workflow/rules/common.smk`** (Snakemake 8+ idioms):
```python
"""Shared configuration, metadata loading, and helper wrappers."""
import os
import pandas as pd
from rules.helpers import (
    get_java_opts as _get_java_opts_impl,
    get_samples as _get_samples_impl,
    get_mutect2_samples as _get_mutect2_samples_impl,
    get_germline_samples as _get_germline_samples_impl,
    get_unique_bam_basenames as _get_unique_bam_basenames_impl,
    has_matched_normal as _has_matched_normal_impl,
    get_scatter_units as _get_scatter_units_impl,
    build_freebayes_params as _build_freebayes_params_impl,
)

# ── Config shortcuts ──────────────────────────────────────────────────
CALLER = config["caller"]
REF = config["ref"]["genome"]
REF_BUILD = config["ref"]["build"]

BAM_FOLDER = config["paths"]["bam_folder"]
OUTPUT_DIR = config["paths"]["output_folder"]
LOG_SUBDIR = config["paths"].get("log_subdir", "logs")
LOG_DIR = os.path.join(OUTPUT_DIR, LOG_SUBDIR)

BAM_EXT = config["bam"]["file_extension"]

SCATTER_MODE = config["scatter"]["mode"]
SCATTER_COUNT = config["scatter"].get("count", 400)
CHROMOSOMES = config["scatter"].get("chromosomes", [])

# GATK resources (optional for FreeBayes-only runs)
GATK_RES = config.get("gatk_resources", {})
PON = GATK_RES.get("panel_of_normals", "")
AF_GNOMAD = GATK_RES.get("af_only_gnomad", "")
COMMON_GNOMAD = GATK_RES.get("common_biallelic_gnomad", "")

# Tool params
MUTECT2_EXTRA = config.get("params", {}).get("mutect2", {}).get("extra", "")
FREEBAYES_PARAMS = config.get("params", {}).get("freebayes", {})
BCFTOOLS_NORM_EXTRA = config.get("params", {}).get("bcftools_norm", {}).get("extra", "")

# ── Output subdirectories ─────────────────────────────────────────────
MUTECT2_DIR = os.path.join(OUTPUT_DIR, "mutect2")
MUTECT2_SCATTER_DIR = os.path.join(MUTECT2_DIR, "scattered")
MUTECT2_MERGED_DIR = os.path.join(MUTECT2_DIR, "merged")
MUTECT2_CONTAM_DIR = os.path.join(MUTECT2_DIR, "contamination")
MUTECT2_FILTERED_DIR = os.path.join(MUTECT2_DIR, "filtered")

FREEBAYES_DIR = os.path.join(OUTPUT_DIR, "freebayes")
FREEBAYES_SCATTER_DIR = os.path.join(FREEBAYES_DIR, "scattered")

INTERVALS_DIR = os.path.join(OUTPUT_DIR, "intervals")

# ── Scatter units ─────────────────────────────────────────────────────
SCATTER_UNITS = _get_scatter_units_impl(SCATTER_MODE, CHROMOSOMES, SCATTER_COUNT)

# ── Samples metadata (already loaded and validated in Snakefile) ──────
# samples_df is available from the Snakefile scope

def get_all_samples():
    return _get_samples_impl(samples_df)

def get_mutect2_samples():
    return _get_mutect2_samples_impl(samples_df)

def get_germline_samples():
    return _get_germline_samples_impl(samples_df)

def get_unique_bams():
    return _get_unique_bam_basenames_impl(samples_df)

def has_normal(sample):
    return _has_matched_normal_impl(samples_df, sample)

# ── Input helper functions ────────────────────────────────────────────
# Use lookup() for simple column retrieval (Snakemake 8+ idiom).
# Keep named functions for complex conditional logic (e.g., optional normals).

def get_tumor_bam(wildcards):
    """Return full path to tumor BAM for a sample."""
    row = samples_df.loc[wildcards.sample]
    return os.path.join(BAM_FOLDER, row["tumor_bam"] + BAM_EXT)

def get_normal_bam(wildcards):
    """Return list with normal BAM path, or empty list if none.

    This cannot use lookup() because it needs conditional empty-list logic.
    """
    row = samples_df.loc[wildcards.sample]
    normal = row.get("normal_bam", "")
    if normal and normal != "." and pd.notna(normal):
        return [os.path.join(BAM_FOLDER, normal + BAM_EXT)]
    return []

def get_normal_name(wildcards):
    """Return '-normal SAMPLE_NAME' arg or empty string."""
    row = samples_df.loc[wildcards.sample]
    normal = row.get("normal_bam", "")
    if normal and normal != "." and pd.notna(normal):
        return f"-normal {normal}"
    return ""

def get_interval_arg(wildcards):
    """Return -L/--targets argument based on scatter mode and unit."""
    unit = wildcards.scatter_unit
    if SCATTER_MODE == "chromosome":
        return f"-L {unit}"
    elif SCATTER_MODE == "interval":
        return f"-L {INTERVALS_DIR}/{unit}.interval_list"
    return ""

def get_java_opts(wildcards, resources):
    """Derive GATK --java-options from allocated resources."""
    return _get_java_opts_impl(resources.mem_mb, resources.tmpdir)

# ── Final output dispatcher ──────────────────────────────────────────
def get_final_outputs():
    """Return all expected final output files based on caller selection.

    Uses collect() + lookup() (Snakemake 8+ idiom) for Mutect2 targets,
    and a static path for the FreeBayes merged output.
    """
    outputs = []
    if CALLER in ("mutect2", "all"):
        outputs.extend(
            collect(
                os.path.join(MUTECT2_FILTERED_DIR, "{item.sample}.filtered.vcf.gz"),
                item=lookup(
                    query="analysis_type in ['tumor_only', 'tumor_normal']",
                    within=samples_df,
                ),
            )
        )
    if CALLER in ("freebayes", "all"):
        freebayes_out = config.get("paths", {}).get(
            "freebayes_merged_vcf",
            os.path.join(FREEBAYES_DIR, "final_merged.vcf.gz"),
        )
        outputs.append(freebayes_out)
    return outputs
```

### Acceptance Criteria
- `snakemake --lint` passes on the Snakefile
- `min_version("8.0")` enforced in Snakefile header
- `validate()` called for both config and samples at startup
- `rule all` uses `localrule: True` (Snakemake 8+ per-rule syntax)
- `get_final_outputs()` uses `collect()` + `lookup()` for Mutect2 targets
- Schema validation catches missing required fields
- `helpers.py` is importable and testable outside Snakemake
- Conda env YAMLs parse correctly with pinned versions

---

## Phase 3: Mutect2 Pipeline Rules

**Goal:** Consolidate all Mutect2 logic into `workflow/rules/mutect2.smk` +
`workflow/rules/scatter.smk`, following GATK best practices.

### Tasks

3.1. **Write `workflow/rules/scatter.smk`** -- shared interval scattering:
- `rule scatter_intervals_by_ns` (conditional on `SCATTER_MODE == "interval"`)
- `rule split_intervals` (conditional on `SCATTER_MODE == "interval"`)
- Both rules output to `INTERVALS_DIR`
- Use `temp()` for intermediate interval files
- Add `benchmark:` directives

3.2. **Write `workflow/rules/mutect2.smk`** with these rules:

| Rule | Source | Key Changes |
|------|--------|-------------|
| `mutect2_call` | `call_variants` from `mutect2_pipeline.smk` | Named input fns, `temp()` outputs, `retries: 2`, `benchmark:`, attempt-based memory |
| `gather_mutect2_vcfs` | `merge_vcfs` from `mutect2_pipeline.smk` | Renamed for clarity |
| `merge_mutect2_stats` | `merge_stats` | Unchanged logic |
| `learn_read_orientation` | `merge_f1r2` | Renamed for clarity |
| `get_pileup_summaries` | `get_pileup_summaries` | Named input fns |
| `calculate_contamination` | `calculate_contamination` | Named input fns, `unpack()` |
| `filter_mutect_calls` | `filter_mutect_calls` | `protected()` output |

3.3. **Key improvements per rule (Snakemake 8+ idioms):**

**`mutect2_call`** -- with `retries:`, attempt-based resource scaling, `ensure()`:
```python
rule mutect2_call:
    """Scatter-parallel Mutect2 somatic variant calling."""
    input:
        tumor_bam=get_tumor_bam,
        normal_bam=get_normal_bam,       # Returns [] if tumor-only
        reference=REF,
        pon=PON,
        gnomad=AF_GNOMAD,
    output:
        vcf=temp(
            ensure(
                os.path.join(MUTECT2_SCATTER_DIR, "{sample}_{scatter_unit}.vcf.gz"),
                non_empty=True,
            )
        ),
        stats=temp(os.path.join(MUTECT2_SCATTER_DIR, "{sample}_{scatter_unit}.vcf.gz.stats")),
        f1r2=temp(os.path.join(MUTECT2_SCATTER_DIR, "{sample}_{scatter_unit}.f1r2.tar.gz")),
    params:
        normal_arg=get_normal_name,
        interval_arg=get_interval_arg,
        java_opts=get_java_opts,
        extra=MUTECT2_EXTRA,
    retries: 2
    resources:
        # Attempt-based scaling (Snakemake 8+): doubles on each retry
        mem_mb=lambda wildcards, attempt: 17600 * attempt,
        runtime=lambda wildcards, attempt: 4320 * attempt,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "mutect2_call/{sample}_{scatter_unit}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/mutect2_call/{sample}_{scatter_unit}.tsv")
    shell:
        r"""
        gatk --java-options '{params.java_opts}' Mutect2 \
            -R {input.reference} \
            -I {input.tumor_bam} \
            {params.normal_arg} \
            --germline-resource {input.gnomad} \
            --panel-of-normals {input.pon} \
            --f1r2-tar-gz {output.f1r2} \
            {params.extra} \
            {params.interval_arg} \
            -O {output.vcf} \
            2> {log}
        """
```

**`filter_mutect_calls`** -- with `ensure()` + `protected()` on final output:
```python
rule filter_mutect_calls:
    """Apply Mutect2 filters using contamination, segmentation, and orientation data."""
    input:
        vcf=os.path.join(MUTECT2_MERGED_DIR, "{sample}.vcf.gz"),
        reference=REF,
        stats=os.path.join(MUTECT2_MERGED_DIR, "{sample}.vcf.gz.stats"),
        contamination=os.path.join(MUTECT2_CONTAM_DIR, "{sample}.contamination.table"),
        segmentation=os.path.join(MUTECT2_CONTAM_DIR, "{sample}.segments.table"),
        ob_priors=os.path.join(MUTECT2_MERGED_DIR, "{sample}_read-orientation-model.tar.gz"),
    output:
        vcf=ensure(
            protected(os.path.join(MUTECT2_FILTERED_DIR, "{sample}.filtered.vcf.gz")),
            non_empty=True,
        ),
    params:
        java_opts=get_java_opts,
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 17600 * attempt,
        runtime=lambda wildcards, attempt: 1440 * attempt,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "filter_mutect/{sample}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/filter_mutect/{sample}.tsv")
    shell:
        r"""
        gatk --java-options '{params.java_opts}' FilterMutectCalls \
            -R {input.reference} \
            -V {input.vcf} \
            --stats {input.stats} \
            --contamination-table {input.contamination} \
            --tumor-segmentation {input.segmentation} \
            --ob-priors {input.ob_priors} \
            -O {output.vcf} \
            2> {log}
        """
```

3.4. **Add wildcard constraints** at the top of `mutect2.smk`:
```python
wildcard_constraints:
    sample="[A-Za-z0-9_-]+",
    scatter_unit="(chr[0-9XYMI]+|[0-9]{4}-scattered|all)",
```

3.5. **Remove manual `os.makedirs()`** calls -- Snakemake handles this.

### Acceptance Criteria
- `snakemake -n --configfile config/config.yaml` produces correct DAG
- Scattered VCFs marked `temp()` + `ensure(non_empty=True)` (Snakemake 8+)
- Final outputs marked `protected()` + `ensure(non_empty=True)`
- All rules have `log:` and `benchmark:` directives
- All rules have `retries:` with attempt-based resource scaling
- Named input functions replace all inline lambdas (keep lambdas only for
  complex conditional logic like optional normal BAM)
- Conda env YAMLs used instead of named environments
- Rule docstrings describe purpose (shown in `--list-rules`)

---

## Phase 4: FreeBayes Pipeline Rules

**Goal:** Port `freebayes_germline_pipeline.smk` to `workflow/rules/freebayes.smk`.

### Tasks

4.1. **Write `workflow/rules/freebayes.smk`** with these rules:

| Rule | Key Changes |
|------|-------------|
| `freebayes_call` | Named input fns, `temp()` output, attempt-based resources (preserved from original), `benchmark:` |
| `merge_freebayes_vcfs` | `protected()` output, `benchmark:` |

4.2. **Shared scatter rules** from `scatter.smk` are reused (no duplication).

4.3. **Port the exponential resource scaling** (this was a good pattern):
```python
rule freebayes_call:
    resources:
        mem_mb=lambda wildcards, attempt: 16384 * (2 ** (attempt - 1)),
    retries: 3
```

4.4. **BAM list handling:** Derive from `samples_df` where
`analysis_type == "germline"` instead of separate `bam_list_file` config.

### Acceptance Criteria
- FreeBayes pipeline produces merged, normalized VCF
- Exponential retry resources preserved
- Scatter rules shared with Mutect2 (not duplicated)
- All rules have `log:`, `benchmark:`, `conda:`

---

## Phase 5: Profiles & Resource Management (Snakemake 8 Executor Plugin)

**Goal:** Externalize all resource allocation to profiles using the Snakemake 8
SLURM executor plugin. Replace the legacy `cluster:` string with
`executor: slurm` and SLURM-specific resource keys.

**Prerequisite:** `pip install snakemake-executor-plugin-slurm`

### Tasks

5.1. **Write `profiles/default/config.v8+.yaml`** (workflow-embedded profile):

This profile ships with the workflow and provides portable defaults.
Note: the `.v8+` suffix tells Snakemake to prefer this over a plain
`config.yaml` when running under Snakemake 8+.

```yaml
# Snakemake 8+ workflow profile for sm-calling
# This is the --workflow-profile (ships with the repo).
# Cluster-specific settings (executor, account, partition) go in
# a separate --profile (e.g., profiles/bih/ or profiles/charite/).

latency-wait: 60
jobs: 20
rerun-incomplete: true
printshellcmds: true
software-deployment-method:
  - conda

# Shared filesystem assumptions (typical HPC with Lustre/GPFS/NFS)
shared-fs-usage:
  - persistence
  - software-deployment
  - sources
  - source-cache
  - input-output

default-resources:
  - mem_mb: 4000
  - runtime: 120

set-threads:
  mutect2_call: 4
  gather_mutect2_vcfs: 4
  merge_mutect2_stats: 4
  learn_read_orientation: 4
  get_pileup_summaries: 4
  calculate_contamination: 4
  filter_mutect_calls: 4
  freebayes_call: 2
  merge_freebayes_vcfs: 12
  scatter_intervals_by_ns: 1
  split_intervals: 4

set-resources:
  mutect2_call:
    - mem_mb: 17600
    - runtime: 4320
  gather_mutect2_vcfs:
    - mem_mb: 17600
    - runtime: 1440
  merge_mutect2_stats:
    - mem_mb: 8000
    - runtime: 720
  learn_read_orientation:
    - mem_mb: 17600
    - runtime: 1440
  get_pileup_summaries:
    - mem_mb: 17600
    - runtime: 4320
  calculate_contamination:
    - mem_mb: 17600
    - runtime: 1440
  filter_mutect_calls:
    - mem_mb: 17600
    - runtime: 1440
  freebayes_call:
    - mem_mb: 16384
    - runtime: 4320
  merge_freebayes_vcfs:
    - mem_mb: 16384
    - runtime: 1440
```

5.2. **Write `profiles/bih/config.v8+.yaml`** (BIH cluster profile):
```yaml
# BIH HPC (cubi) -- uses the SLURM executor plugin
executor: slurm
jobs: 100
latency-wait: 60

default-resources:
  - slurm_partition: "medium"
  - slurm_account: "scholl-lab"
  - mem_mb_per_cpu: 4400
  - runtime: 120
```

5.3. **Write `profiles/charite/config.v8+.yaml`** (Charite cluster profile):
```yaml
# Charité HPC -- uses the SLURM executor plugin
executor: slurm
jobs: 50
latency-wait: 60

default-resources:
  - slurm_partition: "compute"
  - mem_mb_per_cpu: 4400
  - runtime: 120

# Charité-specific: partition overrides for long-running jobs
set-resources:
  mutect2_call:
    - slurm_partition: "long"
  get_pileup_summaries:
    - slurm_partition: "long"
  freebayes_call:
    - slurm_partition: "long"
```

5.4. **Write `profiles/local/config.v8+.yaml`** (local execution fallback):
```yaml
# Local execution (no cluster)
jobs: 4
software-deployment-method:
  - conda
```

5.5. **Remove hardcoded thread/memory values from rules** -- rules declare
attempt-based resource callables (for retry scaling) and minimum thread
counts. Profile `set-resources` and `set-threads` override at invocation.

5.6. **Invocation patterns:**
```bash
# BIH cluster (executor + workflow profile layered):
snakemake --profile profiles/bih --workflow-profile profiles/default

# Charite cluster:
snakemake --profile profiles/charite --workflow-profile profiles/default

# Local development:
snakemake --profile profiles/local --workflow-profile profiles/default -n

# The launcher script (Phase 6) handles this automatically.
```

### Acceptance Criteria
- `snakemake --workflow-profile profiles/default --profile profiles/bih` submits
  SLURM jobs via the executor plugin (no `sbatch` string needed)
- `snakemake --workflow-profile profiles/default --profile profiles/local -n`
  runs locally without SLURM
- Profile values match or exceed current hardcoded resource values
- `config.v8+.yaml` naming ensures Snakemake 8+ compatibility
- No `cluster:` string anywhere (fully migrated to executor plugin)

---

## Phase 6: Launcher & HPC Integration

**Goal:** Single universal launcher script using Snakemake 8 executor plugins
instead of manual SLURM wrapper logic. The launcher auto-detects the cluster
and selects the correct `--profile` (executor config) while always applying
the `--workflow-profile` (resource config).

### Tasks

6.1. **Port `scripts/run_snakemake.sh`** from sm-alignment with Snakemake 8 updates:
- Default SNAKEFILE: `workflow/Snakefile`
- Default CONFIG: `config/config.yaml`
- Cluster auto-detection sets `CLUSTER_PROFILE`:
  - `bih` → `--profile profiles/bih`
  - `charite` → `--profile profiles/charite`
  - `local` → `--profile profiles/local`
- Always applies: `--workflow-profile profiles/default`
- **No more `--use-conda`** -- handled by `software-deployment-method` in profile
- **No more `--cluster "sbatch ..."`** -- handled by `executor: slurm` in profile
- Keep TMPDIR management with trap cleanup
- Keep conda activation logic for Charite

```bash
#!/bin/bash
#SBATCH --job-name=sm_calling
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm_logs/%x-%j.log

set -euo pipefail

detect_cluster() {
    local fqdn
    fqdn=$(hostname -f 2>/dev/null || hostname)
    if [[ -d "/etc/xdg/snakemake/cubi-v1" ]] || [[ "$fqdn" =~ cubi|bihealth ]]; then
        echo "bih"
    elif [[ -f "/etc/profile.d/conda.sh" ]] || [[ "$fqdn" =~ charite|\.sc- ]]; then
        echo "charite"
    else
        echo "local"
    fi
}

CLUSTER=$(detect_cluster)
CONFIG_FILE="${1:-config/config.yaml}"
shift 1 2>/dev/null || true

# Conda activation (Charite requires explicit sourcing)
if [[ "$CLUSTER" == "charite" ]] && [[ -f /etc/profile.d/conda.sh ]]; then
    source /etc/profile.d/conda.sh
fi
conda activate snakemake 2>/dev/null || true

# TMPDIR setup
if [[ "$CLUSTER" == "bih" ]]; then
    BASE_TMPDIR="${HOME}/scratch/tmp"
else
    BASE_TMPDIR="${TMPDIR:-/tmp}/snakemake"
fi
mkdir -p "${BASE_TMPDIR}"
TMPDIR=$(mktemp -d "${BASE_TMPDIR}/sm.XXXXXX")
export TMPDIR
trap 'rm -rf "${TMPDIR}"' EXIT

mkdir -p slurm_logs

echo "=== sm-calling Launch ==="
echo "  Cluster:    ${CLUSTER}"
echo "  Config:     ${CONFIG_FILE}"
echo "  Profile:    profiles/${CLUSTER}"
echo "  TMPDIR:     ${TMPDIR}"
echo "  Extra args: $*"
echo "  Start:      $(date)"
echo "========================="

snakemake \
    -s workflow/Snakefile \
    --configfile "${CONFIG_FILE}" \
    --workflow-profile profiles/default \
    --profile "profiles/${CLUSTER}" \
    "$@"

echo "=== Finished: $(date) ==="
```

6.2. **Usage:**
```bash
# Default (auto-detects cluster, uses config/config.yaml):
sbatch scripts/run_snakemake.sh

# Custom config:
sbatch scripts/run_snakemake.sh config/config_exomes.yaml

# Dry-run (local):
bash scripts/run_snakemake.sh config/config.yaml -n

# Extra Snakemake flags:
sbatch scripts/run_snakemake.sh config/config.yaml --forcerun mutect2_call
```

6.3. Optionally port `scripts/generate_config.py` from sm-alignment, adapted
for calling-specific metadata (tumor/normal BAM columns, analysis type).

### Acceptance Criteria
- Single launcher works on BIH, Charite, and local
- Auto-detects cluster and selects correct `--profile` (executor plugin)
- Uses `--workflow-profile` for resource allocation (portable)
- No `--use-conda`, no `--cluster` flags (pure Snakemake 8)
- Replaces all 7 existing launcher scripts
- TMPDIR setup and cleanup works correctly

---

## Phase 7: Testing & CI/CD

**Goal:** Add tests and automated quality checks.

### Tasks

7.1. **Write `pyproject.toml`** (tool config):
```toml
[tool.ruff]
line-length = 99
select = ["E", "W", "F", "I", "UP", "B", "SIM", "RUF"]

[tool.mypy]
python_version = "3.10"
warn_return_any = true
warn_unused_configs = true

[tool.snakefmt]
line_length = 99

[tool.pytest.ini_options]
testpaths = ["tests"]
markers = ["dryrun: marks tests as dry-run integration tests"]
```

7.2. **Write `Makefile`**:
```makefile
.PHONY: lint format test test-unit test-dryrun

lint:
    ruff check workflow/rules/helpers.py scripts/
    snakefmt --check workflow/

format:
    ruff format workflow/rules/helpers.py scripts/
    snakefmt workflow/

test:
    pytest tests/ -v

test-unit:
    pytest tests/ -v -m "not dryrun"

test-dryrun:
    pytest tests/ -v -m dryrun
```

7.3. **Write `tests/conftest.py`** with fixtures for sample DataFrames and
config dicts.

7.4. **Write `tests/test_helpers.py`**:
- `test_get_java_opts` -- verifies 80/20 split
- `test_get_mutect2_samples` -- filters correctly
- `test_get_germline_samples` -- filters correctly
- `test_get_unique_bam_basenames` -- deduplicates, ignores "."
- `test_has_matched_normal` -- True/False cases
- `test_get_scatter_units_chromosome` -- returns chromosome list
- `test_get_scatter_units_interval` -- returns 0000-0399 formatted
- `test_get_scatter_units_none` -- returns ["all"]
- `test_build_freebayes_params` -- flag vs key-value

7.5. **Write `tests/test_schema_validation.py`**:
- Valid config passes
- Missing `caller` fails
- Invalid `scatter.mode` fails
- Valid samples.tsv passes
- Missing required columns fails

7.6. **Write `tests/test_dryrun.py`** (marked `@pytest.mark.dryrun`):
- Dry-run with Mutect2 caller produces expected targets
- Dry-run with FreeBayes caller produces expected targets
- Dry-run with "all" caller produces both

7.7. **Write `.github/workflows/main.yaml`** for CI:
- Job: lint (ruff + snakefmt)
- Job: test-unit (pytest, no dryrun)
- Job: test-dryrun (requires conda, slower)

### Acceptance Criteria
- `make lint` passes
- `make test-unit` passes (all helper tests green)
- `make test-dryrun` passes with stub config + samples
- CI runs on push to main and PRs

---

## Phase 8: Documentation & Cleanup

**Goal:** Update all documentation, add CLAUDE.md, finalize deprecation.

### Tasks

8.1. **Rewrite `README.md`** covering:
- Quick start (3 commands to run)
- Configuration reference
- Sample sheet format
- Supported callers (Mutect2, FreeBayes)
- HPC cluster setup (BIH, Charite)
- Development (lint, test, format)

8.2. **Write `CLAUDE.md`** (AI assistant guidance):
- Architecture overview
- Key design decisions
- File structure
- How to add a new caller
- Testing instructions

8.3. **Write `config/README.md`** (configuration instructions, required by catalog).

8.4. **Write `.snakemake-workflow-catalog.yml`**:
```yaml
usage:
  software-stack-deployment:
    conda: true
    apptainer: false
  mandatory-flags:
    desc: ""
    flags: ""
```

8.5. **Add deprecation notice** in `deprecated/README.md`:
```markdown
# Deprecated Files

These files are from the pre-refactoring layout. They are preserved for
reference but are no longer maintained.

**Removal target:** 3 months after refactoring merge.

## Migration

See the main README.md for the new workflow structure.
```

8.6. **Clean up `.gitignore`** to cover all generated artifacts.

### Acceptance Criteria
- README has complete setup and usage instructions
- CLAUDE.md accurately describes the architecture
- Workflow catalog YAML present and valid
- Deprecated files have clear removal timeline

---

## Migration Guide

### For Existing Users

1. **Snakemake version:** Upgrade to Snakemake 8+ and install the SLURM plugin:
   ```bash
   conda create -n snakemake -c conda-forge -c bioconda snakemake>=8.0
   pip install snakemake-executor-plugin-slurm
   ```

2. **Config migration:** Map your existing settings:
   ```
   Old: config.yaml (root)                    → New: config/config.yaml
   Old: config/config_mutect2.yaml            → Merged into config/config.yaml
   Old: config/config_freebayes_germline.yaml → Merged into config/config.yaml
   ```

3. **Metadata migration:**
   ```
   Old: calling_metadata.tsv columns:
     individual1, individual2, analysis, sample1, sample2,
     bam1_file_basename, bam2_file_basename

   New: config/samples.tsv columns:
     sample, tumor_bam, normal_bam, analysis_type
   ```

4. **Invocation change (Snakemake 8 style):**
   ```bash
   # Old (4 separate submissions + manual orchestration):
   sbatch scripts/run_mutect2_calling.sh
   sbatch scripts/run_merge_mutect2_calls.sh
   sbatch scripts/run_gatk_calculate_contamination.sh
   sbatch scripts/run_gatk_filter_mutect_calls.sh

   # New (single submission, auto-detects cluster):
   sbatch scripts/run_snakemake.sh

   # Or explicitly with config:
   sbatch scripts/run_snakemake.sh config/config.yaml

   # Or directly (without launcher wrapper):
   snakemake --profile profiles/bih --workflow-profile profiles/default \
       --configfile config/config.yaml
   ```

5. **CLI flag changes (Snakemake 7 → 8):**
   ```bash
   # Old:
   snakemake --use-conda --cluster "sbatch ..." --profile cubi-v1

   # New:
   snakemake --sdm conda --profile profiles/bih --workflow-profile profiles/default
   # (executor: slurm is in the profile, no --cluster needed)
   ```

6. **Output paths:**
   ```
   Old: results/variant_calls/, results/variant_merge/, etc.
   New: results/exomes/variant_calls/mutect2/scattered/, .../merged/, etc.
   ```

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Breaking existing runs | Medium | High | Keep deprecated/ with old files; parallel operation period |
| Snakemake 8 not installed on cluster | Medium | High | Verify Snakemake version on both BIH and Charite before starting; provide install instructions |
| SLURM executor plugin compatibility | Medium | High | Test `snakemake-executor-plugin-slurm` on both clusters early (Phase 5); fall back to `cluster:` string in Charite profile if needed |
| `lookup()`/`collect()` API changes | Low | Medium | Pin minimum Snakemake version; these APIs are stable since 8.0 |
| Schema too strict | Low | Medium | Start with minimal required fields; iterate |
| Conda env version conflicts | Low | Medium | Test pinned versions on both clusters before merge |
| Resource allocation changes | Medium | Low | Profile values match current hardcoded values exactly |
| Metadata format change | High | Medium | Provide conversion script or document mapping |
| `ensure()` false positives | Low | Low | Only use `non_empty=True`; skip checksum validation |

### Recommended Execution Order

Phases 1-2 should be done first (foundation). Phases 3-4 can be done in
parallel. Phases 5-6 depend on 3-4. Phase 7 can start after Phase 2
(test helpers early). Phase 8 is last.

```
Phase 1 ─→ Phase 2 ─┬→ Phase 3 ─┬→ Phase 5 ─→ Phase 6 ─→ Phase 8
                     │           │
                     ├→ Phase 4 ─┘
                     │
                     └→ Phase 7 (can start early, expand later)
```

### Estimated Effort

| Phase | Effort | Complexity |
|-------|--------|------------|
| Phase 1: Directory restructure | Small | Low |
| Phase 2: Common infrastructure | Medium | Medium |
| Phase 3: Mutect2 rules | Large | High |
| Phase 4: FreeBayes rules | Medium | Medium |
| Phase 5: Profiles | Small | Low |
| Phase 6: Launcher | Small | Low |
| Phase 7: Testing | Medium | Medium |
| Phase 8: Documentation | Medium | Low |

---

## Appendix A: Patterns Adopted from sm-alignment

| Pattern | sm-alignment Implementation | sm-calling Adaptation |
|---------|----------------------------|----------------------|
| Main Snakefile | 38 lines: validate + include + rule all | Same pattern, adds caller dispatch |
| common.smk | Config shortcuts + samples_df + helper wrappers | Same, adds caller-specific helpers |
| helpers.py | Pure Python, no Snakemake imports | Same, adds BAM resolution + scatter logic |
| Schema validation | `validate(config, ...)` at startup | Same |
| Conda env YAMLs | `workflow/envs/` with pinned versions | Same |
| Profile-based resources | `profiles/default/` + `profiles/charite/` | Same structure, calling-specific values |
| Universal launcher | `scripts/run_snakemake.sh` with cluster detection | Ported directly |
| temp() intermediates | All BAMs except final marked temp() | Scattered VCFs + stats marked temp() |
| Naming conventions | `{sample}` wildcard from `project_sample` | `{sample}` from `sample` column |
| Log structure | `logs/{rule_name}.{tool}.{sample}.log` | `logs/{rule_name}/{sample}_{scatter_unit}.log` |
| Java opts | 80/20 heap split from resources.mem_mb | Same function, imported from helpers |
| Tests | pytest with conftest fixtures | Same structure |
| Code quality | ruff + mypy + snakefmt via Makefile | Same tooling |

## Appendix B: Snakemake 8+ Features Used (Beyond sm-alignment)

These features are adopted in sm-calling but were NOT used in the sm-alignment
refactoring. They represent the next evolution of the lab's workflow standards.

| Feature | Purpose | Where Used |
|---------|---------|------------|
| `min_version("8.0")` | Enforce Snakemake 8+ at startup | `workflow/Snakefile` |
| `executor: slurm` | Native SLURM plugin (no `cluster:` string) | `profiles/bih/`, `profiles/charite/` |
| `config.v8+.yaml` | Version-specific profile naming | All profile directories |
| `localrule: True` | Per-rule local execution (replaces global `localrules:`) | `rule all` |
| `lookup()` | Declarative metadata lookups | `common.smk` target generation |
| `collect()` | Declarative target aggregation | `get_final_outputs()` |
| `ensure(non_empty=True)` | Output validation (catches truncated files) | `mutect2_call`, `filter_mutect_calls` |
| `shared-fs-usage` | Granular shared filesystem assumptions | `profiles/default/` |
| Attempt-based `resources:` | Auto-scaling memory/runtime on retry | All GATK rules |
| `retries:` | Per-rule automatic retry | All variant calling rules |
| `slurm_partition` resource | SLURM-specific partition targeting | `profiles/charite/` |
| `slurm_account` resource | SLURM-specific account targeting | `profiles/bih/` |
| `--rerun-triggers` | Control re-execution triggers | Documented in profiles |

### Snakemake 8 Dependency

```
snakemake >= 8.0
snakemake-executor-plugin-slurm (for cluster execution)
```

Install via:
```bash
pip install snakemake snakemake-executor-plugin-slurm
# or
conda install -c conda-forge -c bioconda snakemake>=8.0
pip install snakemake-executor-plugin-slurm
```
