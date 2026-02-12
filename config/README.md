# Configuration

## `config.yaml`

Unified pipeline configuration. All paths are relative to the repository root.

### Required Fields

| Field | Type | Description |
|-------|------|-------------|
| `caller` | string | `"mutect2"`, `"freebayes"`, or `"all"` |
| `ref.genome` | string | Path to reference genome FASTA |
| `ref.build` | string | `"GRCh37"` or `"GRCh38"` |
| `paths.samples` | string | Path to samples TSV |
| `paths.bam_folder` | string | Directory with input BAMs |
| `paths.output_folder` | string | Root output directory |
| `bam.file_extension` | string | BAM file suffix (e.g., `.merged.dedup.bqsr.bam`) |
| `scatter.mode` | string | `"chromosome"`, `"interval"`, or `"none"` |

### Required for Mutect2

| Field | Description |
|-------|-------------|
| `gatk_resources.panel_of_normals` | Path to Panel of Normals VCF |
| `gatk_resources.af_only_gnomad` | Path to gnomAD allele frequency VCF |
| `gatk_resources.common_biallelic_gnomad` | Path to gnomAD common biallelic VCF |

### Mutect2 Parameters (`params.mutect2`)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `genotype_germline_sites` | boolean | `true` | Emit germline sites (required for PureCN) |
| `genotype_pon_sites` | boolean | `true` | Emit PoN sites (required for PureCN) |
| `annotations` | array of strings | `[]` | Extra `--annotation` flags for Mutect2 |
| `annotation_groups` | array of strings | `[]` | Extra `--annotation-group` flags for Mutect2 |
| `extra` | string | `""` | Passthrough for any other Mutect2 flags |

### QC Parameters (`params.bcftools_stats`)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `extra` | string | `""` | Extra flags for `bcftools stats` |

### PureCN Settings (`purecn`)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | boolean | `false` | Enable PureCN copy number analysis |
| `genome` | string | `"hg38"` | PureCN genome identifier (`"hg19"` or `"hg38"`) |
| `intervals_bed` | string | `""` | BED file with capture bait coordinates |
| `normaldb` | string | `""` | Pre-built `normalDB.rds` (skip NormalDB.R if provided) |
| `mapping_bias` | string | `""` | Pre-built `mapping_bias.rds` |
| `snp_blacklist` | string | `""` | Optional simple repeats BED for SNP filtering |
| `extra` | string | `""` | Extra PureCN.R arguments |
| `seed` | integer | `123` | Random seed for PureCN |
| `postoptimize` | boolean | `true` | Run PureCN post-optimization |

## `samples.tsv`

Tab-separated sample metadata. Required columns:

| Column | Description |
|--------|-------------|
| `sample` | Unique sample/analysis identifier |
| `tumor_bam` | BAM file basename (without extension) |
| `normal_bam` | Matched normal BAM basename, or `.` if none |
| `analysis_type` | `tumor_only`, `tumor_normal`, or `germline` |
