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

## `samples.tsv`

Tab-separated sample metadata. Required columns:

| Column | Description |
|--------|-------------|
| `sample` | Unique sample/analysis identifier |
| `tumor_bam` | BAM file basename (without extension) |
| `normal_bam` | Matched normal BAM basename, or `.` if none |
| `analysis_type` | `tumor_only`, `tumor_normal`, or `germline` |
