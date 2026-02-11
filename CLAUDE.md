# CLAUDE.md

## Project Overview

sm-calling is a Snakemake 8+ variant calling pipeline supporting Mutect2
(somatic) and FreeBayes (germline). It runs on BIH and Charite HPC clusters
via the SLURM executor plugin.

## Architecture

- `workflow/Snakefile` -- entry point: validates config + samples, includes rules
- `workflow/rules/common.smk` -- config shortcuts, samples_df, input helper functions
- `workflow/rules/helpers.py` -- pure Python (no Snakemake imports), unit-testable
- `workflow/rules/scatter.smk` -- interval scattering (shared by both callers)
- `workflow/rules/mutect2.smk` -- full Mutect2 best-practice pipeline
- `workflow/rules/freebayes.smk` -- FreeBayes germline pipeline

## Key Design Decisions

- **Snakemake 8+**: `min_version("8.0")`, `localrule: True`, `ensure()`,
  `retries:` with attempt-based resource scaling
- **Profile layering**: `--workflow-profile profiles/default` (resources) +
  `--profile profiles/{cluster}` (executor)
- **SLURM executor plugin**: `executor: slurm` in profiles, no `--cluster` string
- **Schema validation**: `validate()` for both config and samples at startup
- **Conda env YAMLs**: Pinned versions in `workflow/envs/`, not named envs
- **Pure helpers**: `helpers.py` has no Snakemake imports for testability

## Adding a New Caller

1. Create `workflow/rules/newcaller.smk`
2. Add caller-specific input functions to `common.smk`
3. Update `get_final_outputs()` in `common.smk` to handle the new caller
4. Add `include: "rules/newcaller.smk"` to `workflow/Snakefile`
5. Add conda env YAML to `workflow/envs/`
6. Update `config.schema.yaml` caller enum
7. Add tests to `tests/`

## Testing

```bash
pytest tests/ -v             # all tests
pytest tests/ -v -m "not dryrun"  # unit tests only
```

## File Naming Conventions

- Rules: `{sample}_{scatter_unit}` wildcards
- Outputs: `{output_dir}/{caller}/{stage}/{sample}.{ext}`
- Logs: `{output_dir}/logs/{rule_name}/{sample}_{scatter_unit}.log`
- Benchmarks: `{output_dir}/logs/benchmarks/{rule_name}/{sample}_{scatter_unit}.tsv`
