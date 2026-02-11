# Deprecated Files

These files are from the pre-refactoring layout. They are preserved for
reference but are no longer maintained.

**Removal target:** 3 months after refactoring merge.

## Structure

- `scripts/snakemake/` -- Old standalone `.smk` files (7 workflows)
- `scripts/launchers/` -- Old SLURM launcher scripts (7 scripts)
- `configs/` -- Old configuration files (3 configs)

## Migration

See the main [README.md](../README.md) for the new workflow structure.

The new layout uses:
- `workflow/Snakefile` as the single entry point
- `config/config.yaml` as the unified configuration
- `config/samples.tsv` as the unified sample sheet
- `scripts/run_snakemake.sh` as the universal launcher
