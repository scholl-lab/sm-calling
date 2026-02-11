"""Pure helper functions for the sm-calling workflow.

All functions are free of Snakemake imports so they can be unit-tested.
"""

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
    return bool(normal) and normal != "." and pd.notna(normal)


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


def build_freebayes_params(extra: str) -> str:
    """Return the FreeBayes extra params string (passthrough from config)."""
    return extra
