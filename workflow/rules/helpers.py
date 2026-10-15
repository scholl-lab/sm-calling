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


def get_tumor_bam_basenames(samples_df: pd.DataFrame, sample: str) -> list[str]:
    """Return tumor BAM basenames for a sample (splitting ';'-delimited multi-region entries)."""
    row = samples_df.loc[samples_df["sample"] == sample].iloc[0]
    raw = str(row.get("tumor_bam", ""))
    return [b.strip() for b in raw.split(";") if b.strip() and b.strip() != "."]


def get_unique_bam_basenames(samples_df: pd.DataFrame) -> list[str]:
    """Return all unique BAM basenames across tumor and normal columns."""
    bams = set()
    for _, row in samples_df.iterrows():
        if row["tumor_bam"] and row["tumor_bam"] != ".":
            for b in str(row["tumor_bam"]).split(";"):
                b = b.strip()
                if b and b != ".":
                    bams.add(b)
        if (
            row.get("analysis_type") == "tumor_normal"
            and pd.notna(row.get("normal_bam"))
            and row["normal_bam"] not in ("", ".")
        ):
            bams.add(row["normal_bam"])
    return sorted(bams)


def has_matched_normal(samples_df: pd.DataFrame, sample: str) -> bool:
    """Check whether a sample has a matched normal BAM."""
    row = samples_df.loc[samples_df["sample"] == sample].iloc[0]
    if row.get("analysis_type") != "tumor_normal":
        return False
    normal = row.get("normal_bam", "")
    return bool(normal) and normal != "." and pd.notna(normal)


def validate_samples_df(samples_df: pd.DataFrame) -> None:
    """Validate sample sheet invariants beyond per-row JSON schema checks."""
    if samples_df["sample"].duplicated().any() or samples_df.index.duplicated().any():
        dups = sorted(
            samples_df.loc[samples_df["sample"].duplicated(), "sample"].astype(str).unique()
        )
        raise ValueError(f"Duplicate sample IDs in samples.tsv: {dups}")

    for _, row in samples_df.iterrows():
        sample = row.get("sample", "")
        tumor = row.get("tumor_bam", "")
        normal = row.get("normal_bam", "")
        analysis_type = row.get("analysis_type", "")

        if not isinstance(sample, str) or not isinstance(tumor, str):
            raise ValueError(
                "Sample and BAM identifiers must be strings (read samples.tsv with dtype=str)"
            )
        tumor_parts = [b.strip() for b in tumor.split(";")]
        if not sample or sample == "." or not tumor or any(not b or b == "." for b in tumor_parts):
            raise ValueError(f"Sample '{sample}' has empty or '.' sample/tumor_bam")

        has_normal_bam = bool(normal) and normal != "." and pd.notna(normal)
        if has_normal_bam and not isinstance(normal, str):
            raise ValueError("normal_bam must be a string (read samples.tsv with dtype=str)")
        if analysis_type == "tumor_normal" and not has_normal_bam:
            raise ValueError(
                f"Sample '{sample}' has analysis_type='tumor_normal' but missing normal_bam"
            )
        if analysis_type != "tumor_normal" and has_normal_bam:
            raise ValueError(
                f"Sample '{sample}' has analysis_type='{analysis_type}' "
                f"but normal_bam='{normal}' is provided"
            )


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


def build_interval_arg(
    mode: str,
    unit: str,
    intervals_dir: str = "",
    target_bed: str = "",
) -> str:
    """Build GATK -L interval arguments with optional target_bed intersection."""
    if mode == "chromosome":
        if target_bed:
            return f"-L {unit} -L {target_bed} --interval-set-rule INTERSECTION"
        return f"-L {unit}"
    elif mode == "interval":
        if target_bed:
            return (
                f"-L {intervals_dir}/{unit}.interval_list "
                f"-L {target_bed} --interval-set-rule INTERSECTION"
            )
        return f"-L {intervals_dir}/{unit}.interval_list"
    elif target_bed:
        return f"-L {target_bed}"
    return ""


def build_mutect2_extra_args(
    genotype_germline_sites: bool = True,
    genotype_pon_sites: bool = True,
    disable_adaptive_pruning: bool = True,
    annotations: list[str] | None = None,
    annotation_groups: list[str] | None = None,
    extra: str = "",
) -> str:
    """Build the complete Mutect2 extra arguments string.

    Combines dedicated config fields with the passthrough ``extra`` string.
    """
    args: list[str] = []
    if genotype_germline_sites:
        args.append("--genotype-germline-sites true")
    if genotype_pon_sites:
        args.append("--genotype-pon-sites true")
    if disable_adaptive_pruning:
        args.append("--disable-adaptive-pruning true")
    for ann in annotations or []:
        args.append(f"--annotation {ann}")
    for grp in annotation_groups or []:
        args.append(f"--annotation-group {grp}")
    if extra:
        args.append(extra)
    return " ".join(args)


def get_normal_bam_basenames(samples_df: pd.DataFrame) -> list[str]:
    """Return sorted basenames of normal BAMs (for PoN and PureCN NormalDB)."""
    has_in_pon = "in_pon" in samples_df.columns
    normals: set[str] = set()
    for _, row in samples_df.iterrows():
        in_pon = str(row.get("in_pon", "")).strip().lower() in (
            "yes",
            "true",
            "1",
        )
        if has_in_pon and not in_pon:
            continue
        analysis_type = row.get("analysis_type")
        normal = row.get("normal_bam", "")
        if analysis_type == "tumor_normal" and normal and normal != "." and pd.notna(normal):
            normals.add(normal)
        elif has_in_pon and in_pon and analysis_type == "germline":
            for b in str(row.get("tumor_bam", "")).split(";"):
                b = b.strip()
                if b and b != ".":
                    normals.add(b)
    return sorted(normals)


def build_freebayes_params(extra: str) -> str:
    """Return the FreeBayes extra params string (passthrough from config)."""
    return extra
