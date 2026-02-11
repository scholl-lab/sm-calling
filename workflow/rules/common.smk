"""Shared configuration shortcuts, metadata accessors, and helper wrappers."""

import os

import pandas as pd

from rules.helpers import (
    get_java_opts as _get_java_opts_impl,
    get_mutect2_samples as _get_mutect2_samples_impl,
    get_germline_samples as _get_germline_samples_impl,
    get_unique_bam_basenames as _get_unique_bam_basenames_impl,
    has_matched_normal as _has_matched_normal_impl,
    get_scatter_units as _get_scatter_units_impl,
)


# -- Config shortcuts --------------------------------------------------------
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

INTERVALS_DIR = config["paths"].get("intervals_dir", os.path.join(OUTPUT_DIR, "intervals"))

# GATK resources (optional for FreeBayes-only runs)
GATK_RES = config.get("gatk_resources", {})
PON = GATK_RES.get("panel_of_normals", "")
AF_GNOMAD = GATK_RES.get("af_only_gnomad", "")
COMMON_GNOMAD = GATK_RES.get("common_biallelic_gnomad", "")

# Tool params
MUTECT2_EXTRA = config.get("params", {}).get("mutect2", {}).get("extra", "")
FREEBAYES_EXTRA = config.get("params", {}).get("freebayes", {}).get("extra", "")
BCFTOOLS_NORM_EXTRA = config.get("params", {}).get("bcftools_norm", {}).get("extra", "")


# -- Output subdirectories ---------------------------------------------------
MUTECT2_DIR = os.path.join(OUTPUT_DIR, "mutect2")
MUTECT2_SCATTER_DIR = os.path.join(MUTECT2_DIR, "scattered")
MUTECT2_MERGED_DIR = os.path.join(MUTECT2_DIR, "merged")
MUTECT2_CONTAM_DIR = os.path.join(MUTECT2_DIR, "contamination")
MUTECT2_FILTERED_DIR = os.path.join(MUTECT2_DIR, "filtered")

FREEBAYES_DIR = os.path.join(OUTPUT_DIR, "freebayes")
FREEBAYES_SCATTER_DIR = os.path.join(FREEBAYES_DIR, "scattered")


# -- Scatter units -----------------------------------------------------------
SCATTER_UNITS = _get_scatter_units_impl(SCATTER_MODE, CHROMOSOMES, SCATTER_COUNT)


# -- Sample accessors (wrappers around helpers.py) ---------------------------
def get_mutect2_sample_list():
    return _get_mutect2_samples_impl(samples_df)


def get_germline_sample_list():
    return _get_germline_samples_impl(samples_df)


def get_unique_bams():
    return _get_unique_bam_basenames_impl(samples_df)


def has_normal(sample):
    return _has_matched_normal_impl(samples_df, sample)


# -- Input helper functions --------------------------------------------------
# Use lookup() for simple column retrieval (Snakemake 8+ idiom).
# Keep named functions for complex conditional logic (e.g., optional normals).


def get_tumor_bam(wildcards):
    """Return full path to tumor BAM for a sample."""
    row = samples_df.loc[wildcards.sample]
    return os.path.join(BAM_FOLDER, row["tumor_bam"] + BAM_EXT)


def get_normal_bam(wildcards):
    """Return list with normal BAM path, or empty list if none.

    Cannot use lookup() because of conditional empty-list logic.
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
    """Return -L argument based on scatter mode and unit."""
    unit = wildcards.scatter_unit
    if SCATTER_MODE == "chromosome":
        return f"-L {unit}"
    elif SCATTER_MODE == "interval":
        return f"-L {INTERVALS_DIR}/{unit}.interval_list"
    return ""


def get_freebayes_region_arg(wildcards):
    """Return FreeBayes region argument based on scatter mode."""
    unit = wildcards.scatter_unit
    if SCATTER_MODE == "chromosome":
        return f"-r {unit}"
    elif SCATTER_MODE == "interval":
        return f"--targets {INTERVALS_DIR}/{unit}.bed"
    return ""


def get_java_opts(wildcards, resources):
    """Derive GATK --java-options from allocated resources."""
    return _get_java_opts_impl(resources.mem_mb, resources.tmpdir)


def get_tumor_pileup(wildcards):
    """Return pileup table path for the tumor BAM of a sample."""
    row = samples_df.loc[wildcards.sample]
    return os.path.join(
        MUTECT2_CONTAM_DIR,
        row["tumor_bam"] + ".getpileupsummaries.table",
    )


def get_normal_pileup(wildcards):
    """Return list with normal pileup path, or empty list if none."""
    row = samples_df.loc[wildcards.sample]
    normal = row.get("normal_bam", "")
    if normal and normal != "." and pd.notna(normal):
        return [
            os.path.join(
                MUTECT2_CONTAM_DIR,
                normal + ".getpileupsummaries.table",
            )
        ]
    return []


# -- Final output dispatcher -------------------------------------------------
def get_final_outputs():
    """Return all expected final output files based on caller selection."""
    outputs = []
    if CALLER in ("mutect2", "all"):
        for sample in get_mutect2_sample_list():
            outputs.append(os.path.join(MUTECT2_FILTERED_DIR, f"{sample}.filtered.vcf.gz"))
    if CALLER in ("freebayes", "all"):
        outputs.append(os.path.join(FREEBAYES_DIR, "final_merged.vcf.gz"))
    return outputs
