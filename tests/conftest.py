"""Shared pytest fixtures for sm-calling tests."""

import pandas as pd
import pytest


@pytest.fixture
def samples_df():
    """Sample DataFrame matching config/samples.tsv schema."""
    return pd.DataFrame(
        {
            "sample": ["IND001_To", "IND002_TN", "IND003_G"],
            "tumor_bam": ["IND001.tumor", "IND002.tumor", "IND003"],
            "normal_bam": [".", "IND002.normal", "."],
            "analysis_type": ["tumor_only", "tumor_normal", "germline"],
        }
    ).set_index("sample", drop=False)


@pytest.fixture
def config_dict():
    """Minimal valid configuration dictionary."""
    return {
        "caller": "mutect2",
        "ref": {
            "genome": "/ref/GRCh38.fna",
            "build": "GRCh38",
        },
        "gatk_resources": {
            "panel_of_normals": "/res/pon.vcf.gz",
            "af_only_gnomad": "/res/gnomad.vcf.gz",
            "common_biallelic_gnomad": "/res/gnomad_common.vcf.gz",
        },
        "paths": {
            "samples": "config/samples.tsv",
            "bam_folder": "/bams",
            "output_folder": "/output",
            "log_subdir": "logs",
        },
        "bam": {
            "file_extension": ".merged.dedup.bqsr.bam",
        },
        "scatter": {
            "mode": "chromosome",
            "count": 400,
            "chromosomes": [f"chr{i}" for i in range(1, 23)]
            + ["chrX", "chrY", "chrM"],
        },
        "params": {
            "mutect2": {"extra": "--genotype-germline-sites true"},
            "freebayes": {"extra": "--min-coverage 20"},
            "bcftools_norm": {"extra": "-m-any --force"},
        },
    }
