"""Unit tests for workflow/rules/helpers.py."""

import sys
import os

# Add workflow/rules to the path so helpers.py can be imported directly
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "rules"))

from helpers import (
    build_freebayes_params,
    build_mutect2_extra_args,
    get_germline_samples,
    get_java_opts,
    get_mutect2_samples,
    get_normal_bam_basenames,
    get_scatter_units,
    get_unique_bam_basenames,
    has_matched_normal,
)


class TestGetJavaOpts:
    def test_standard_allocation(self):
        result = get_java_opts(17600, "/tmp")
        assert "-Xms3520m" in result
        assert "-Xmx14080m" in result
        assert "-Djava.io.tmpdir=/tmp" in result

    def test_small_allocation(self):
        result = get_java_opts(1000, "/scratch")
        assert "-Xms200m" in result
        assert "-Xmx800m" in result


class TestGetMutect2Samples:
    def test_filters_correctly(self, samples_df):
        result = get_mutect2_samples(samples_df)
        assert result == ["IND001_To", "IND002_TN"]

    def test_excludes_germline(self, samples_df):
        result = get_mutect2_samples(samples_df)
        assert "IND003_G" not in result


class TestGetGermlineSamples:
    def test_filters_correctly(self, samples_df):
        result = get_germline_samples(samples_df)
        assert result == ["IND003_G"]

    def test_excludes_somatic(self, samples_df):
        result = get_germline_samples(samples_df)
        assert "IND001_To" not in result
        assert "IND002_TN" not in result


class TestGetUniqueBamBasenames:
    def test_deduplicates(self, samples_df):
        result = get_unique_bam_basenames(samples_df)
        assert "IND001.tumor" in result
        assert "IND002.tumor" in result
        assert "IND002.normal" in result
        assert "IND003" in result
        assert len(result) == 4

    def test_ignores_dot(self, samples_df):
        result = get_unique_bam_basenames(samples_df)
        assert "." not in result


class TestHasMatchedNormal:
    def test_tumor_normal(self, samples_df):
        assert has_matched_normal(samples_df, "IND002_TN") is True

    def test_tumor_only(self, samples_df):
        assert has_matched_normal(samples_df, "IND001_To") is False

    def test_germline(self, samples_df):
        assert has_matched_normal(samples_df, "IND003_G") is False


class TestGetScatterUnits:
    def test_chromosome_mode(self):
        chroms = ["chr1", "chr2", "chrX"]
        result = get_scatter_units("chromosome", chroms)
        assert result == chroms

    def test_interval_mode(self):
        result = get_scatter_units("interval", scatter_count=5)
        assert result == [
            "0000-scattered",
            "0001-scattered",
            "0002-scattered",
            "0003-scattered",
            "0004-scattered",
        ]

    def test_none_mode(self):
        result = get_scatter_units("none")
        assert result == ["all"]

    def test_chromosome_mode_empty(self):
        result = get_scatter_units("chromosome", None)
        assert result == []


class TestBuildMutect2ExtraArgs:
    def test_defaults(self):
        result = build_mutect2_extra_args()
        assert "--genotype-germline-sites true" in result
        assert "--genotype-pon-sites true" in result

    def test_disabled(self):
        result = build_mutect2_extra_args(
            genotype_germline_sites=False,
            genotype_pon_sites=False,
        )
        assert result == ""

    def test_annotations(self):
        result = build_mutect2_extra_args(
            genotype_germline_sites=False,
            genotype_pon_sites=False,
            annotations=["OrientationBiasReadCounts", "StrandBiasBySample"],
            annotation_groups=["AS_StandardAnnotation"],
        )
        assert "--annotation OrientationBiasReadCounts" in result
        assert "--annotation StrandBiasBySample" in result
        assert "--annotation-group AS_StandardAnnotation" in result

    def test_with_extra(self):
        result = build_mutect2_extra_args(
            genotype_germline_sites=True,
            genotype_pon_sites=False,
            extra="--max-reads-per-alignment-start 50",
        )
        assert "--genotype-germline-sites true" in result
        assert "--genotype-pon-sites" not in result
        assert "--max-reads-per-alignment-start 50" in result

    def test_all_combined(self):
        result = build_mutect2_extra_args(
            genotype_germline_sites=True,
            genotype_pon_sites=True,
            annotations=["OrientationBiasReadCounts"],
            annotation_groups=["AS_StandardAnnotation"],
            extra="--custom-flag",
        )
        parts = result.split()
        # Verify ordering: germline, pon, annotations, annotation-groups, extra
        germline_idx = parts.index("--genotype-germline-sites")
        pon_idx = parts.index("--genotype-pon-sites")
        ann_idx = parts.index("--annotation")
        grp_idx = parts.index("--annotation-group")
        custom_idx = parts.index("--custom-flag")
        assert germline_idx < pon_idx < ann_idx < grp_idx < custom_idx


class TestGetNormalBamBasenames:
    def test_extracts_normals(self, samples_df):
        result = get_normal_bam_basenames(samples_df)
        assert result == ["IND002.normal"]

    def test_ignores_dots(self, samples_df):
        result = get_normal_bam_basenames(samples_df)
        assert "." not in result

    def test_empty_when_no_normals(self):
        import pandas as pd

        df = pd.DataFrame(
            {
                "sample": ["S1"],
                "tumor_bam": ["S1.tumor"],
                "normal_bam": ["."],
                "analysis_type": ["tumor_only"],
            }
        ).set_index("sample", drop=False)
        result = get_normal_bam_basenames(df)
        assert result == []


class TestBuildFreebayesParams:
    def test_passthrough(self):
        result = build_freebayes_params("--min-coverage 20 --standard-filters")
        assert result == "--min-coverage 20 --standard-filters"

    def test_empty(self):
        result = build_freebayes_params("")
        assert result == ""
