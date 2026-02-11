"""Unit tests for workflow/rules/helpers.py."""

import sys
import os

# Add workflow/rules to the path so helpers.py can be imported directly
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "rules"))

from helpers import (
    build_freebayes_params,
    get_germline_samples,
    get_java_opts,
    get_mutect2_samples,
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


class TestBuildFreebayesParams:
    def test_passthrough(self):
        result = build_freebayes_params("--min-coverage 20 --standard-filters")
        assert result == "--min-coverage 20 --standard-filters"

    def test_empty(self):
        result = build_freebayes_params("")
        assert result == ""
