"""Unit tests for scripts/generate_config.py."""

import os
import sys

import pandas as pd
import pytest

# Add scripts/ to path so generate_config can be imported
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from generate_config import (
    _build_config_yaml,
    _find_gatk_resource_vcfs,
    _find_genome_fastas,
    _infer_output_folder,
    build_samples_dataframe,
    discover_bam_files,
    discover_gatk_resources,
    discover_reference_data,
    extract_patient_id,
    find_best_normal,
    format_size,
    generate_rows,
    generate_sample_name,
    group_and_pair,
    guess_role,
    write_samples_tsv,
)

# ---------------------------------------------------------------------------
# TestDiscoverBamFiles
# ---------------------------------------------------------------------------


class TestDiscoverBamFiles:
    def test_finds_bams_with_extension(self, tmp_path):
        (tmp_path / "sample1.merged.dedup.bqsr.bam").write_bytes(b"x" * 100)
        (tmp_path / "sample2.merged.dedup.bqsr.bam").write_bytes(b"x" * 200)
        result = discover_bam_files(str(tmp_path), ".merged.dedup.bqsr.bam")
        assert len(result) == 2
        basenames = [r["basename"] for r in result]
        assert "sample1" in basenames
        assert "sample2" in basenames

    def test_ignores_wrong_extension(self, tmp_path):
        (tmp_path / "sample1.merged.dedup.bqsr.bam").write_bytes(b"x" * 100)
        (tmp_path / "sample2.sorted.bam").write_bytes(b"x" * 200)
        result = discover_bam_files(str(tmp_path), ".merged.dedup.bqsr.bam")
        assert len(result) == 1
        assert result[0]["basename"] == "sample1"

    def test_detects_bai_index(self, tmp_path):
        (tmp_path / "sample1.merged.dedup.bqsr.bam").write_bytes(b"x" * 100)
        (tmp_path / "sample1.merged.dedup.bqsr.bam.bai").write_bytes(b"x")
        result = discover_bam_files(str(tmp_path), ".merged.dedup.bqsr.bam")
        assert result[0]["has_index"] is True

    def test_no_index_detected(self, tmp_path):
        (tmp_path / "sample1.merged.dedup.bqsr.bam").write_bytes(b"x" * 100)
        result = discover_bam_files(str(tmp_path), ".merged.dedup.bqsr.bam")
        assert result[0]["has_index"] is False

    def test_returns_size(self, tmp_path):
        (tmp_path / "sample1.merged.dedup.bqsr.bam").write_bytes(b"x" * 12345)
        result = discover_bam_files(str(tmp_path), ".merged.dedup.bqsr.bam")
        assert result[0]["size_bytes"] == 12345

    def test_empty_directory(self, tmp_path):
        result = discover_bam_files(str(tmp_path), ".merged.dedup.bqsr.bam")
        assert result == []

    def test_nonexistent_directory(self):
        result = discover_bam_files("/nonexistent/path", ".bam")
        assert result == []


# ---------------------------------------------------------------------------
# TestFormatSize
# ---------------------------------------------------------------------------


class TestFormatSize:
    def test_gigabytes(self):
        assert format_size(45_300_000_000) == "45.3 GB"

    def test_megabytes(self):
        assert format_size(512_000_000) == "512.0 MB"

    def test_kilobytes(self):
        assert format_size(1500) == "1.5 KB"

    def test_bytes(self):
        assert format_size(500) == "500 B"


# ---------------------------------------------------------------------------
# TestInferOutputFolder
# ---------------------------------------------------------------------------


class TestInferOutputFolder:
    def test_replaces_last_component(self):
        assert _infer_output_folder("results/exomes/bqsr") == "results/exomes/variant_calls"

    def test_relative_parent(self):
        assert _infer_output_folder("../results/A5297/bqsr") == "../results/A5297/variant_calls"

    def test_absolute_path(self):
        result = _infer_output_folder("/data/projects/results/cohort1/bqsr")
        assert result == "/data/projects/results/cohort1/variant_calls"

    def test_no_backslashes_on_windows(self):
        result = _infer_output_folder("results\\exomes\\bqsr")
        assert "\\" not in result


# ---------------------------------------------------------------------------
# TestGuessRole
# ---------------------------------------------------------------------------


class TestGuessRole:
    # Tumor patterns
    def test_dash_T(self):
        role, reason = guess_role("APA1-T")
        assert role == "tumor"

    def test_dash_T1(self):
        role, _ = guess_role("APA1-T1")
        assert role == "tumor"

    def test_underscore_tumor(self):
        role, _ = guess_role("SAMPLE_tumor")
        assert role == "tumor"

    def test_dot_FFPE(self):
        role, _ = guess_role("SAMPLE.FFPE")
        assert role == "tumor"

    def test_suffix_met(self):
        role, _ = guess_role("SAMPLE-met")
        assert role == "tumor"

    def test_suffix_primary(self):
        role, _ = guess_role("SAMPLE_primary")
        assert role == "tumor"

    # Normal patterns
    def test_dash_N(self):
        role, _ = guess_role("APA1-N")
        assert role == "normal"

    def test_suffix_N1(self):
        role, _ = guess_role("A5297_DNA_02_STREAM_P1_N1")
        assert role == "normal"

    def test_suffix_normal(self):
        role, _ = guess_role("SAMPLE_normal")
        assert role == "normal"

    def test_suffix_blood(self):
        role, _ = guess_role("SAMPLE-blood")
        assert role == "normal"

    def test_suffix_germline(self):
        role, _ = guess_role("SAMPLE.germline")
        assert role == "normal"

    # Lesion patterns (tumor)
    def test_suffix_L1(self):
        role, reason = guess_role("A5297_DNA_01_STREAM_P1_L1")
        assert role == "tumor"
        assert "lesion" in reason

    def test_suffix_L(self):
        role, _ = guess_role("SAMPLE-L")
        assert role == "tumor"

    def test_suffix_L2(self):
        role, _ = guess_role("SAMPLE_L2")
        assert role == "tumor"

    # Unknown
    def test_unknown_no_match(self):
        role, reason = guess_role("SAMPLE_DNA_01_STREAM_P1")
        assert role == "unknown"
        assert "no pattern" in reason

    # No false positives
    def test_no_false_positive_ント(self):
        """Internal 'N' should not match."""
        role, _ = guess_role("INTERNAL_NAME")
        assert role == "unknown"

    def test_no_false_positive_T_in_middle(self):
        """Internal 'T' should not match."""
        role, _ = guess_role("TOTALLY_FINE")
        assert role == "unknown"

    def test_case_insensitive_tumor(self):
        role, _ = guess_role("SAMPLE-Tumor")
        assert role == "tumor"

    def test_case_insensitive_normal(self):
        role, _ = guess_role("SAMPLE_Normal")
        assert role == "normal"


# ---------------------------------------------------------------------------
# TestExtractPatientId
# ---------------------------------------------------------------------------


class TestExtractPatientId:
    def test_strips_T_suffix(self):
        assert extract_patient_id("APA1-T") == "APA1"

    def test_strips_N_suffix(self):
        assert extract_patient_id("APA1-N") == "APA1"

    def test_strips_tumor_suffix(self):
        assert extract_patient_id("SAMPLE_tumor") == "SAMPLE"

    def test_strips_normal_suffix(self):
        assert extract_patient_id("SAMPLE_normal") == "SAMPLE"

    def test_strips_N1_suffix(self):
        assert extract_patient_id("A5297_DNA_02_STREAM_P1_N1") == "A5297_DNA_02_STREAM_P1"

    def test_strips_L1_suffix(self):
        assert extract_patient_id("A5297_DNA_01_STREAM_P1_L1") == "A5297_DNA_01_STREAM_P1"

    def test_no_suffix_passthrough(self):
        assert extract_patient_id("SAMPLE_DNA_01_STREAM_P1") == "SAMPLE_DNA_01_STREAM_P1"

    def test_strips_blood_suffix(self):
        assert extract_patient_id("SAMPLE-blood") == "SAMPLE"

    def test_strips_FFPE(self):
        assert extract_patient_id("SAMPLE.FFPE") == "SAMPLE"


# ---------------------------------------------------------------------------
# TestFindBestNormal
# ---------------------------------------------------------------------------


class TestFindBestNormal:
    def test_same_patient_priority(self):
        normals = ["APA1-N", "APA2-N"]
        result = find_best_normal("APA1-T", normals, patient_groups={"APA1": {}})
        assert result == "APA1-N"

    def test_string_similarity(self):
        normals = ["A5297_DNA_02_STREAM_P1_N1"]
        result = find_best_normal("A5297_DNA_01_STREAM_P1_L1", normals)
        assert result == "A5297_DNA_02_STREAM_P1_N1"

    def test_no_match_below_threshold(self):
        normals = ["COMPLETELY_DIFFERENT"]
        result = find_best_normal("APA1-T", normals)
        assert result is None

    def test_empty_normals(self):
        result = find_best_normal("APA1-T", [])
        assert result is None

    def test_best_of_multiple(self):
        normals = ["APA1-N", "APA2-N", "APA3-N"]
        result = find_best_normal("APA1-T", normals)
        assert result == "APA1-N"

    def test_patient_groups_preferred(self):
        normals = ["APA2-N", "APA1-N"]
        groups = {
            "APA1": {"tumors": ["APA1-T"], "normals": ["APA1-N"]},
            "APA2": {"tumors": [], "normals": ["APA2-N"]},
        }
        result = find_best_normal("APA1-T", normals, patient_groups=groups)
        assert result == "APA1-N"


# ---------------------------------------------------------------------------
# TestGroupAndPair
# ---------------------------------------------------------------------------


class TestGroupAndPair:
    def test_simple_pair(self):
        entries = [
            {"basename": "APA1-T", "role": "tumor"},
            {"basename": "APA1-N", "role": "normal"},
        ]
        groups, pairings = group_and_pair(entries)
        assert len(pairings) == 1
        assert pairings[0]["tumor"] == "APA1-T"
        assert pairings[0]["normal"] == "APA1-N"

    def test_tumor_only(self):
        entries = [
            {"basename": "APA1-T", "role": "tumor"},
        ]
        _, pairings = group_and_pair(entries)
        assert len(pairings) == 1
        assert pairings[0]["normal"] is None

    def test_multiple_patients(self):
        entries = [
            {"basename": "APA1-T", "role": "tumor"},
            {"basename": "APA1-N", "role": "normal"},
            {"basename": "APA2-T", "role": "tumor"},
            {"basename": "APA2-N", "role": "normal"},
        ]
        _, pairings = group_and_pair(entries)
        assert len(pairings) == 2

        # APA1-T should be paired with APA1-N
        apa1 = next(p for p in pairings if p["tumor"] == "APA1-T")
        assert apa1["normal"] == "APA1-N"

        # APA2-T should be paired with APA2-N
        apa2 = next(p for p in pairings if p["tumor"] == "APA2-T")
        assert apa2["normal"] == "APA2-N"

    def test_skipped_entries_excluded(self):
        entries = [
            {"basename": "APA1-T", "role": "tumor"},
            {"basename": "SKIP_ME", "role": "skip"},
        ]
        _, pairings = group_and_pair(entries)
        assert len(pairings) == 1

    def test_complex_names_pair(self):
        entries = [
            {"basename": "A5297_DNA_01_STREAM_P1_L1", "role": "tumor"},
            {"basename": "A5297_DNA_02_STREAM_P1_N1", "role": "normal"},
        ]
        _, pairings = group_and_pair(entries)
        assert len(pairings) == 1
        # These should pair via string similarity
        assert pairings[0]["normal"] == "A5297_DNA_02_STREAM_P1_N1"

    def test_patient_groups_populated(self):
        entries = [
            {"basename": "APA1-T", "role": "tumor"},
            {"basename": "APA1-N", "role": "normal"},
        ]
        groups, _ = group_and_pair(entries)
        assert "APA1" in groups
        assert "APA1-T" in groups["APA1"]["tumors"]
        assert "APA1-N" in groups["APA1"]["normals"]

    def test_multiple_tumors_one_normal(self):
        entries = [
            {"basename": "APA1-T1", "role": "tumor"},
            {"basename": "APA1-T2", "role": "tumor"},
            {"basename": "APA1-N", "role": "normal"},
        ]
        _, pairings = group_and_pair(entries)
        assert len(pairings) == 2
        # Same-patient normal is shared across all tumors from that patient
        normals_used = [p["normal"] for p in pairings if p["normal"] is not None]
        assert len(normals_used) == 2
        assert all(n == "APA1-N" for n in normals_used)

    def test_no_tumors(self):
        entries = [
            {"basename": "APA1-N", "role": "normal"},
        ]
        _, pairings = group_and_pair(entries)
        assert len(pairings) == 0


# ---------------------------------------------------------------------------
# TestGenerateRows
# ---------------------------------------------------------------------------


class TestGenerateRows:
    @pytest.fixture
    def simple_pairings(self):
        return [
            {"tumor": "APA1-T", "normal": "APA1-N", "patient_id": "APA1"},
        ]

    @pytest.fixture
    def simple_entries(self):
        return [
            {"basename": "APA1-T", "role": "tumor"},
            {"basename": "APA1-N", "role": "normal"},
        ]

    def test_both_mode(self, simple_pairings, simple_entries):
        rows = generate_rows(simple_pairings, simple_entries, "both")
        assert len(rows) == 2
        types = [r["analysis_type"] for r in rows]
        assert "tumor_only" in types
        assert "tumor_normal" in types

    def test_tumor_only_mode(self, simple_pairings, simple_entries):
        rows = generate_rows(simple_pairings, simple_entries, "tumor_only")
        assert len(rows) == 1
        assert rows[0]["analysis_type"] == "tumor_only"
        assert rows[0]["normal_bam"] == "."

    def test_tumor_normal_mode(self, simple_pairings, simple_entries):
        rows = generate_rows(simple_pairings, simple_entries, "tumor_normal")
        assert len(rows) == 1
        assert rows[0]["analysis_type"] == "tumor_normal"
        assert rows[0]["normal_bam"] == "APA1-N"

    def test_germline_mode(self, simple_pairings, simple_entries):
        rows = generate_rows(simple_pairings, simple_entries, "germline")
        assert len(rows) == 2
        for row in rows:
            assert row["analysis_type"] == "germline"

    def test_unpaired_tumor_in_tn_mode(self):
        pairings = [{"tumor": "APA1-T", "normal": None, "patient_id": "APA1"}]
        entries = [{"basename": "APA1-T", "role": "tumor"}]
        rows = generate_rows(pairings, entries, "tumor_normal")
        assert len(rows) == 0  # No TN row without normal

    def test_mixed_patients(self):
        pairings = [
            {"tumor": "APA1-T", "normal": "APA1-N", "patient_id": "APA1"},
            {"tumor": "APA2-T", "normal": None, "patient_id": "APA2"},
        ]
        entries = [
            {"basename": "APA1-T", "role": "tumor"},
            {"basename": "APA1-N", "role": "normal"},
            {"basename": "APA2-T", "role": "tumor"},
        ]
        rows = generate_rows(pairings, entries, "both")
        # APA1: To + TN, APA2: To only
        assert len(rows) == 3


# ---------------------------------------------------------------------------
# TestGenerateSampleName
# ---------------------------------------------------------------------------


class TestGenerateSampleName:
    def test_tumor_normal_suffix(self):
        assert generate_sample_name("APA1", "tumor_normal", 0, 1) == "APA1_TN"

    def test_tumor_only_suffix(self):
        assert generate_sample_name("APA1", "tumor_only", 0, 1) == "APA1_To"

    def test_germline_suffix(self):
        assert generate_sample_name("APA1", "germline", 0, 1) == "APA1_G"

    def test_disambiguation_index(self):
        assert generate_sample_name("APA1", "tumor_normal", 0, 2) == "APA1_TN_1"
        assert generate_sample_name("APA1", "tumor_normal", 1, 2) == "APA1_TN_2"

    def test_single_no_index(self):
        name = generate_sample_name("APA1", "tumor_only", 0, 1)
        assert "_1" not in name


# ---------------------------------------------------------------------------
# TestGatkResourceDiscovery
# ---------------------------------------------------------------------------


class TestGatkResourceDiscovery:
    def test_finds_pon(self, tmp_path):
        (tmp_path / "1000g_pon.hg38.vcf.gz").write_bytes(b"x")
        result = _find_gatk_resource_vcfs(tmp_path, "panel_of_normals", ["*pon*.vcf*"])
        assert len(result) == 1
        assert "1000g_pon" in result[0]

    def test_excludes_tbi(self, tmp_path):
        (tmp_path / "1000g_pon.hg38.vcf.gz").write_bytes(b"x")
        (tmp_path / "1000g_pon.hg38.vcf.gz.tbi").write_bytes(b"x")
        result = _find_gatk_resource_vcfs(tmp_path, "panel_of_normals", ["*pon*.vcf*"])
        assert len(result) == 1
        assert not result[0].endswith(".tbi")

    def test_af_gnomad_excludes_common(self, tmp_path):
        (tmp_path / "af-only-gnomad.hg38.vcf.gz").write_bytes(b"x")
        (tmp_path / "af-only-gnomad.hg38.common_biallelic.vcf.gz").write_bytes(b"x")
        result = _find_gatk_resource_vcfs(tmp_path, "af_only_gnomad", ["*af-only-gnomad*.vcf*"])
        assert len(result) == 1
        # Check the filename only, not the full path (temp dirs may contain "common")
        from pathlib import Path

        assert "common" not in Path(result[0]).name

    def test_common_biallelic_found(self, tmp_path):
        (tmp_path / "af-only-gnomad.hg38.common_biallelic.vcf.gz").write_bytes(b"x")
        result = _find_gatk_resource_vcfs(
            tmp_path, "common_biallelic_gnomad", ["*common_biallelic*.vcf*"]
        )
        assert len(result) == 1

    def test_empty_dir(self, tmp_path):
        result = _find_gatk_resource_vcfs(tmp_path, "panel_of_normals", ["*pon*.vcf*"])
        assert result == []

    def test_nonexistent_dir(self):
        from pathlib import Path

        result = _find_gatk_resource_vcfs(Path("/nonexistent"), "panel_of_normals", ["*pon*.vcf*"])
        assert result == []


class TestDiscoverGatkResources:
    def test_finds_all_resources(self, tmp_path):
        (tmp_path / "1000g_pon.hg38.vcf.gz").write_bytes(b"x")
        (tmp_path / "af-only-gnomad.hg38.vcf.gz").write_bytes(b"x")
        (tmp_path / "af-only-gnomad.hg38.common_biallelic.vcf.gz").write_bytes(b"x")
        result = discover_gatk_resources(tmp_path)
        assert result["panel_of_normals"] != ""
        assert result["af_only_gnomad"] != ""
        assert result["common_biallelic_gnomad"] != ""

    def test_missing_resources_logged(self, tmp_path):
        result = discover_gatk_resources(tmp_path)
        assert any("not found" in line for line in result["search_log"])


# ---------------------------------------------------------------------------
# TestReferenceDiscovery
# ---------------------------------------------------------------------------


class TestReferenceDiscovery:
    def test_finds_fasta(self, tmp_path):
        (tmp_path / "genome.fna").write_bytes(b"x")
        (tmp_path / "genome.fna.fai").write_bytes(b"x")
        result = discover_reference_data(tmp_path)
        assert result["genome"] != ""

    def test_detects_grch38(self, tmp_path):
        (tmp_path / "GRCh38.fna").write_bytes(b"x")
        result = discover_reference_data(tmp_path)
        assert result["build"] == "GRCh38"

    def test_detects_grch37(self, tmp_path):
        (tmp_path / "hs37d5.fasta").write_bytes(b"x")
        result = discover_reference_data(tmp_path)
        assert result["build"] == "GRCh37"

    def test_no_genome_found(self, tmp_path):
        result = discover_reference_data(tmp_path)
        assert result["genome"] == ""
        assert any("No reference genome" in line for line in result["search_log"])


class TestFindGenomeFastas:
    def test_finds_multiple_extensions(self, tmp_path):
        (tmp_path / "genome.fna").write_bytes(b"x")
        (tmp_path / "genome.fa").write_bytes(b"x")
        result = _find_genome_fastas(tmp_path)
        assert len(result) == 2

    def test_checks_fai(self, tmp_path):
        (tmp_path / "genome.fna").write_bytes(b"x")
        (tmp_path / "genome.fna.fai").write_bytes(b"x")
        result = _find_genome_fastas(tmp_path)
        assert result[0]["has_fai"] is True

    def test_skips_hidden(self, tmp_path):
        (tmp_path / ".hidden.fna").write_bytes(b"x")
        result = _find_genome_fastas(tmp_path)
        assert len(result) == 0

    def test_nonexistent_dir(self):
        from pathlib import Path

        result = _find_genome_fastas(Path("/nonexistent"))
        assert result == []


# ---------------------------------------------------------------------------
# TestBuildConfigYaml
# ---------------------------------------------------------------------------


class TestBuildConfigYaml:
    @pytest.fixture
    def ref_data(self):
        return {"genome": "/ref/GRCh38.fna", "build": "GRCh38"}

    @pytest.fixture
    def gatk_resources(self):
        return {
            "panel_of_normals": "/res/pon.vcf.gz",
            "af_only_gnomad": "/res/gnomad.vcf.gz",
            "common_biallelic_gnomad": "/res/common.vcf.gz",
        }

    def test_contains_caller(self, ref_data, gatk_resources):
        content = _build_config_yaml(
            "mutect2",
            ref_data,
            gatk_resources,
            "/bams",
            "/output",
            "config/samples.tsv",
            ".bam",
            "chromosome",
        )
        assert 'caller: "mutect2"' in content

    def test_contains_ref(self, ref_data, gatk_resources):
        content = _build_config_yaml(
            "mutect2",
            ref_data,
            gatk_resources,
            "/bams",
            "/output",
            "config/samples.tsv",
            ".bam",
            "chromosome",
        )
        assert "genome:" in content
        assert "GRCh38" in content

    def test_contains_gatk_resources(self, ref_data, gatk_resources):
        content = _build_config_yaml(
            "mutect2",
            ref_data,
            gatk_resources,
            "/bams",
            "/output",
            "config/samples.tsv",
            ".bam",
            "chromosome",
        )
        assert "panel_of_normals:" in content
        assert "af_only_gnomad:" in content
        assert "common_biallelic_gnomad:" in content

    def test_contains_scatter(self, ref_data, gatk_resources):
        content = _build_config_yaml(
            "mutect2",
            ref_data,
            gatk_resources,
            "/bams",
            "/output",
            "config/samples.tsv",
            ".bam",
            "chromosome",
        )
        assert 'mode: "chromosome"' in content
        assert "chr1" in content
        assert "chrY" in content

    def test_edit_me_placeholders(self):
        content = _build_config_yaml(
            "mutect2",
            {"genome": "", "build": "GRCh38"},
            {"panel_of_normals": "", "af_only_gnomad": "", "common_biallelic_gnomad": ""},
            "/bams",
            "/output",
            "config/samples.tsv",
            ".bam",
            "chromosome",
        )
        assert "EDIT_ME" in content


# ---------------------------------------------------------------------------
# TestBuildSamplesDataframe
# ---------------------------------------------------------------------------


class TestBuildSamplesDataframe:
    def test_correct_columns(self):
        rows = [
            {
                "sample": "S1_TN",
                "tumor_bam": "S1-T",
                "normal_bam": "S1-N",
                "analysis_type": "tumor_normal",
                "patient_id": "S1",
            },
        ]
        df = build_samples_dataframe(rows)
        assert list(df.columns) == ["sample", "tumor_bam", "normal_bam", "analysis_type"]

    def test_sorted_by_sample(self):
        rows = [
            {
                "sample": "B_TN",
                "tumor_bam": "B-T",
                "normal_bam": "B-N",
                "analysis_type": "tumor_normal",
                "patient_id": "B",
            },
            {
                "sample": "A_To",
                "tumor_bam": "A-T",
                "normal_bam": ".",
                "analysis_type": "tumor_only",
                "patient_id": "A",
            },
        ]
        df = build_samples_dataframe(rows)
        assert df.iloc[0]["sample"] == "A_To"
        assert df.iloc[1]["sample"] == "B_TN"

    def test_dot_sentinel_preserved(self):
        rows = [
            {
                "sample": "S1_To",
                "tumor_bam": "S1-T",
                "normal_bam": ".",
                "analysis_type": "tumor_only",
                "patient_id": "S1",
            },
        ]
        df = build_samples_dataframe(rows)
        assert df.iloc[0]["normal_bam"] == "."

    def test_empty_rows(self):
        df = build_samples_dataframe([])
        assert len(df) == 0
        assert list(df.columns) == ["sample", "tumor_bam", "normal_bam", "analysis_type"]


# ---------------------------------------------------------------------------
# Integration-style tests: end-to-end row generation
# ---------------------------------------------------------------------------


class TestEndToEndGeneration:
    """Test the full pipeline from BAM entries to sample rows."""

    def test_simple_pair_both_mode(self):
        entries = [
            {"basename": "APA1-T", "role": "tumor"},
            {"basename": "APA1-N", "role": "normal"},
        ]
        _, pairings = group_and_pair(entries)
        rows = generate_rows(pairings, entries, "both")
        df = build_samples_dataframe(rows)

        assert len(df) == 2
        to_row = df[df["analysis_type"] == "tumor_only"].iloc[0]
        assert to_row["tumor_bam"] == "APA1-T"
        assert to_row["normal_bam"] == "."

        tn_row = df[df["analysis_type"] == "tumor_normal"].iloc[0]
        assert tn_row["tumor_bam"] == "APA1-T"
        assert tn_row["normal_bam"] == "APA1-N"

    def test_complex_names(self):
        entries = [
            {"basename": "A5297_DNA_01_STREAM_P1_L1", "role": "tumor"},
            {"basename": "A5297_DNA_02_STREAM_P1_N1", "role": "normal"},
        ]
        _, pairings = group_and_pair(entries)
        rows = generate_rows(pairings, entries, "both")
        df = build_samples_dataframe(rows)

        assert len(df) == 2
        tn_rows = df[df["analysis_type"] == "tumor_normal"]
        assert len(tn_rows) == 1
        assert tn_rows.iloc[0]["normal_bam"] == "A5297_DNA_02_STREAM_P1_N1"

    def test_germline_mode_all_bams(self):
        entries = [
            {"basename": "APA1-T", "role": "tumor"},
            {"basename": "APA1-N", "role": "normal"},
            {"basename": "APA2-T", "role": "tumor"},
        ]
        _, pairings = group_and_pair(entries)
        rows = generate_rows(pairings, entries, "germline")
        df = build_samples_dataframe(rows)

        assert len(df) == 3
        assert all(df["analysis_type"] == "germline")

    def test_multiple_patients_both(self):
        entries = [
            {"basename": "APA1-T", "role": "tumor"},
            {"basename": "APA1-N", "role": "normal"},
            {"basename": "APA2-T", "role": "tumor"},
            {"basename": "APA2-N", "role": "normal"},
        ]
        _, pairings = group_and_pair(entries)
        rows = generate_rows(pairings, entries, "both")
        df = build_samples_dataframe(rows)

        # 2 patients x (1 To + 1 TN) = 4 rows
        assert len(df) == 4
        assert len(df[df["analysis_type"] == "tumor_only"]) == 2
        assert len(df[df["analysis_type"] == "tumor_normal"]) == 2

    def test_tumor_no_normal(self):
        entries = [
            {"basename": "APA1-T", "role": "tumor"},
        ]
        _, pairings = group_and_pair(entries)
        rows = generate_rows(pairings, entries, "both")
        df = build_samples_dataframe(rows)

        # Only tumor_only row (no normal for TN)
        assert len(df) == 1
        assert df.iloc[0]["analysis_type"] == "tumor_only"

    def test_skipped_bams_excluded(self):
        entries = [
            {"basename": "APA1-T", "role": "tumor"},
            {"basename": "APA1-N", "role": "normal"},
            {"basename": "SKIP_ME", "role": "skip"},
        ]
        _, pairings = group_and_pair(entries)
        rows = generate_rows(pairings, entries, "both")
        df = build_samples_dataframe(rows)

        bam_names = list(df["tumor_bam"]) + list(df["normal_bam"])
        assert "SKIP_ME" not in bam_names


# ---------------------------------------------------------------------------
# Test write/output functions
# ---------------------------------------------------------------------------


class TestWriteSamplesTsv:
    def test_dry_run_no_file(self, tmp_path, capsys):
        df = pd.DataFrame(
            {
                "sample": ["S1"],
                "tumor_bam": ["T1"],
                "normal_bam": ["."],
                "analysis_type": ["tumor_only"],
            }
        )
        out = tmp_path / "samples.tsv"
        write_samples_tsv(df, out, dry_run=True)
        assert not out.exists()
        assert "dry-run" in capsys.readouterr().out

    def test_writes_file(self, tmp_path):
        df = pd.DataFrame(
            {
                "sample": ["S1"],
                "tumor_bam": ["T1"],
                "normal_bam": ["."],
                "analysis_type": ["tumor_only"],
            }
        )
        out = tmp_path / "samples.tsv"
        write_samples_tsv(df, out, force=True)
        assert out.exists()
        content = out.read_text()
        assert "S1" in content
        assert "tumor_only" in content

    def test_creates_parent_dirs(self, tmp_path):
        df = pd.DataFrame(
            {
                "sample": ["S1"],
                "tumor_bam": ["T1"],
                "normal_bam": ["."],
                "analysis_type": ["tumor_only"],
            }
        )
        out = tmp_path / "subdir" / "samples.tsv"
        write_samples_tsv(df, out, force=True)
        assert out.exists()

    def test_tsv_format(self, tmp_path):
        df = pd.DataFrame(
            {
                "sample": ["S1", "S2"],
                "tumor_bam": ["T1", "T2"],
                "normal_bam": [".", "N2"],
                "analysis_type": ["tumor_only", "tumor_normal"],
            }
        )
        out = tmp_path / "samples.tsv"
        write_samples_tsv(df, out, force=True)
        # Read back and verify tab-separated
        loaded = pd.read_csv(out, sep="\t")
        assert list(loaded.columns) == ["sample", "tumor_bam", "normal_bam", "analysis_type"]
        assert len(loaded) == 2
