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
        assert "--disable-adaptive-pruning true" in result

    def test_disabled(self):
        result = build_mutect2_extra_args(
            genotype_germline_sites=False,
            genotype_pon_sites=False,
            disable_adaptive_pruning=False,
        )
        assert result == ""

    def test_disable_adaptive_pruning_on(self):
        result = build_mutect2_extra_args(
            genotype_germline_sites=False,
            genotype_pon_sites=False,
            disable_adaptive_pruning=True,
        )
        assert result == "--disable-adaptive-pruning true"

    def test_disable_adaptive_pruning_off(self):
        result = build_mutect2_extra_args(
            genotype_germline_sites=False,
            genotype_pon_sites=False,
            disable_adaptive_pruning=False,
        )
        assert "--disable-adaptive-pruning" not in result

    def test_annotations(self):
        result = build_mutect2_extra_args(
            genotype_germline_sites=False,
            genotype_pon_sites=False,
            disable_adaptive_pruning=False,
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
            disable_adaptive_pruning=True,
            annotations=["OrientationBiasReadCounts"],
            annotation_groups=["AS_StandardAnnotation"],
            extra="--custom-flag",
        )
        parts = result.split()
        # Verify ordering: germline, pon, pruning, annotations, annotation-groups, extra
        germline_idx = parts.index("--genotype-germline-sites")
        pon_idx = parts.index("--genotype-pon-sites")
        pruning_idx = parts.index("--disable-adaptive-pruning")
        ann_idx = parts.index("--annotation")
        grp_idx = parts.index("--annotation-group")
        custom_idx = parts.index("--custom-flag")
        assert germline_idx < pon_idx < pruning_idx < ann_idx < grp_idx < custom_idx


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


class _DefaultNamespace(dict):
    def __getattr__(self, key):
        return self.get(key, "")


def _extract_rule_shell(smk_relpath: str, rule_name: str) -> str:
    import re
    from pathlib import Path

    text = (Path(__file__).parent.parent / smk_relpath).read_text()
    pattern = rf"rule {rule_name}:.*?shell:\s*r?\"\"\"(.*?)\"\"\""
    match = re.search(pattern, text, re.S)
    assert match, f"shell block not found for {rule_name}"
    return match.group(1)


class TestMutect2CallSmResolution:
    @staticmethod
    def _run_mutect2_shell(tmp_path, tumor_headers, normal_header=None):
        import subprocess

        bin_dir = tmp_path / "bin"
        bin_dir.mkdir(exist_ok=True)
        gatk_args_file = tmp_path / "gatk_args.txt"

        samtools_bin = bin_dir / "samtools"
        samtools_bin.write_text('#!/usr/bin/env bash\nshift 2\ncat "$1"\n')
        samtools_bin.chmod(0o755)

        gatk_bin = bin_dir / "gatk"
        gatk_bin.write_text(f'#!/usr/bin/env bash\nprintf "%s\\n" "$@" > "{gatk_args_file}"\n')
        gatk_bin.chmod(0o755)

        if isinstance(tumor_headers, str):
            tumor_headers = {"tumor_file.bam": tumor_headers}
        tumor_paths = []
        for fname, hdr in tumor_headers.items():
            p = tmp_path / fname
            p.write_text(hdr)
            tumor_paths.append(str(p))

        normal_path = ""
        if normal_header is not None:
            np = tmp_path / "normal_file.bam"
            np.write_text(normal_header)
            normal_path = str(np)

        log_path = tmp_path / "mutect2.log"
        shell_tpl = _extract_rule_shell("workflow/rules/mutect2.smk", "mutect2_call")
        script = shell_tpl.format(
            input=_DefaultNamespace(
                reference="ref.fa",
                tumor_bam=" ".join(tumor_paths),
                normal_bam=normal_path,
                gnomad="gnomad.vcf.gz",
                pon="pon.vcf.gz",
            ),
            output=_DefaultNamespace(f1r2="f1r2.tar.gz", vcf="out.vcf.gz"),
            params=_DefaultNamespace(
                java_opts="-Xmx4g",
                tumor_input_arg=" ".join(f"-I {p}" for p in tumor_paths),
                normal_input_arg=f"-I {normal_path}" if normal_path else "",
                pon_arg="--panel-of-normals pon.vcf.gz",
                extra="",
                interval_arg="-L chr1",
            ),
            threads=4,
            log=str(log_path),
        )
        env = {**os.environ, "PATH": f"{bin_dir}:{os.environ['PATH']}"}
        res = subprocess.run(
            ["bash", "-euo", "pipefail", "-c", script],
            env=env,
            capture_output=True,
            text=True,
        )
        gatk_args = gatk_args_file.read_text().splitlines() if gatk_args_file.exists() else []
        return res, gatk_args

    def test_resolves_normal_sm_when_filename_differs(self, tmp_path):
        res, gatk_args = self._run_mutect2_shell(
            tmp_path,
            tumor_headers={"APA1_T.bam": "@RG\tID:1\tSM:TRUE_TUMOR_SM\n"},
            normal_header="@RG\tID:2\tSM:TRUE_NORMAL_SM\n",
        )
        assert res.returncode == 0
        idx = gatk_args.index("-normal")
        assert gatk_args[idx + 1] == "TRUE_NORMAL_SM"

    def test_fails_when_tumor_and_normal_sm_collide(self, tmp_path):
        res, _ = self._run_mutect2_shell(
            tmp_path,
            tumor_headers={"APA1_T.bam": "@RG\tID:1\tSM:SAME_SM\n"},
            normal_header="@RG\tID:2\tSM:SAME_SM\n",
        )
        assert res.returncode != 0

    def test_fails_on_missing_or_multisample_sm(self, tmp_path):
        res_empty, _ = self._run_mutect2_shell(
            tmp_path,
            tumor_headers={"APA1_T.bam": "@RG\tID:1\n"},
        )
        assert res_empty.returncode != 0

        res_multi, _ = self._run_mutect2_shell(
            tmp_path,
            tumor_headers={"APA1_T.bam": "@RG\tID:1\tSM:S1\n@RG\tID:2\tSM:S2\n"},
        )
        assert res_multi.returncode != 0

