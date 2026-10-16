"""Unit tests for workflow/rules/helpers.py."""

import os
import sys

import pytest

# Add workflow and workflow/rules to the path so helpers.py and rules.helpers can be imported
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "rules"))

from helpers import (
    build_freebayes_params,
    build_interval_arg,
    build_mutect2_extra_args,
    build_pon_interval_arg,
    get_germline_samples,
    get_java_opts,
    get_mutect2_samples,
    get_normal_bam_basenames,
    get_scatter_units,
    get_tumor_bam_basenames,
    get_unique_bam_basenames,
    has_matched_normal,
    validate_samples_df,
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

    def test_tumor_only_ignores_normal_bam(self):
        import pandas as pd

        df = pd.DataFrame(
            {
                "sample": ["S1"],
                "tumor_bam": ["S1.T"],
                "normal_bam": ["S1.N"],
                "analysis_type": ["tumor_only"],
            }
        ).set_index("sample", drop=False)
        assert has_matched_normal(df, "S1") is False


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

    def test_chromosome_mode_filters_by_target_bed_preserving_order(self):
        from unittest.mock import mock_open, patch

        bed_content = "# header\nchr1\t100\t200\nchr2\t300\t400\n"
        with patch("builtins.open", mock_open(read_data=bed_content)):
            result = get_scatter_units(
                "chromosome",
                ["chr2", "chr1", "chrM"],
                target_bed="exome.bed",
            )
        assert result == ["chr2", "chr1"]


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

    def test_fails_when_multiple_tumor_bams_share_duplicate_sm(self, tmp_path):
        res, _ = self._run_mutect2_shell(
            tmp_path,
            tumor_headers={
                "T1.bam": "@RG\tID:1\tSM:SAME_TUMOR_SM\n",
                "T2.bam": "@RG\tID:2\tSM:SAME_TUMOR_SM\n",
            },
            normal_header="@RG\tID:3\tSM:NORMAL_SM\n",
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


class TestValidateSamplesDf:
    def test_valid_samples_df(self, samples_df):
        validate_samples_df(samples_df)

    def test_duplicate_sample_rejected(self):
        import pandas as pd
        import pytest

        df = pd.DataFrame(
            {
                "sample": ["S1", "S1"],
                "tumor_bam": ["S1.T1", "S1.T2"],
                "normal_bam": [".", "."],
                "analysis_type": ["tumor_only", "tumor_only"],
            }
        ).set_index("sample", drop=False)
        with pytest.raises(ValueError, match="Duplicate sample"):
            validate_samples_df(df)

    def test_numeric_ids_read_as_str_and_non_str_rejected(self, tmp_path):
        import pandas as pd
        import pytest

        tsv = tmp_path / "samples.tsv"
        tsv.write_text(
            "sample\ttumor_bam\tnormal_bam\tanalysis_type\n1001\t2001\t3001\ttumor_normal\n"
        )
        df_str = pd.read_table(tsv, dtype=str, keep_default_na=False).set_index(
            "sample", drop=False, verify_integrity=True
        )
        validate_samples_df(df_str)
        assert df_str.loc["1001", "tumor_bam"] == "2001"
        assert has_matched_normal(df_str, "1001") is True

        df_int = pd.read_table(tsv).set_index("sample", drop=False)
        with pytest.raises(ValueError, match="strings"):
            validate_samples_df(df_int)

    def test_tumor_normal_missing_normal_rejected(self):
        import pandas as pd
        import pytest

        df = pd.DataFrame(
            {
                "sample": ["S1"],
                "tumor_bam": ["S1.T"],
                "normal_bam": ["."],
                "analysis_type": ["tumor_normal"],
            }
        ).set_index("sample", drop=False)
        with pytest.raises(ValueError, match="missing normal_bam"):
            validate_samples_df(df)

    def test_tumor_only_with_normal_rejected(self):
        import pandas as pd
        import pytest

        df = pd.DataFrame(
            {
                "sample": ["S1"],
                "tumor_bam": ["S1.T"],
                "normal_bam": ["S1.N"],
                "analysis_type": ["tumor_only"],
            }
        ).set_index("sample", drop=False)
        with pytest.raises(ValueError, match="tumor_only"):
            validate_samples_df(df)


class TestScatterIntervalsAndTargetBed:
    def test_build_interval_arg_with_target_bed(self):
        assert (
            build_interval_arg("chromosome", "chr1", target_bed="targets.bed")
            == "-L chr1 -L targets.bed --interval-set-rule INTERSECTION"
        )
        assert build_interval_arg("chromosome", "chr1") == "-L chr1"
        assert (
            build_interval_arg("interval", "0000-scattered", "/iv", "targets.bed")
            == "-L /iv/0000-scattered.interval_list -L targets.bed --interval-set-rule INTERSECTION"
        )
        assert build_interval_arg("none", "all", target_bed="targets.bed") == "-L targets.bed"

    def test_split_intervals_dag_and_bed_conversion(self, tmp_path, config_dict, samples_df):
        import subprocess
        from pathlib import Path

        cfg = dict(config_dict)
        cfg["scatter"] = {"mode": "interval", "count": 2, "target_bed": "targets.bed"}
        ns = {"config": cfg, "samples_df": samples_df}
        common_code = (Path(__file__).parent.parent / "workflow/rules/common.smk").read_text()
        exec(common_code, ns)

        wc = _DefaultNamespace(scatter_unit="0000-scattered")
        assert ns["get_interval_input"](wc) == [
            "/output/intervals/0000-scattered.interval_list",
            "targets.bed",
        ]
        assert ns["get_freebayes_interval_input"](wc) == [
            "/output/intervals/0000-scattered.bed",
        ]

        target_bed_path = tmp_path / "targets.bed"
        target_bed_path.write_text("chr1\t100\t200\nchr2\t300\t400\n")
        cfg_chr = dict(config_dict)
        cfg_chr["scatter"] = {
            "mode": "chromosome",
            "chromosomes": ["chr1", "chr2", "chrM"],
            "target_bed": str(target_bed_path),
        }
        ns_chr = {"config": cfg_chr, "samples_df": samples_df}
        exec(common_code, ns_chr)
        assert ns_chr["SCATTER_UNITS"] == ["chr1", "chr2"]
        wc_chr = _DefaultNamespace(scatter_unit="chr1")
        assert ns_chr["get_freebayes_interval_input"](wc_chr) == [str(target_bed_path)]
        region_arg = ns_chr["get_freebayes_region_arg"](wc_chr)
        assert "-r " not in region_arg
        assert "--targets <(awk" in region_arg

        captured_targets = tmp_path / "captured_targets.bed"
        res = subprocess.run(
            [
                "bash",
                "-euo",
                "pipefail",
                "-c",
                f'capture() {{ cat "$2" > "{captured_targets}"; }}; capture {region_arg}',
            ],
            capture_output=True,
            text=True,
        )
        assert res.returncode == 0
        assert captured_targets.read_text() == "chr1\t100\t200\n"

        split_shell = _extract_rule_shell("workflow/rules/scatter.smk", "split_intervals")
        assert 'gatk IntervalListToBed -I "$interval" -O "$bed"' in split_shell
        assert "sed '" not in split_shell

        mutect2_text = (Path(__file__).parent.parent / "workflow/rules/mutect2.smk").read_text()
        assert "intervals=get_interval_input" in mutect2_text
        assert "--interval-set-rule INTERSECTION" in mutect2_text
        assert "360 * attempt" in mutect2_text

        freebayes_text = (
            Path(__file__).parent.parent / "workflow/rules/freebayes.smk"
        ).read_text()
        assert "intervals=get_freebayes_interval_input" in freebayes_text


class TestProfilesAndResources:
    def test_bih_and_default_profiles_no_mem_conflict_or_static_override(self):
        from pathlib import Path

        import yaml

        root = Path(__file__).parent.parent
        bih = yaml.safe_load((root / "profiles/bih/config.v8+.yaml").read_text())
        charite = yaml.safe_load((root / "profiles/charite/config.v8+.yaml").read_text())
        default = yaml.safe_load((root / "profiles/default/config.v8+.yaml").read_text())

        assert bih["default-resources"]["slurm_account"] == "hpc-ag-scholl"
        assert "mem_mb_per_cpu" not in bih["default-resources"]
        assert "mem_mb_per_cpu" not in charite["default-resources"]
        assert default["default-resources"]["mem_mb"] == 4000
        assert "set-resources" not in default

    def test_rule_runtimes_within_7day_partition_limit_and_hmm_threads(self, tmp_path):
        import re
        from pathlib import Path

        root = Path(__file__).parent.parent
        for smk in [
            "workflow/rules/mutect2.smk",
            "workflow/rules/freebayes.smk",
            "workflow/rules/purecn.smk",
        ]:
            text = (root / smk).read_text()
            for m in re.finditer(
                r"retries:\s*(\d+).*?runtime=lambda wildcards, attempt:\s*(\d+)\s*\*\s*attempt",
                text,
                re.S,
            ):
                retries = int(m.group(1))
                base_runtime = int(m.group(2))
                assert base_runtime * (retries + 1) <= 10080

        res, gatk_args = TestMutect2CallSmResolution._run_mutect2_shell(
            tmp_path,
            tumor_headers={"APA1_T.bam": "@RG\tID:1\tSM:T1\n"},
        )
        assert res.returncode == 0
        idx = gatk_args.index("--native-pair-hmm-threads")
        assert gatk_args[idx + 1] == "4"

    def test_profile_runtime_has_explicit_time_unit(self):
        import re
        from pathlib import Path

        import yaml

        root = Path(__file__).parent.parent
        profiles = sorted((root / "profiles").glob("*/config*.yaml"))
        assert profiles
        checked = 0
        for path in profiles:
            data = yaml.safe_load(path.read_text()) or {}
            default_res = data.get("default-resources") or {}
            if "runtime" in default_res:
                val = default_res["runtime"]
                assert isinstance(val, str) and re.fullmatch(r"\d+[smhd]", val), (
                    f"{path}: default-resources.runtime must have an explicit unit, got {val!r}"
                )
                checked += 1
            set_res = data.get("set-resources") or {}
            for rule, resources in set_res.items():
                if isinstance(resources, dict) and "runtime" in resources:
                    val = resources["runtime"]
                    assert isinstance(val, str) and re.fullmatch(r"\d+[smhd]", val), (
                        f"{path}: set-resources.{rule}.runtime must have an explicit unit, got {val!r}"
                    )
                    checked += 1
        assert checked >= 3


class TestFilterMutectCallsNormAndExtra:
    def test_extra_passthrough_and_bcftools_norm_without_atomize(
        self, tmp_path, config_dict, samples_df
    ):
        import shutil
        import subprocess
        from pathlib import Path

        import yaml

        cfg = dict(config_dict)
        cfg["params"] = {
            **config_dict["params"],
            "filter_mutect_calls": {"extra": "--min-median-mapping-quality 25"},
        }
        ns = {"config": cfg, "samples_df": samples_df}
        root = Path(__file__).parent.parent
        exec((root / "workflow/rules/common.smk").read_text(), ns)
        assert ns["FILTER_MUTECT_CALLS_EXTRA"] == "--min-median-mapping-quality 25"

        gatk_env = yaml.safe_load((root / "workflow/envs/gatk.yaml").read_text())
        assert "bcftools=1.21" in gatk_env["dependencies"]
        assert "htslib=1.21" in gatk_env["dependencies"]

        real_bcftools = shutil.which("bcftools")
        if not real_bcftools:
            for candidate in (
                "/home/bernt-popp/miniforge3/envs/vntyper/bin/bcftools",
                "/home/bernt-popp/miniforge3/envs/varlens-tools/bin/bcftools",
            ):
                if Path(candidate).exists():
                    real_bcftools = candidate
                    break
        if not real_bcftools:
            pytest.skip("bcftools binary not available")

        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(">chr1\nACGTACGTACGT\n")

        raw_vcf_content = (
            "##fileformat=VCFv4.2\n"
            '##FILTER=<ID=PASS,Description="All filters passed">\n'
            '##FILTER=<ID=FAIL,Description="Fail">\n'
            '##FILTER=<ID=weak_evidence,Description="Weak evidence">\n'
            '##FILTER=<ID=strand_bias,Description="Strand bias">\n'
            '##FILTER=<ID=clustered_events,Description="Clustered events">\n'
            '##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description="Filter status">\n'
            '##INFO=<ID=TLOD,Number=A,Type=Float,Description="Tumor LOD">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
            "##contig=<ID=chr1,length=12>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n"
            "chr1\t2\t.\tC\tG,T\t.\tPASS\tAS_FilterStatus=SITE|weak_evidence;TLOD=12.3,1.1\tGT:AD\t0/1/2:10,8,1\n"
            "chr1\t5\t.\tA\tC,G\t.\tclustered_events\tAS_FilterStatus=SITE|weak_evidence,strand_bias;TLOD=9.0,0.8\tGT:AD\t0/1/2:10,7,1\n"
            "chr1\t8\t.\tT\tA,C\t.\tFAIL\tAS_FilterStatus=weak_evidence|strand_bias;TLOD=1.2,1.5\tGT:AD\t0/1/2:10,1,1\n"
        )
        seed_vcf = tmp_path / "seed.vcf.gz"
        subprocess.run(
            [real_bcftools, "view", "-Oz", "-o", str(seed_vcf)],
            input=raw_vcf_content,
            text=True,
            check=True,
        )

        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        gatk_args_file = tmp_path / "gatk_args.txt"
        bcftools_args_file = tmp_path / "bcftools_args.txt"

        gatk_bin = bin_dir / "gatk"
        gatk_bin.write_text(
            "#!/usr/bin/env bash\n"
            f'printf "%s\\n" "$@" > "{gatk_args_file}"\n'
            "while [ $# -gt 0 ]; do\n"
            '  if [ "$1" = "-O" ]; then cp "' + str(seed_vcf) + '" "$2"; break; fi\n'
            "  shift\n"
            "done\n"
        )
        gatk_bin.chmod(0o755)

        bcftools_bin = bin_dir / "bcftools"
        bcftools_bin.write_text(
            "#!/usr/bin/env bash\n"
            f'printf "%s\\n" "$@" >> "{bcftools_args_file}"\n'
            f'exec "{real_bcftools}" "$@"\n'
        )
        bcftools_bin.chmod(0o755)

        out_vcf = str(tmp_path / "S1.filtered.vcf.gz")
        shell_tpl = _extract_rule_shell("workflow/rules/mutect2.smk", "filter_mutect_calls")
        script = shell_tpl.format(
            input=_DefaultNamespace(
                reference=str(ref_fa),
                vcf="in.vcf.gz",
                stats="in.vcf.gz.stats",
                contamination="S1.contamination.table",
                segmentation="S1.segments.table",
                ob_priors="S1_read-orientation-model.tar.gz",
            ),
            output=_DefaultNamespace(vcf=out_vcf, tbi=f"{out_vcf}.tbi"),
            params=_DefaultNamespace(
                java_opts="-Xmx4g",
                extra=ns["FILTER_MUTECT_CALLS_EXTRA"],
                contam_args=(
                    "--contamination-table S1.contamination.table "
                    "--tumor-segmentation S1.segments.table"
                ),
            ),
            log=str(tmp_path / "filter.log"),
        )
        env = {**os.environ, "PATH": f"{bin_dir}:{os.environ['PATH']}"}
        res = subprocess.run(
            ["bash", "-euo", "pipefail", "-c", script],
            env=env,
            capture_output=True,
            text=True,
        )
        assert res.returncode == 0, (tmp_path / "filter.log").read_text()
        gatk_args = gatk_args_file.read_text().splitlines()
        assert "--min-median-mapping-quality" in gatk_args
        assert "25" in gatk_args

        bcftools_args = bcftools_args_file.read_text().splitlines()
        assert "norm" in bcftools_args
        assert "-m-any" in bcftools_args
        assert "-W=tbi" in bcftools_args
        assert out_vcf in bcftools_args
        assert "-a" not in bcftools_args
        assert "--atomize" not in bcftools_args
        assert Path(f"{out_vcf}.tbi").exists()

        records = [
            line.split("\t")
            for line in subprocess.run(
                [real_bcftools, "view", "-H", out_vcf],
                capture_output=True,
                text=True,
                check=True,
            ).stdout.splitlines()
        ]
        assert len(records) == 6
        # chr1:2 C>G (SITE -> PASS) and C>T (weak_evidence -> weak_evidence)
        assert (records[0][1], records[0][4], records[0][6]) == ("2", "G", "PASS")
        assert "AS_FilterStatus=SITE;" in records[0][7]
        assert (records[1][1], records[1][4], records[1][6]) == ("2", "T", "weak_evidence")
        assert "AS_FilterStatus=weak_evidence;" in records[1][7]
        # chr1:5 A>C (clustered_events + SITE) and A>G (clustered_events + weak_evidence,strand_bias)
        assert (records[2][1], records[2][4], records[2][6]) == ("5", "C", "clustered_events")
        assert "AS_FilterStatus=SITE;" in records[2][7]
        assert (records[3][1], records[3][4], records[3][6]) == (
            "5",
            "G",
            "clustered_events;weak_evidence;strand_bias",
        )
        assert "AS_FilterStatus=weak_evidence,strand_bias;" in records[3][7]
        # chr1:8 T>A (FAIL + weak_evidence -> weak_evidence) and T>C (FAIL + strand_bias -> strand_bias)
        assert (records[4][1], records[4][4], records[4][6]) == ("8", "A", "weak_evidence")
        assert "AS_FilterStatus=weak_evidence;" in records[4][7]
        assert (records[5][1], records[5][4], records[5][6]) == ("8", "C", "strand_bias")
        assert "AS_FilterStatus=strand_bias;" in records[5][7]

    def test_multiallelic_as_filter_status_split_preserves_allele_and_site_filters(self):
        import re
        import shutil
        import subprocess
        from pathlib import Path

        real_bcftools = shutil.which("bcftools")
        if not real_bcftools:
            for candidate in (
                "/home/bernt-popp/miniforge3/envs/vntyper/bin/bcftools",
                "/home/bernt-popp/miniforge3/envs/varlens-tools/bin/bcftools",
            ):
                if Path(candidate).exists():
                    real_bcftools = candidate
                    break
        if not real_bcftools:
            pytest.skip("bcftools binary not available")

        shell_tpl = _extract_rule_shell("workflow/rules/mutect2.smk", "filter_mutect_calls")
        awk_blocks = re.findall(r"awk -F'\\t' -v OFS='\\t' '(.*?)'", shell_tpl, re.S)
        assert len(awk_blocks) == 2
        pre_awk = awk_blocks[0].replace("{{", "{").replace("}}", "}")
        post_awk = awk_blocks[1].replace("{{", "{").replace("}}", "}")

        raw_vcf_content = (
            "##fileformat=VCFv4.2\n"
            '##FILTER=<ID=PASS,Description="All filters passed">\n'
            '##FILTER=<ID=FAIL,Description="Fail">\n'
            '##FILTER=<ID=weak_evidence,Description="Weak evidence">\n'
            '##FILTER=<ID=strand_bias,Description="Strand bias">\n'
            '##FILTER=<ID=clustered_events,Description="Clustered events">\n'
            '##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description="Filter status">\n'
            '##INFO=<ID=TLOD,Number=A,Type=Float,Description="Tumor LOD">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
            "##contig=<ID=chr1,length=12>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n"
            "chr1\t2\t.\tA\tC,G\t.\tPASS\tAS_FilterStatus=SITE|weak_evidence;TLOD=12.3,1.1\tGT:AD\t0/1/2:10,8,1\n"
            "chr1\t5\t.\tA\tC,G\t.\tclustered_events\tAS_FilterStatus=SITE|weak_evidence,strand_bias;TLOD=9.0,0.8\tGT:AD\t0/1/2:10,7,1\n"
            "chr1\t8\t.\tT\tA,C\t.\tFAIL\tAS_FilterStatus=weak_evidence|strand_bias;TLOD=1.2,1.5\tGT:AD\t0/1/2:10,1,1\n"
        )
        cmd = (
            f"awk -F'\\t' -v OFS='\\t' '{pre_awk}' "
            f'| "{real_bcftools}" norm -m-any -Ov - '
            f"| awk -F'\\t' -v OFS='\\t' '{post_awk}'"
        )
        res = subprocess.run(
            ["bash", "-euo", "pipefail", "-c", cmd],
            input=raw_vcf_content,
            capture_output=True,
            text=True,
            check=True,
        )
        records = [
            line.split("\t") for line in res.stdout.splitlines() if not line.startswith("#")
        ]
        assert len(records) == 6
        assert (records[0][4], records[0][6], records[0][7]) == (
            "C",
            "PASS",
            "AS_FilterStatus=SITE;TLOD=12.3",
        )
        assert (records[1][4], records[1][6], records[1][7]) == (
            "G",
            "weak_evidence",
            "AS_FilterStatus=weak_evidence;TLOD=1.1",
        )
        assert (records[2][4], records[2][6], records[2][7]) == (
            "C",
            "clustered_events",
            "AS_FilterStatus=SITE;TLOD=9",
        )
        assert (records[3][4], records[3][6], records[3][7]) == (
            "G",
            "clustered_events;weak_evidence;strand_bias",
            "AS_FilterStatus=weak_evidence,strand_bias;TLOD=0.8",
        )
        assert (records[4][4], records[4][6], records[4][7]) == (
            "A",
            "weak_evidence",
            "AS_FilterStatus=weak_evidence;TLOD=1.2",
        )
        assert (records[5][4], records[5][6], records[5][7]) == (
            "C",
            "strand_bias",
            "AS_FilterStatus=strand_bias;TLOD=1.5",
        )


class TestCohortPonAndHotspotMasking:
    def test_in_pon_filtering_and_cohort_pon_dag_with_hotspot_masking(self, tmp_path, config_dict):
        import subprocess
        from pathlib import Path

        import pandas as pd

        df = pd.DataFrame(
            {
                "sample": ["S1", "S2", "S3", "G1", "G2"],
                "tumor_bam": ["S1.T", "S2.T", "S3.T", "G1.N", "G2.N"],
                "normal_bam": ["S1.N", "S2.N", ".", ".", "."],
                "analysis_type": [
                    "tumor_normal",
                    "tumor_normal",
                    "tumor_only",
                    "germline",
                    "germline",
                ],
                "in_pon": ["yes", "no", "no", "yes", "no"],
            }
        ).set_index("sample", drop=False)

        assert get_normal_bam_basenames(df) == ["G1.N", "S1.N"]

        root = Path(__file__).parent.parent
        cfg = dict(config_dict)
        cfg["gatk_resources"] = {
            **config_dict["gatk_resources"],
            "panel_of_normals": "",
            "pon_mask_bed": "hotspots_pad15.bed",
        }
        ns = {"config": cfg, "samples_df": df}
        exec((root / "workflow/rules/common.smk").read_text(), ns)

        assert ns["get_pon"]() == "/output/mutect2/pon/pon.masked.vcf.gz"
        assert ns["get_normal_bams"]() == ["G1.N", "S1.N"]

        cfg_unmasked = dict(config_dict)
        cfg_unmasked["gatk_resources"] = {
            **config_dict["gatk_resources"],
            "panel_of_normals": "",
            "pon_mask_bed": "",
        }
        ns_unmasked = {"config": cfg_unmasked, "samples_df": df}
        exec((root / "workflow/rules/common.smk").read_text(), ns_unmasked)
        assert ns_unmasked["get_pon"]() == "/output/mutect2/pon/pon.vcf.gz"

        pon_m2_shell = _extract_rule_shell("workflow/rules/mutect2.smk", "pon_mutect2")
        assert "--max-mnp-distance 0" in pon_m2_shell

        pon_db_shell = _extract_rule_shell("workflow/rules/mutect2.smk", "pon_genomicsdb")
        assert "GenomicsDBImport" in pon_db_shell
        assert "--merge-input-intervals true" in pon_db_shell

        pon_create_shell = _extract_rule_shell("workflow/rules/mutect2.smk", "pon_create")
        assert "CreateSomaticPanelOfNormals" in pon_create_shell
        assert "gendb://{input.db}" in pon_create_shell

        mask_shell = _extract_rule_shell("workflow/rules/mutect2.smk", "pon_mask_hotspots")
        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        bcftools_args_file = tmp_path / "bcftools_mask_args.txt"
        bcftools_bin = bin_dir / "bcftools"
        bcftools_bin.write_text(
            f'#!/usr/bin/env bash\nprintf "%s\\n" "$@" > "{bcftools_args_file}"\n'
        )
        bcftools_bin.chmod(0o755)

        script = mask_shell.format(
            input=_DefaultNamespace(
                vcf="/output/mutect2/pon/pon.vcf.gz",
                tbi="/output/mutect2/pon/pon.vcf.gz.tbi",
                mask_bed="hotspots_pad15.bed",
            ),
            output=_DefaultNamespace(
                vcf="/output/mutect2/pon/pon.masked.vcf.gz",
                tbi="/output/mutect2/pon/pon.masked.vcf.gz.tbi",
            ),
            log=str(tmp_path / "mask.log"),
        )
        env = {**os.environ, "PATH": f"{bin_dir}:{os.environ['PATH']}"}
        res = subprocess.run(
            ["bash", "-euo", "pipefail", "-c", script],
            env=env,
            capture_output=True,
            text=True,
        )
        assert res.returncode == 0
        args = bcftools_args_file.read_text().splitlines()
        assert args[:5] == ["view", "-T", "^hotspots_pad15.bed", "-Oz", "-W=tbi"]

    def test_pon_intervals_without_chromosome_scattering_or_targets(self, config_dict, samples_df):
        from pathlib import Path

        assert (
            build_pon_interval_arg("none", [], "/output/intervals", "")
            == "-L /output/intervals/scattered.interval_list"
        )
        assert (
            build_pon_interval_arg("interval", ["chr1"], "/output/intervals", "")
            == "-L /output/intervals/scattered.interval_list"
        )
        assert (
            build_pon_interval_arg("chromosome", ["chr1", "chr2"], "/output/intervals", "")
            == "-L chr1 -L chr2"
        )
        assert (
            build_pon_interval_arg("none", [], "/output/intervals", "targets.bed")
            == "-L targets.bed"
        )

        root = Path(__file__).parent.parent
        cfg = dict(config_dict)
        cfg["gatk_resources"] = {
            **config_dict["gatk_resources"],
            "panel_of_normals": "",
            "pon_mask_bed": "",
        }
        cfg["scatter"] = {
            "mode": "none",
            "chromosomes": [],
            "target_bed": "",
        }
        ns = {"config": cfg, "samples_df": samples_df}
        exec((root / "workflow/rules/common.smk").read_text(), ns)
        assert ns["get_pon_interval_input"]() == ["/output/intervals/scattered.interval_list"]
        assert ns["get_pon_interval_arg"]() == "-L /output/intervals/scattered.interval_list"


class TestMultiRegionTumorBams:
    def test_semicolon_delimited_tumor_bams_and_filter_mutect_calls(self, tmp_path, config_dict):
        import subprocess
        from pathlib import Path

        import pandas as pd

        df = pd.DataFrame(
            {
                "sample": ["MR1", "S1"],
                "tumor_bam": ["T1;T2", "S1.T"],
                "normal_bam": ["N1", "N1"],
                "analysis_type": ["tumor_normal", "tumor_normal"],
            }
        ).set_index("sample", drop=False)

        validate_samples_df(df)
        assert get_tumor_bam_basenames(df, "MR1") == ["T1", "T2"]
        assert get_unique_bam_basenames(df) == ["N1", "S1.T", "T1", "T2"]

        cfg = dict(config_dict)
        cfg["purecn"] = {"enabled": True, "genome": "hg38"}
        ns = {"config": cfg, "samples_df": df}
        root = Path(__file__).parent.parent
        exec((root / "workflow/rules/common.smk").read_text(), ns)

        wc = _DefaultNamespace(sample="MR1")
        assert ns["get_tumor_bam"](wc) == [
            "/bams/T1.merged.dedup.bqsr.bam",
            "/bams/T2.merged.dedup.bqsr.bam",
        ]
        assert (
            ns["get_tumor_input_arg"](wc)
            == "-I /bams/T1.merged.dedup.bqsr.bam -I /bams/T2.merged.dedup.bqsr.bam"
        )
        assert ns["get_contamination_tables"](wc) == [
            "/output/mutect2/contamination/MR1/T1.contamination.table",
            "/output/mutect2/contamination/MR1/T2.contamination.table",
        ]
        assert ns["get_segmentation_tables"](wc) == [
            "/output/mutect2/contamination/MR1/T1.segments.table",
            "/output/mutect2/contamination/MR1/T2.segments.table",
        ]

        final_outs = ns["get_final_outputs"]()
        assert "/output/purecn/results/S1/S1.csv" in final_outs
        assert "/output/purecn/results/MR1/MR1.csv" not in final_outs

        res, gatk_args = TestMutect2CallSmResolution._run_mutect2_shell(
            tmp_path,
            tumor_headers={
                "T1.bam": "@RG\tID:1\tSM:T1\n",
                "T2.bam": "@RG\tID:2\tSM:T2\n",
            },
            normal_header="@RG\tID:3\tSM:N1\n",
        )
        assert res.returncode == 0
        assert gatk_args.count("-I") == 3

        bin_dir = tmp_path / "fbin"
        bin_dir.mkdir()
        gatk_args_file = tmp_path / "filter_gatk_args.txt"
        for tool in ("gatk", "bcftools"):
            p = bin_dir / tool
            p.write_text(
                f'#!/usr/bin/env bash\n[ "{tool}" = "gatk" ] && printf "%s\\n" "$@" > "{gatk_args_file}" || true\n'
            )
            p.chmod(0o755)

        out_vcf = str(tmp_path / "MR1.filtered.vcf.gz")
        shell_tpl = _extract_rule_shell("workflow/rules/mutect2.smk", "filter_mutect_calls")
        script = shell_tpl.format(
            input=_DefaultNamespace(
                reference="ref.fa",
                vcf="in.vcf.gz",
                stats="in.vcf.gz.stats",
                contamination=" ".join(ns["get_contamination_tables"](wc)),
                segmentation=" ".join(ns["get_segmentation_tables"](wc)),
                ob_priors="MR1_read-orientation-model.tar.gz",
            ),
            output=_DefaultNamespace(vcf=out_vcf, tbi=f"{out_vcf}.tbi"),
            params=_DefaultNamespace(
                java_opts="-Xmx4g",
                contam_args=ns["get_contamination_args"](wc),
                extra="",
            ),
            log=str(tmp_path / "filter_mr.log"),
        )
        env = {**os.environ, "PATH": f"{bin_dir}:{os.environ['PATH']}"}
        res_f = subprocess.run(
            ["bash", "-euo", "pipefail", "-c", script],
            env=env,
            capture_output=True,
            text=True,
        )
        assert res_f.returncode == 0
        f_args = gatk_args_file.read_text().splitlines()
        assert f_args.count("--contamination-table") == 2
        assert f_args.count("--tumor-segmentation") == 2


@pytest.mark.dryrun
class TestSnakemakeDryRun:
    def test_snakemake_dryrun_interval_pon_and_underscored_multiregion(self, tmp_path):
        import shutil
        import subprocess
        from pathlib import Path

        import yaml

        snakemake_bin = shutil.which("snakemake") or "/tmp/sm-venv/bin/snakemake"
        if not Path(snakemake_bin).exists():
            pytest.skip("snakemake binary not available")

        root = Path(__file__).parent.parent
        bam_dir = tmp_path / "bams"
        bam_dir.mkdir()
        for bn in (
            "APA_1_T1",
            "APA_1_T2",
            "APA_1_N",
            "APA_2_T",
            "APA_2_N",
            "GERM_1_N",
            "GERM_2_N",
        ):
            (bam_dir / f"{bn}.merged.dedup.bqsr.bam").write_text("")

        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(">chr1\nACGT\n")
        gnomad_vcf = tmp_path / "gnomad.vcf.gz"
        gnomad_vcf.write_text("")
        common_gnomad_vcf = tmp_path / "common_gnomad.vcf.gz"
        common_gnomad_vcf.write_text("")
        target_bed = tmp_path / "targets.bed"
        target_bed.write_text("chr1\t0\t4\n")
        mask_bed = tmp_path / "hotspots_pad15.bed"
        mask_bed.write_text("chr1\t1\t2\n")

        samples_tsv = tmp_path / "samples.tsv"
        samples_tsv.write_text(
            "sample\ttumor_bam\tnormal_bam\tanalysis_type\tin_pon\n"
            "APA_1\tAPA_1_T1;APA_1_T2\tAPA_1_N\ttumor_normal\tyes\n"
            "APA_2\tAPA_2_T\tAPA_2_N\ttumor_normal\tno\n"
            "GERM_1\tGERM_1_N\t.\tgermline\tyes\n"
            "GERM_2\tGERM_2_N\t.\tgermline\tno\n"
        )

        cfg = {
            "caller": "all",
            "paths": {
                "samples": str(samples_tsv),
                "bam_folder": str(bam_dir),
                "output_folder": str(tmp_path / "out"),
                "log_subdir": "logs",
            },
            "ref": {
                "genome": str(ref_fa),
                "build": "GRCh38",
            },
            "bam": {
                "file_extension": ".merged.dedup.bqsr.bam",
            },
            "gatk_resources": {
                "af_only_gnomad": str(gnomad_vcf),
                "panel_of_normals": "",
                "pon_mask_bed": str(mask_bed),
                "common_biallelic_gnomad": str(common_gnomad_vcf),
            },
            "scatter": {
                "mode": "interval",
                "count": 2,
                "chromosomes": ["chr1"],
                "target_bed": str(target_bed),
            },
        }
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(yaml.safe_dump(cfg))

        res = subprocess.run(
            [
                snakemake_bin,
                "-s",
                str(root / "workflow/Snakefile"),
                "--configfile",
                str(cfg_file),
                "--workflow-profile",
                str(root / "profiles/default"),
                "--profile",
                str(root / "profiles/bih"),
                "--executor",
                "dryrun",
                "-n",
                "-p",
            ],
            cwd=str(tmp_path),
            capture_output=True,
            text=True,
        )
        assert res.returncode == 0, res.stderr + "\n" + res.stdout
        out = res.stdout

        assert "rule split_intervals:" in out
        assert "rule pon_mutect2:" in out
        assert "APA_1_N.pon.vcf.gz" in out
        assert "GERM_1_N.pon.vcf.gz" in out
        assert "APA_2_N.pon.vcf.gz" not in out
        assert "GERM_2_N.pon.vcf.gz" not in out
        assert "rule pon_genomicsdb:" in out
        assert "rule pon_create:" in out
        assert "rule pon_mask_hotspots:" in out
        assert "APA_1/APA_1_T1.contamination.table" in out
        assert "APA_1/APA_1_T2.contamination.table" in out
        assert "APA_2/APA_2_T.contamination.table" in out
        assert "--matched " in out
        assert "APA_1_N.getpileupsummaries.table" in out
        assert "APA_2_N.getpileupsummaries.table" in out

        # Verify PoN in mode=none with chromosomes=[], target_bed="", panel_of_normals=""
        cfg_none = {
            **cfg,
            "scatter": {
                "mode": "none",
                "chromosomes": [],
                "target_bed": "",
            },
        }
        cfg_none_file = tmp_path / "config_none.yaml"
        cfg_none_file.write_text(yaml.safe_dump(cfg_none))
        res_none = subprocess.run(
            [
                snakemake_bin,
                "-s",
                str(root / "workflow/Snakefile"),
                "--configfile",
                str(cfg_none_file),
                "--executor",
                "dryrun",
                "-n",
                "-p",
            ],
            cwd=str(tmp_path),
            capture_output=True,
            text=True,
        )
        assert res_none.returncode == 0, res_none.stderr + "\n" + res_none.stdout
        assert "rule scatter_intervals_by_ns:" in res_none.stdout
        assert "scattered.interval_list" in res_none.stdout

        # Verify chromosome mode with target_bed omits configured chromosomes absent from BED
        cfg_chr = {
            **cfg,
            "scatter": {
                "mode": "chromosome",
                "chromosomes": ["chr1", "chrM"],
                "target_bed": str(target_bed),
            },
        }
        cfg_chr_file = tmp_path / "config_chr.yaml"
        cfg_chr_file.write_text(yaml.safe_dump(cfg_chr))
        res_chr = subprocess.run(
            [
                snakemake_bin,
                "-s",
                str(root / "workflow/Snakefile"),
                "--configfile",
                str(cfg_chr_file),
                "--executor",
                "dryrun",
                "-n",
                "-p",
            ],
            cwd=str(tmp_path),
            capture_output=True,
            text=True,
        )
        assert res_chr.returncode == 0, res_chr.stderr + "\n" + res_chr.stdout
        assert "APA_1_chr1.vcf.gz" in res_chr.stdout
        assert "freebayes_chr1.vcf.gz" in res_chr.stdout
        assert "chrM" not in res_chr.stdout
