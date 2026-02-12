"""Tests for config and samples schema validation."""

import copy

import pytest
import yaml

try:
    import jsonschema
except ImportError:
    jsonschema = None

SCHEMA_DIR = "workflow/schemas"


def _load_schema(name):
    with open(f"{SCHEMA_DIR}/{name}") as f:
        return yaml.safe_load(f)


@pytest.fixture
def config_schema():
    return _load_schema("config.schema.yaml")


@pytest.fixture
def samples_schema():
    return _load_schema("samples.schema.yaml")


@pytest.mark.skipif(jsonschema is None, reason="jsonschema not installed")
class TestConfigSchema:
    def test_valid_config(self, config_dict, config_schema):
        jsonschema.validate(config_dict, config_schema)

    def test_missing_caller(self, config_dict, config_schema):
        broken = copy.deepcopy(config_dict)
        del broken["caller"]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, config_schema)

    def test_invalid_caller(self, config_dict, config_schema):
        broken = copy.deepcopy(config_dict)
        broken["caller"] = "invalid"
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, config_schema)

    def test_invalid_scatter_mode(self, config_dict, config_schema):
        broken = copy.deepcopy(config_dict)
        broken["scatter"]["mode"] = "invalid"
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, config_schema)

    def test_missing_ref(self, config_dict, config_schema):
        broken = copy.deepcopy(config_dict)
        del broken["ref"]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, config_schema)


    def test_mutect2_genotype_germline_sites_boolean(self, config_dict, config_schema):
        cfg = copy.deepcopy(config_dict)
        cfg["params"]["mutect2"]["genotype_germline_sites"] = True
        jsonschema.validate(cfg, config_schema)

    def test_mutect2_genotype_germline_sites_rejects_string(
        self, config_dict, config_schema
    ):
        cfg = copy.deepcopy(config_dict)
        cfg["params"]["mutect2"]["genotype_germline_sites"] = "yes"
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(cfg, config_schema)

    def test_mutect2_annotations_array(self, config_dict, config_schema):
        cfg = copy.deepcopy(config_dict)
        cfg["params"]["mutect2"]["annotations"] = [
            "OrientationBiasReadCounts",
            "StrandBiasBySample",
        ]
        jsonschema.validate(cfg, config_schema)

    def test_purecn_enabled_boolean(self, config_dict, config_schema):
        cfg = copy.deepcopy(config_dict)
        cfg["purecn"] = {"enabled": True, "genome": "hg38"}
        jsonschema.validate(cfg, config_schema)

    def test_purecn_genome_validates_enum(self, config_dict, config_schema):
        cfg = copy.deepcopy(config_dict)
        cfg["purecn"] = {"enabled": True, "genome": "invalid"}
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(cfg, config_schema)

    def test_purecn_seed_integer(self, config_dict, config_schema):
        cfg = copy.deepcopy(config_dict)
        cfg["purecn"] = {"enabled": False, "genome": "hg38", "seed": 42}
        jsonschema.validate(cfg, config_schema)


@pytest.mark.skipif(jsonschema is None, reason="jsonschema not installed")
class TestSamplesSchema:
    def test_valid_sample(self, samples_schema):
        sample = {
            "sample": "IND001_To",
            "tumor_bam": "IND001.tumor",
            "normal_bam": ".",
            "analysis_type": "tumor_only",
        }
        jsonschema.validate(sample, samples_schema)

    def test_missing_sample(self, samples_schema):
        broken = {
            "tumor_bam": "IND001.tumor",
            "analysis_type": "tumor_only",
        }
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, samples_schema)

    def test_invalid_analysis_type(self, samples_schema):
        broken = {
            "sample": "IND001",
            "tumor_bam": "IND001.tumor",
            "analysis_type": "invalid",
        }
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, samples_schema)

    def test_missing_tumor_bam(self, samples_schema):
        broken = {
            "sample": "IND001",
            "analysis_type": "tumor_only",
        }
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, samples_schema)
