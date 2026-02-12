#!/usr/bin/env python3
"""Generate config/samples.tsv and config/config.yaml for the sm-calling pipeline.

Scans a BAM directory, classifies tumor/normal roles via filename heuristics,
pairs tumors with normals, and generates the pipeline config files.

Two modes:
  Interactive:  python scripts/generate_config.py          (guided wizard)
  Flags:        python scripts/generate_config.py --bam-folder /path/to/bams

Usage:
    python scripts/generate_config.py
    python scripts/generate_config.py --bam-folder /data/bams
    python scripts/generate_config.py --bam-folder /data/bams --dry-run
    python scripts/generate_config.py --bam-folder /data/bams --config-template --caller mutect2
"""

from __future__ import annotations

import argparse
import re
import sys
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any

try:
    import pandas as pd
except ImportError:
    sys.exit(
        "Error: pandas is required but not installed.\n"
        "Install it with: pip install pandas\n"
        "(pandas is already available in Snakemake conda environments.)"
    )


# ---------------------------------------------------------------------------
# Section A: Constants & Patterns
# ---------------------------------------------------------------------------

DEFAULT_BAM_EXTENSION = ".merged.dedup.bqsr.bam"

# Patterns that identify tumor BAMs (applied to the basename without extension)
TUMOR_PATTERNS = [
    (
        re.compile(r"[._-](?:tumor|tumour|FFPE|met|metastatic|primary|tum)$", re.I),
        "suffix '{match}'",
    ),
    (
        re.compile(r"[._-]T[A-Z]?\d*$"),
        "suffix '{match}'",
    ),
]

# Patterns that identify normal BAMs
NORMAL_PATTERNS = [
    (
        re.compile(r"[._-](?:normal|blood|buffy|germline|pbmc|ctrl)$", re.I),
        "suffix '{match}'",
    ),
    (
        re.compile(r"[._-]N\d*$"),
        "suffix '{match}'",
    ),
]

# Combined pattern for stripping role suffixes to extract patient ID
ROLE_SUFFIX_RE = re.compile(
    r"[._-](?:tumor|tumour|FFPE|met|metastatic|primary|tum|"
    r"normal|blood|buffy|germline|pbmc|ctrl|"
    r"T[A-Z]?\d*|N\d*)$",
    re.I,
)

# Reference genome file extensions
GENOME_EXTENSIONS = ("*.fna", "*.fa", "*.fasta")

# GATK resource search patterns
GATK_RESOURCE_PATTERNS = {
    "panel_of_normals": ["*pon*.vcf*", "*1000g_pon*.vcf*"],
    "af_only_gnomad": ["*af-only-gnomad*.vcf*"],
    "common_biallelic_gnomad": ["*common_biallelic*.vcf*"],
}

# Standard search directories for reference data (relative to project root)
REF_SEARCH_DIRS = [
    "resources/ref",
    "resources/ref/GRCh38",
    "analysis/ref/GRCh38",
    "../resources/ref/GRCh38",
    "../resources/ref",
]

# GATK resource bundle search directories
GATK_RESOURCE_SEARCH_DIRS = [
    "resources/gatk_bundle/hg38",
    "resources/gatk_bundle",
    "analysis/GATK_resource_bundle",
    "../resources/gatk_bundle/hg38",
    "../resources/gatk_bundle",
]

# Well-known shared locations on BIH and Charite HPC
SHARED_REF_DIRS = [
    "/data/cephfs-1/work/groups/scholl/shared/ref/GRCh38",
    "/data/cephfs-1/work/groups/scholl/shared/ref",
]

SHARED_GATK_DIRS = [
    "/data/cephfs-1/work/projects/apa-sequencing/analysis/GATK_resource_bundle",
]

# Default chromosomes for scatter mode
DEFAULT_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


# ---------------------------------------------------------------------------
# Section B: BAM Discovery
# ---------------------------------------------------------------------------


def discover_bam_files(bam_folder: str | Path, bam_extension: str) -> list[dict]:
    """Scan a directory for BAM files matching the given extension.

    Args:
        bam_folder: Directory to scan.
        bam_extension: File extension to match (e.g., '.merged.dedup.bqsr.bam').

    Returns:
        List of dicts with keys: basename, path, size_bytes, has_index.
        Sorted alphabetically by basename.
    """
    bam_dir = Path(bam_folder)
    if not bam_dir.is_dir():
        return []

    results = []
    all_files = sorted(bam_dir.iterdir())
    bam_files = [f for f in all_files if f.is_file() and f.name.endswith(bam_extension)]

    for bam_path in bam_files:
        basename = bam_path.name[: -len(bam_extension)]
        # Check for .bai index (both conventions: file.bam.bai and file.bai)
        has_index = (
            (bam_dir / (bam_path.name + ".bai")).is_file()
            or (bam_dir / (basename + bam_extension[:-4] + ".bai")).is_file()
            or (bam_dir / (bam_path.name[:-4] + ".bai")).is_file()
        )
        results.append(
            {
                "basename": basename,
                "path": str(bam_path),
                "size_bytes": bam_path.stat().st_size,
                "has_index": has_index,
            }
        )

    return results


def format_size(size_bytes: int) -> str:
    """Format file size in human-readable form.

    Args:
        size_bytes: File size in bytes.

    Returns:
        Formatted string like '45.3 GB', '512 MB', '1.2 KB'.
    """
    if size_bytes >= 1_000_000_000:
        return f"{size_bytes / 1_000_000_000:.1f} GB"
    if size_bytes >= 1_000_000:
        return f"{size_bytes / 1_000_000:.1f} MB"
    if size_bytes >= 1_000:
        return f"{size_bytes / 1_000:.1f} KB"
    return f"{size_bytes} B"


# ---------------------------------------------------------------------------
# Section C: Role Detection + Patient Grouping
# ---------------------------------------------------------------------------


def guess_role(basename: str) -> tuple[str, str]:
    """Guess whether a BAM basename is tumor, normal, or unknown.

    Args:
        basename: BAM file basename (without extension).

    Returns:
        Tuple of (role, reason) where role is 'tumor', 'normal', or 'unknown'.
    """
    for pattern, reason_tpl in TUMOR_PATTERNS:
        m = pattern.search(basename)
        if m:
            return ("tumor", reason_tpl.format(match=m.group()))

    for pattern, reason_tpl in NORMAL_PATTERNS:
        m = pattern.search(basename)
        if m:
            return ("normal", reason_tpl.format(match=m.group()))

    return ("unknown", "no pattern matched")


def extract_patient_id(basename: str) -> str:
    """Extract patient ID by stripping the role suffix from a BAM basename.

    Args:
        basename: BAM file basename (without extension).

    Returns:
        Patient ID string. If no role suffix found, returns the basename as-is.
    """
    return ROLE_SUFFIX_RE.sub("", basename)


def find_best_normal(
    tumor_basename: str,
    available_normals: list[str],
    patient_groups: dict[str, dict[str, list[str]]] | None = None,
) -> str | None:
    """Find the best matching normal BAM for a tumor BAM.

    Priority:
        1. Same patient_id (from patient_groups)
        2. Highest string similarity (SequenceMatcher)

    Args:
        tumor_basename: Basename of the tumor BAM.
        available_normals: List of normal BAM basenames.
        patient_groups: Optional patient grouping dict.

    Returns:
        Best matching normal basename, or None if no match above threshold.
    """
    if not available_normals:
        return None

    # Priority 1: same patient ID
    if patient_groups is not None:
        tumor_pid = extract_patient_id(tumor_basename)
        for normal in available_normals:
            normal_pid = extract_patient_id(normal)
            if tumor_pid == normal_pid:
                return normal

    # Priority 2: string similarity
    best, best_score = None, 0.0
    for normal in available_normals:
        score = SequenceMatcher(None, tumor_basename, normal).ratio()
        if score > best_score:
            best, best_score = normal, score

    return best if best_score > 0.5 else None


def group_and_pair(bam_entries: list[dict]) -> tuple[dict[str, dict[str, list[str]]], list[dict]]:
    """Group BAMs by patient ID and pair tumors with normals.

    Args:
        bam_entries: List of dicts with at least 'basename' and 'role' keys.

    Returns:
        Tuple of (patient_groups, pairings).
        patient_groups: {patient_id: {tumors: [...], normals: [...]}}.
        pairings: [{tumor: basename, normal: basename|None, patient_id: str}, ...].
    """
    # Group by patient ID
    patient_groups: dict[str, dict[str, list[str]]] = {}
    for entry in bam_entries:
        if entry.get("role") == "skip":
            continue
        pid = extract_patient_id(entry["basename"])
        if pid not in patient_groups:
            patient_groups[pid] = {"tumors": [], "normals": []}
        role = entry.get("role", "unknown")
        if role == "tumor":
            patient_groups[pid]["tumors"].append(entry["basename"])
        elif role == "normal":
            patient_groups[pid]["normals"].append(entry["basename"])

    # Pair tumors with normals
    pairings: list[dict] = []
    used_normals: set[str] = set()

    # All normals (for cross-patient fallback)
    all_normals = [
        e["basename"] for e in bam_entries if e.get("role") == "normal"
    ]

    for pid, group in patient_groups.items():
        for tumor in group["tumors"]:
            # Try same-patient normals first
            normal = None
            for n in group["normals"]:
                if n not in used_normals:
                    normal = n
                    break

            # If no same-patient normal, try cross-patient via similarity
            if normal is None:
                remaining = [n for n in all_normals if n not in used_normals]
                normal = find_best_normal(tumor, remaining, patient_groups)

            if normal is not None:
                used_normals.add(normal)

            pairings.append(
                {
                    "tumor": tumor,
                    "normal": normal,
                    "patient_id": pid,
                }
            )

    return patient_groups, pairings


def generate_sample_name(patient_id: str, analysis_type: str, index: int, total: int) -> str:
    """Generate a sample name from patient ID and analysis type.

    Args:
        patient_id: Patient identifier.
        analysis_type: One of 'tumor_only', 'tumor_normal', 'germline'.
        index: 0-based index for disambiguation.
        total: Total samples of this type for this patient.

    Returns:
        Sample name like 'APA1_TN', 'APA1_To', 'APA1_G', or 'APA1_TN_2'.
    """
    suffix_map = {
        "tumor_only": "To",
        "tumor_normal": "TN",
        "germline": "G",
    }
    suffix = suffix_map.get(analysis_type, analysis_type)
    name = f"{patient_id}_{suffix}"
    if total > 1:
        name = f"{name}_{index + 1}"
    return name


def generate_rows(
    pairings: list[dict],
    bam_entries: list[dict],
    analysis_mode: str,
) -> list[dict]:
    """Generate sample rows for the given analysis mode.

    Args:
        pairings: From group_and_pair().
        bam_entries: Original BAM entries (for germline mode).
        analysis_mode: One of 'both', 'tumor_only', 'tumor_normal', 'germline'.

    Returns:
        List of dicts with keys: sample, tumor_bam, normal_bam, analysis_type, patient_id.
    """
    rows: list[dict] = []

    if analysis_mode == "germline":
        # One row per BAM (use tumor_bam column for the BAM, regardless of role)
        patient_counts: dict[str, int] = {}
        for entry in bam_entries:
            if entry.get("role") == "skip":
                continue
            pid = extract_patient_id(entry["basename"])
            patient_counts[pid] = patient_counts.get(pid, 0) + 1

        patient_indices: dict[str, int] = {}
        for entry in bam_entries:
            if entry.get("role") == "skip":
                continue
            pid = extract_patient_id(entry["basename"])
            idx = patient_indices.get(pid, 0)
            patient_indices[pid] = idx + 1
            sample = generate_sample_name(pid, "germline", idx, patient_counts[pid])
            rows.append(
                {
                    "sample": sample,
                    "tumor_bam": entry["basename"],
                    "normal_bam": ".",
                    "analysis_type": "germline",
                    "patient_id": pid,
                }
            )
        return rows

    # Count per-patient for disambiguation
    patient_type_counts: dict[str, dict[str, int]] = {}
    for pairing in pairings:
        pid = pairing["patient_id"]
        if pid not in patient_type_counts:
            patient_type_counts[pid] = {"tumor_only": 0, "tumor_normal": 0}
        if analysis_mode in ("both", "tumor_only"):
            patient_type_counts[pid]["tumor_only"] += 1
        if analysis_mode in ("both", "tumor_normal") and pairing["normal"] is not None:
            patient_type_counts[pid]["tumor_normal"] += 1

    patient_type_indices: dict[str, dict[str, int]] = {}

    for pairing in pairings:
        pid = pairing["patient_id"]
        if pid not in patient_type_indices:
            patient_type_indices[pid] = {"tumor_only": 0, "tumor_normal": 0}

        # Tumor-only row
        if analysis_mode in ("both", "tumor_only"):
            idx = patient_type_indices[pid]["tumor_only"]
            total = patient_type_counts[pid]["tumor_only"]
            patient_type_indices[pid]["tumor_only"] = idx + 1
            sample = generate_sample_name(pid, "tumor_only", idx, total)
            rows.append(
                {
                    "sample": sample,
                    "tumor_bam": pairing["tumor"],
                    "normal_bam": ".",
                    "analysis_type": "tumor_only",
                    "patient_id": pid,
                }
            )

        # Tumor-normal row
        if analysis_mode in ("both", "tumor_normal") and pairing["normal"] is not None:
            idx = patient_type_indices[pid]["tumor_normal"]
            total = patient_type_counts[pid]["tumor_normal"]
            patient_type_indices[pid]["tumor_normal"] = idx + 1
            sample = generate_sample_name(pid, "tumor_normal", idx, total)
            rows.append(
                {
                    "sample": sample,
                    "tumor_bam": pairing["tumor"],
                    "normal_bam": pairing["normal"],
                    "analysis_type": "tumor_normal",
                    "patient_id": pid,
                }
            )

    return rows


# ---------------------------------------------------------------------------
# Section D: Reference + GATK Resource Discovery
# ---------------------------------------------------------------------------


def _find_genome_fastas(search_dir: Path) -> list[dict[str, Any]]:
    """Find reference genome FASTA files in a directory.

    Args:
        search_dir: Directory to search.

    Returns:
        List of dicts with keys: path, has_fai, name.
    """
    results: list[dict[str, Any]] = []
    if not search_dir.is_dir():
        return results
    for ext in GENOME_EXTENSIONS:
        for fasta in search_dir.glob(ext):
            if fasta.name.startswith("."):
                continue
            has_fai = (fasta.parent / (fasta.name + ".fai")).is_file()
            results.append(
                {
                    "path": str(fasta),
                    "has_fai": has_fai,
                    "name": fasta.name,
                }
            )
    return results


def _find_gatk_resource_vcfs(
    search_dir: Path,
    resource_key: str,
    patterns: list[str],
) -> list[str]:
    """Find GATK resource VCF files matching glob patterns.

    Args:
        search_dir: Directory to search.
        resource_key: Resource identifier (for exclusion logic).
        patterns: Glob patterns to match.

    Returns:
        List of file path strings (deduplicated).
    """
    results: list[str] = []
    if not search_dir.is_dir():
        return results
    for pattern in patterns:
        for vcf in search_dir.glob(pattern):
            if vcf.name.endswith(".tbi") or vcf.name.endswith(".idx"):
                continue
            # For af_only_gnomad, exclude files with 'common' in the name
            if resource_key == "af_only_gnomad" and "common" in vcf.name.lower():
                continue
            results.append(str(vcf))
    return sorted(set(results))


def discover_reference_data(
    ref_dir: Path | None,
    project_root: Path | None = None,
) -> dict[str, Any]:
    """Scan for reference genome FASTA files.

    Args:
        ref_dir: Explicit reference directory to search first.
        project_root: Project root for relative path searches.

    Returns:
        Dict with keys: genome, build, search_log.
    """
    search_log: list[str] = []
    genomes: list[dict[str, Any]] = []

    ref_search: list[Path] = []
    if ref_dir:
        ref_search.append(ref_dir)
    if project_root:
        for rel in REF_SEARCH_DIRS:
            ref_search.append(project_root / rel)
    for shared in SHARED_REF_DIRS:
        ref_search.append(Path(shared))

    for search_path in ref_search:
        if not search_path.is_dir():
            continue
        found = _find_genome_fastas(search_path)
        if found:
            search_log.append(f"  Found {len(found)} genome(s) in {search_path}")
            for g in found:
                status = "FAI" if g["has_fai"] else "no index"
                search_log.append(f"    {g['name']}  [{status}]")
            genomes.extend(found)

    genome_path = ""
    build = "GRCh38"
    if genomes:
        best = sorted(genomes, key=lambda g: (g["has_fai"], bool(g["path"])), reverse=True)[0]
        genome_path = str(Path(best["path"]).resolve()) if best["path"] else ""
        name_lower = best["name"].lower()
        if "grch37" in name_lower or "hg19" in name_lower or "hs37" in name_lower:
            build = "GRCh37"

    if not genomes:
        search_log.append("  No reference genome found")

    return {
        "genome": genome_path,
        "build": build,
        "search_log": search_log,
    }


def discover_gatk_resources(
    resource_dir: Path | None,
    project_root: Path | None = None,
) -> dict[str, Any]:
    """Scan for GATK resource VCFs (PoN, gnomAD).

    Args:
        resource_dir: Explicit resource directory to search first.
        project_root: Project root for relative path searches.

    Returns:
        Dict with keys: panel_of_normals, af_only_gnomad,
        common_biallelic_gnomad, search_log.
    """
    search_log: list[str] = []
    found_resources: dict[str, str] = {
        "panel_of_normals": "",
        "af_only_gnomad": "",
        "common_biallelic_gnomad": "",
    }

    search_dirs: list[Path] = []
    if resource_dir:
        search_dirs.append(resource_dir)
    if project_root:
        for rel in GATK_RESOURCE_SEARCH_DIRS:
            search_dirs.append(project_root / rel)
    for shared in SHARED_GATK_DIRS:
        search_dirs.append(Path(shared))

    for search_path in search_dirs:
        if not search_path.is_dir():
            continue
        for key, patterns in GATK_RESOURCE_PATTERNS.items():
            if found_resources[key]:
                continue  # already found
            matches = _find_gatk_resource_vcfs(search_path, key, patterns)
            if matches:
                found_resources[key] = matches[0]
                search_log.append(f"  Found {key}: {Path(matches[0]).name}")

    for key, path in found_resources.items():
        if not path:
            search_log.append(f"  {key}: not found")

    return {**found_resources, "search_log": search_log}


# ---------------------------------------------------------------------------
# Section E: Config Template
# ---------------------------------------------------------------------------


def _resolve_path(p: str) -> str:
    """Resolve a path string to canonical form, skipping placeholders."""
    if not p or "EDIT_ME" in p:
        return p
    return str(Path(p).resolve()).replace("\\", "/")


def _build_config_yaml(
    caller: str,
    ref_data: dict[str, Any],
    gatk_resources: dict[str, Any],
    bam_folder: str,
    output_folder: str,
    samples_path: str,
    bam_extension: str,
    scatter_mode: str,
) -> str:
    """Build config.yaml content matching config.schema.yaml.

    Uses string formatting (not yaml.dump) to preserve comments.

    Args:
        caller: Variant caller choice.
        ref_data: From discover_reference_data().
        gatk_resources: From discover_gatk_resources().
        bam_folder: BAM input directory.
        output_folder: Pipeline output directory.
        samples_path: Path to samples.tsv.
        bam_extension: BAM file extension.
        scatter_mode: Scatter mode (chromosome/interval/none).

    Returns:
        YAML config file content as a string.
    """
    genome = ref_data.get("genome", "") or "EDIT_ME: /path/to/reference.fna"
    build = ref_data.get("build", "GRCh38")

    pon = gatk_resources.get("panel_of_normals", "") or "EDIT_ME: /path/to/pon.vcf.gz"
    af_gnomad = (
        gatk_resources.get("af_only_gnomad", "") or "EDIT_ME: /path/to/af-only-gnomad.vcf.gz"
    )
    common_gnomad = (
        gatk_resources.get("common_biallelic_gnomad", "")
        or "EDIT_ME: /path/to/common_biallelic.vcf.gz"
    )

    genome = _resolve_path(genome)
    pon = _resolve_path(pon)
    af_gnomad = _resolve_path(af_gnomad)
    common_gnomad = _resolve_path(common_gnomad)
    bam_folder = _resolve_path(bam_folder)

    # Build chromosomes list
    chrom_lines = "\n".join(f'    - "{c}"' for c in DEFAULT_CHROMOSOMES)

    return f"""\
# sm-calling unified configuration
# Generated by: scripts/generate_config.py
# Review all paths below before running the pipeline.

# --- Caller selection ---
caller: "{caller}"  # "mutect2" | "freebayes" | "all"

# --- Reference genome ---
ref:
  genome: "{genome}"
  build: "{build}"

# --- GATK resources (required for Mutect2) ---
gatk_resources:
  panel_of_normals: "{pon}"
  af_only_gnomad: "{af_gnomad}"
  common_biallelic_gnomad: "{common_gnomad}"

# --- Paths ---
paths:
  samples: "{samples_path}"
  bam_folder: "{bam_folder}"
  output_folder: "{output_folder}"
  log_subdir: "logs"
  intervals_dir: "analysis/intervals"

# --- BAM file settings ---
bam:
  file_extension: "{bam_extension}"

# --- Scatter settings ---
scatter:
  mode: "{scatter_mode}"  # "chromosome" | "interval" | "none"
  count: 400
  chromosomes:
{chrom_lines}

# --- Tool-specific parameters (passthrough) ---
params:
  mutect2:
    extra: "--genotype-germline-sites true --genotype-pon-sites true"
  freebayes:
    extra: "--min-coverage 20 --limit-coverage 500 --use-best-n-alleles 4 --standard-filters"
  bcftools_norm:
    extra: "-m-any --force -a --atom-overlaps ."
"""


def generate_config_template(
    config_output: Path,
    caller: str,
    ref_data: dict[str, Any],
    gatk_resources: dict[str, Any],
    bam_folder: str,
    output_folder: str,
    samples_path: str,
    bam_extension: str,
    scatter_mode: str,
    dry_run: bool = False,
    force: bool = False,
) -> None:
    """Write config.yaml with overwrite protection.

    Args:
        config_output: Path for the output config file.
        caller: Variant caller choice.
        ref_data: From discover_reference_data().
        gatk_resources: From discover_gatk_resources().
        bam_folder: BAM input directory.
        output_folder: Pipeline output directory.
        samples_path: Path to samples.tsv.
        bam_extension: BAM file extension.
        scatter_mode: Scatter mode.
        dry_run: If True, print but don't write.
        force: If True, overwrite without asking.
    """
    content = _build_config_yaml(
        caller=caller,
        ref_data=ref_data,
        gatk_resources=gatk_resources,
        bam_folder=bam_folder,
        output_folder=output_folder,
        samples_path=samples_path,
        bam_extension=bam_extension,
        scatter_mode=scatter_mode,
    )

    if dry_run:
        print(f"\n--- Config template ({config_output}) ---")
        print(content)
        return

    if config_output.is_file() and not force:
        response = input(f"Config file already exists: {config_output}. Overwrite? [y/N] ")
        if response.lower() not in ("y", "yes"):
            print(f"Skipped writing {config_output}")
            return

    config_output.parent.mkdir(parents=True, exist_ok=True)
    with open(config_output, "w", encoding="utf-8") as fh:
        fh.write(content)
    print(f"Written: {config_output}")


# ---------------------------------------------------------------------------
# Section F: Output
# ---------------------------------------------------------------------------


def build_samples_dataframe(rows: list[dict]) -> pd.DataFrame:
    """Build a samples DataFrame from generated rows.

    Args:
        rows: List of row dicts with sample, tumor_bam, normal_bam, analysis_type keys.

    Returns:
        DataFrame with columns: sample, tumor_bam, normal_bam, analysis_type.
        Sorted by sample name.
    """
    if not rows:
        return pd.DataFrame(columns=["sample", "tumor_bam", "normal_bam", "analysis_type"])

    df = pd.DataFrame(rows)[["sample", "tumor_bam", "normal_bam", "analysis_type"]]
    return df.sort_values("sample").reset_index(drop=True)


def print_bam_table(bam_entries: list[dict]) -> None:
    """Pretty-print BAM discovery table with sizes and role guesses.

    Args:
        bam_entries: List of BAM entry dicts with basename, size_bytes, role, reason keys.
    """
    print(f"\nFound {len(bam_entries)} BAMs:\n")
    # Calculate column widths
    max_name = max((len(e["basename"]) for e in bam_entries), default=10)
    max_size = max((len(format_size(e["size_bytes"])) for e in bam_entries), default=6)

    fmt = f"  {{num:>3}}  {{name:<{max_name}}}  {{size:>{max_size}}}   {{role:<8}} {{reason}}"
    for i, entry in enumerate(bam_entries, 1):
        role = entry.get("role", "unknown")
        reason = entry.get("reason", "")
        display_role = role if role != "unknown" else "?"
        print(
            fmt.format(
                num=i,
                name=entry["basename"],
                size=format_size(entry["size_bytes"]),
                role=display_role,
                reason=f"({reason})" if reason else "",
            )
        )
    print()


def print_samples_table(df: pd.DataFrame) -> None:
    """Pretty-print the proposed/final samples.tsv with summary counts.

    Args:
        df: Samples DataFrame.
    """
    if df.empty:
        print("  (no samples)")
        return

    print("\nProposed samples.tsv:\n")

    # Column widths
    col_widths: dict[str, int] = {}
    for col in df.columns:
        max_val = df[col].astype(str).str.len().max()
        col_widths[col] = max(len(col), max_val)

    # Header
    header = "  #   " + "  ".join(col.ljust(col_widths[col]) for col in df.columns)
    print(header)

    # Rows
    for i, (_, row) in enumerate(df.iterrows(), 1):
        line = f"  {i:<3} " + "  ".join(
            str(row[col]).ljust(col_widths[col]) for col in df.columns
        )
        print(line)

    # Summary
    counts = df["analysis_type"].value_counts()
    summary_parts = [f"{count} {atype}" for atype, count in counts.items()]
    print(f"\n  Summary: {', '.join(summary_parts)} ({len(df)} rows)")
    print()


def write_samples_tsv(
    df: pd.DataFrame,
    output_path: Path,
    dry_run: bool = False,
    force: bool = False,
) -> None:
    """Write the samples DataFrame to a TSV file.

    Args:
        df: Samples DataFrame.
        output_path: Path for the output TSV.
        dry_run: If True, do not write.
        force: If True, overwrite without asking.
    """
    if dry_run:
        print(f"[dry-run] Would write: {output_path} ({len(df)} rows)")
        return

    if output_path.is_file() and not force:
        response = input(f"Output file already exists: {output_path}. Overwrite? [y/N] ")
        if response.lower() not in ("y", "yes"):
            print(f"Skipped writing {output_path}")
            return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    print(f"Written: {output_path} ({len(df)} rows)")


# ---------------------------------------------------------------------------
# Section G: Interactive + CLI
# ---------------------------------------------------------------------------


def _prompt(prompt: str, default: str = "") -> str:
    """Prompt the user for input with an optional default value."""
    if default:
        result = input(f"{prompt} [{default}]: ").strip()
        return result if result else default
    return input(f"{prompt}: ").strip()


def _prompt_yn(prompt: str, default: bool = True) -> bool:
    """Prompt for yes/no with a default."""
    suffix = "[Y/n]" if default else "[y/N]"
    result = input(f"{prompt} {suffix}: ").strip().lower()
    if not result:
        return default
    return result in ("y", "yes")


def _prompt_path(prompt: str, default: str = "", must_exist: bool = True) -> str:
    """Prompt for a filesystem path, re-prompting on invalid input."""
    while True:
        raw = _prompt(prompt, default)
        if not raw:
            if not must_exist:
                return ""
            print("  Path cannot be empty. Please try again.")
            continue
        p = Path(raw).expanduser()
        if must_exist and not p.exists():
            print(f"  Path does not exist: {p}")
            retry = _prompt_yn("  Try again?", default=True)
            if not retry:
                return str(p)
            continue
        return str(p)


def _prompt_choice(prompt: str, choices: list[str], default: str = "") -> str:
    """Prompt with numbered choices.

    Args:
        prompt: Prompt text.
        choices: List of choice strings.
        default: Default choice value.

    Returns:
        Selected choice string.
    """
    for i, choice in enumerate(choices, 1):
        marker = " <-- current" if choice == default else ""
        print(f"    [{i}] {choice}{marker}")
    while True:
        raw = _prompt(f"  {prompt}", default=default)
        # Accept by number
        if raw.isdigit():
            idx = int(raw) - 1
            if 0 <= idx < len(choices):
                return choices[idx]
        # Accept by value
        if raw in choices:
            return raw
        print(f"  Invalid choice. Enter 1-{len(choices)} or a value.")


def prompt_unresolved_roles(bam_entries: list[dict]) -> list[dict]:
    """Prompt user to classify BAMs with unknown roles.

    Args:
        bam_entries: BAM entry dicts. Entries with role='unknown' will be prompted.

    Returns:
        Updated bam_entries list with roles assigned.
    """
    unknowns = [e for e in bam_entries if e.get("role") == "unknown"]
    if not unknowns:
        return bam_entries

    print(f"  {len(unknowns)} BAM(s) need manual role assignment.\n")
    for entry in unknowns:
        size = format_size(entry["size_bytes"])
        while True:
            raw = input(
                f"  {entry['basename']} ({size}) -- T(umor) / N(ormal) / S(kip)? "
            ).strip().upper()
            if raw in ("T", "TUMOR"):
                entry["role"] = "tumor"
                entry["reason"] = "manual"
                break
            if raw in ("N", "NORMAL"):
                entry["role"] = "normal"
                entry["reason"] = "manual"
                break
            if raw in ("S", "SKIP"):
                entry["role"] = "skip"
                entry["reason"] = "skipped"
                break
            print("  Please enter T, N, or S.")

    return bam_entries


def _edit_row(df: pd.DataFrame, row_idx: int, bam_entries: list[dict]) -> pd.DataFrame:
    """Edit a single row in the samples table interactively.

    Args:
        df: Current samples DataFrame.
        row_idx: 0-based row index to edit.
        bam_entries: BAM entries for listing available BAMs.

    Returns:
        Updated DataFrame.
    """
    row = df.iloc[row_idx]
    available_bams = [e["basename"] for e in bam_entries if e.get("role") != "skip"]

    print(
        f"\n  Row {row_idx + 1}: {row['sample']}  "
        f"tumor={row['tumor_bam']}  normal={row['normal_bam']}"
    )

    while True:
        print("    [1] Change tumor BAM")
        print("    [2] Change normal BAM")
        print("    [3] Change analysis type")
        print("    [4] Delete row")
        print("    [5] Done editing this row")
        choice = _prompt("  Select", default="5")

        if choice == "1":
            print("  Available BAMs:")
            for i, bam in enumerate(available_bams, 1):
                print(f"    [{i}] {bam}")
            sel = _prompt("  Select BAM")
            if sel.isdigit() and 1 <= int(sel) <= len(available_bams):
                df.at[df.index[row_idx], "tumor_bam"] = available_bams[int(sel) - 1]
            elif sel in available_bams:
                df.at[df.index[row_idx], "tumor_bam"] = sel
        elif choice == "2":
            print("  Available BAMs (enter '.' for none):")
            for i, bam in enumerate(available_bams, 1):
                print(f"    [{i}] {bam}")
            sel = _prompt("  Select BAM", default=".")
            if sel == ".":
                df.at[df.index[row_idx], "normal_bam"] = "."
            elif sel.isdigit() and 1 <= int(sel) <= len(available_bams):
                df.at[df.index[row_idx], "normal_bam"] = available_bams[int(sel) - 1]
            elif sel in available_bams:
                df.at[df.index[row_idx], "normal_bam"] = sel
        elif choice == "3":
            types = ["tumor_only", "tumor_normal", "germline"]
            new_type = _prompt_choice("Analysis type", types, default=row["analysis_type"])
            df.at[df.index[row_idx], "analysis_type"] = new_type
        elif choice == "4":
            df = df.drop(df.index[row_idx]).reset_index(drop=True)
            print("  Row deleted.")
            break
        elif choice == "5":
            break

    return df


def review_loop(
    df: pd.DataFrame,
    bam_entries: list[dict],
    pairings: list[dict],
    analysis_mode: str,
) -> pd.DataFrame | None:
    """Interactive review loop: Accept / Edit / Mode change / Quit.

    Args:
        df: Proposed samples DataFrame.
        bam_entries: BAM entry dicts.
        pairings: From group_and_pair().
        analysis_mode: Current analysis mode.

    Returns:
        Final DataFrame, or None if user quits.
    """
    while True:
        print_samples_table(df)
        choice = _prompt("  [A]ccept / [E]dit rows / [M]ode change / [Q]uit", default="A")
        choice = choice.upper()

        if choice in ("A", "ACCEPT"):
            return df

        if choice in ("Q", "QUIT"):
            return None

        if choice in ("M", "MODE"):
            modes = ["tumor_only", "tumor_normal", "both", "germline"]
            new_mode = _prompt_choice("Analysis mode", modes, default=analysis_mode)
            analysis_mode = new_mode
            rows = generate_rows(pairings, bam_entries, analysis_mode)
            df = build_samples_dataframe(rows)
            continue

        if choice in ("E", "EDIT"):
            raw = _prompt("  Edit which row(s)? (comma-separated, e.g. 1,3)")
            try:
                indices = [int(x.strip()) - 1 for x in raw.split(",")]
            except ValueError:
                print("  Invalid input.")
                continue
            for idx in sorted(indices, reverse=True):
                if 0 <= idx < len(df):
                    df = _edit_row(df, idx, bam_entries)
                else:
                    print(f"  Row {idx + 1} out of range.")
            continue


def interactive_mode() -> None:
    """Guided wizard for generating sm-calling config files."""
    print()
    print("=" * 60)
    print("  sm-calling -- Config Generator")
    print("=" * 60)
    print()

    # Phase 1: Scan + Propose
    bam_folder = _prompt_path(
        "BAM folder", default="results/exomes/bqsr"
    )
    bam_ext = _prompt("BAM extension", default=DEFAULT_BAM_EXTENSION)

    bam_entries = discover_bam_files(bam_folder, bam_ext)

    if not bam_entries:
        # Check if there are any .bam files at all
        bam_dir = Path(bam_folder)
        if bam_dir.is_dir():
            all_bams = [f for f in bam_dir.iterdir() if f.name.endswith(".bam")]
            if all_bams:
                print(
                    f"\n  No BAMs matching '*{bam_ext}' found, "
                    f"but {len(all_bams)} .bam file(s) exist."
                )
                print("  Try a different extension?")
            else:
                print(f"\n  No BAM files found in {bam_folder}")
        else:
            print(f"\n  Directory does not exist: {bam_folder}")
        return

    # Guess roles
    for entry in bam_entries:
        role, reason = guess_role(entry["basename"])
        entry["role"] = role
        entry["reason"] = reason

    print_bam_table(bam_entries)

    # Prompt for unresolved
    bam_entries = prompt_unresolved_roles(bam_entries)

    # Pair and generate rows
    _patient_groups, pairings = group_and_pair(bam_entries)
    analysis_mode = "both"
    rows = generate_rows(pairings, bam_entries, analysis_mode)
    df = build_samples_dataframe(rows)

    # Phase 2: Review
    df = review_loop(df, bam_entries, pairings, analysis_mode)
    if df is None:
        print("Quit without writing.")
        return

    # Phase 3: Write
    samples_output = _prompt("Output samples.tsv", default="config/samples.tsv")
    gen_config = _prompt_yn("Also generate config.yaml?", default=True)

    caller = "mutect2"
    scatter_mode = "chromosome"
    ref_data: dict[str, Any] = {"genome": "", "build": "GRCh38", "search_log": []}
    gatk_resources: dict[str, Any] = {
        "panel_of_normals": "",
        "af_only_gnomad": "",
        "common_biallelic_gnomad": "",
        "search_log": [],
    }
    config_output = "config/config.yaml"

    if gen_config:
        config_output = _prompt("Output config.yaml path", default="config/config.yaml")
        caller = _prompt_choice(
            "Caller", ["mutect2", "freebayes", "all"], default="mutect2"
        )
        scatter_mode = _prompt_choice(
            "Scatter mode", ["chromosome", "interval", "none"], default="chromosome"
        )

        # Reference data discovery
        print("\nScanning for reference data...")
        project_root = Path(bam_folder).resolve().parent
        ref_dir_str = _prompt_path(
            "  Reference data directory (leave empty to auto-scan)",
            must_exist=False,
        )
        ref_dir = Path(ref_dir_str).resolve() if ref_dir_str else None
        ref_data = discover_reference_data(ref_dir, project_root)
        for line in ref_data.get("search_log", []):
            print(line)

        # GATK resources
        print("\nScanning for GATK resources...")
        gatk_dir_str = _prompt_path(
            "  GATK resource directory (leave empty to auto-scan)",
            must_exist=False,
        )
        gatk_dir = Path(gatk_dir_str).resolve() if gatk_dir_str else None
        gatk_resources = discover_gatk_resources(gatk_dir, project_root)
        for line in gatk_resources.get("search_log", []):
            print(line)

    # Write files
    output_path = Path(samples_output)
    write_samples_tsv(df, output_path)

    if gen_config:
        output_folder = _prompt(
            "Output folder for pipeline results",
            default="results/exomes/variant_calls",
        )
        generate_config_template(
            config_output=Path(config_output),
            caller=caller,
            ref_data=ref_data,
            gatk_resources=gatk_resources,
            bam_folder=bam_folder,
            output_folder=output_folder,
            samples_path=samples_output,
            bam_extension=bam_ext,
            scatter_mode=scatter_mode,
        )

    print("\nDone.")


def main() -> None:
    """Main entry point: parse arguments, discover files, generate output."""
    parser = argparse.ArgumentParser(
        description="Generate config/samples.tsv and config/config.yaml for sm-calling.\n"
        "Run without arguments for an interactive guided wizard.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  %(prog)s                                                  # interactive wizard
  %(prog)s --bam-folder /data/bams                          # auto-detect + write
  %(prog)s --bam-folder /data/bams --dry-run                # preview only
  %(prog)s --bam-folder /data/bams --config-template --caller mutect2
  %(prog)s --bam-folder /data/bams --analysis-mode tumor_only --force
        """,
    )
    parser.add_argument(
        "--bam-folder",
        help="Directory containing input BAM files (triggers non-interactive mode)",
    )
    parser.add_argument(
        "--bam-extension",
        default=DEFAULT_BAM_EXTENSION,
        help=f"BAM file extension (default: {DEFAULT_BAM_EXTENSION})",
    )
    parser.add_argument(
        "--analysis-mode",
        choices=["auto", "tumor_only", "tumor_normal", "both", "germline"],
        default="both",
        help="Analysis mode (default: both). 'auto' tries heuristic pairing.",
    )
    parser.add_argument(
        "--caller",
        choices=["mutect2", "freebayes", "all"],
        default="mutect2",
        help="Variant caller for config.yaml (default: mutect2)",
    )
    parser.add_argument(
        "--output",
        default="config/samples.tsv",
        help="Output samples.tsv path (default: config/samples.tsv)",
    )
    parser.add_argument(
        "--config-template",
        action="store_true",
        help="Also generate a config.yaml file",
    )
    parser.add_argument(
        "--config-output",
        default="config/config.yaml",
        help="Output config.yaml path (default: config/config.yaml)",
    )
    parser.add_argument(
        "--ref-dir",
        help="Reference genome directory",
    )
    parser.add_argument(
        "--gatk-resource-dir",
        help="GATK resource VCF directory",
    )
    parser.add_argument(
        "--scatter-mode",
        choices=["chromosome", "interval", "none"],
        default="chromosome",
        help="Scatter mode (default: chromosome)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show output without writing files",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing files without asking",
    )

    args = parser.parse_args()

    # No --bam-folder -> interactive mode
    if not args.bam_folder:
        interactive_mode()
        return

    bam_folder = Path(args.bam_folder).resolve()

    # Discover BAMs
    bam_entries = discover_bam_files(str(bam_folder), args.bam_extension)

    if not bam_entries:
        bam_dir = Path(bam_folder)
        if bam_dir.is_dir():
            all_bams = [f for f in bam_dir.iterdir() if f.name.endswith(".bam")]
            if all_bams:
                sys.exit(
                    f"Error: No BAMs matching '*{args.bam_extension}' in {bam_folder}\n"
                    f"Found {len(all_bams)} .bam file(s). Try --bam-extension."
                )
        sys.exit(f"Error: No BAM files found in {bam_folder}")

    # Guess roles
    for entry in bam_entries:
        role, reason = guess_role(entry["basename"])
        entry["role"] = role
        entry["reason"] = reason

    # Check for unknowns in CLI mode
    unknowns = [e for e in bam_entries if e["role"] == "unknown"]
    if unknowns and args.analysis_mode != "germline":
        unknown_names = ", ".join(e["basename"] for e in unknowns)
        print(
            f"Warning: {len(unknowns)} BAM(s) with unrecognized role: {unknown_names}",
            file=sys.stderr,
        )
        print("  These will be skipped. Use interactive mode for manual assignment.", file=sys.stderr)
        # In CLI mode, skip unknowns
        for e in unknowns:
            e["role"] = "skip"

    # Analysis mode
    analysis_mode = args.analysis_mode
    if analysis_mode == "auto":
        analysis_mode = "both"

    # Pair and generate
    _patient_groups, pairings = group_and_pair(bam_entries)

    # Check for ambiguous pairings in CLI auto mode
    if args.analysis_mode == "auto":
        tumors = [e for e in bam_entries if e["role"] == "tumor"]
        normals = [e for e in bam_entries if e["role"] == "normal"]
        if len(tumors) > 1 and len(normals) > 1:
            # Check if pairing is ambiguous (multiple normals available per tumor)
            for pairing in pairings:
                if pairing["normal"] is None and normals:
                    sys.exit(
                        "Error: Ambiguous tumor-normal pairing detected.\n"
                        "Use interactive mode or specify --analysis-mode explicitly."
                    )

    rows = generate_rows(pairings, bam_entries, analysis_mode)
    df = build_samples_dataframe(rows)

    print_bam_table(bam_entries)
    print_samples_table(df)

    # Write samples.tsv
    write_samples_tsv(df, Path(args.output), dry_run=args.dry_run, force=args.force)

    # Generate config template
    if args.config_template:
        ref_dir = Path(args.ref_dir).resolve() if args.ref_dir else None
        gatk_dir = Path(args.gatk_resource_dir).resolve() if args.gatk_resource_dir else None
        project_root = bam_folder.parent

        print("Scanning for reference data...")
        ref_data = discover_reference_data(ref_dir, project_root)
        for line in ref_data.get("search_log", []):
            print(line)

        print("Scanning for GATK resources...")
        gatk_resources = discover_gatk_resources(gatk_dir, project_root)
        for line in gatk_resources.get("search_log", []):
            print(line)
        print()

        generate_config_template(
            config_output=Path(args.config_output),
            caller=args.caller,
            ref_data=ref_data,
            gatk_resources=gatk_resources,
            bam_folder=str(bam_folder),
            output_folder="results/exomes/variant_calls",
            samples_path=args.output,
            bam_extension=args.bam_extension,
            scatter_mode=args.scatter_mode,
            dry_run=args.dry_run,
            force=args.force,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
