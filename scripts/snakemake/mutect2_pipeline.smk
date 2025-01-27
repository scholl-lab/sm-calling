##############################################################################
# mutect2_pipeline.smk
##############################################################################
import os
import functools
import csv
import yaml

##############################################################################
# 1) Load the config file
##############################################################################
configfile: "config.yaml"

# We'll assume your config.yaml has these top-level keys (edit if needed):
#   final_bam_folder
#   output_folder
#   reference_unpacked
#   panel_of_normals
#   af_only_gnomad
#   common_biallelic_gnomad
#   final_bam_file_extension
#   mutect_scatter_by_chromosome
#   log_dir_sub
#   metadata_file

# Extract user-defined paths and references
FINAL_BAM_FOLDER         = config["final_bam_folder"]
OUTPUT_FOLDER            = config["output_folder"]
REFERENCE_UNPACKED       = config["reference_unpacked"]
PANEL_OF_NORMALS         = config["panel_of_normals"]
AF_ONLY_GNOMAD           = config["af_only_gnomad"]
COMMON_BIALLELIC_GNOMAD  = config["common_biallelic_gnomad"]
FINAL_BAM_EXTENSION      = config.get("final_bam_file_extension", ".bam")
SCATTER_BY_CHROM         = config.get("mutect_scatter_by_chromosome", False)
LOG_SUBFOLDER            = config["log_dir_sub"]
METADATA_FILE            = config["metadata_file"]

##############################################################################
# 2) Optional environment-based scratch directory (e.g. for cluster jobs)
##############################################################################
SCRATCH_DIR = os.environ.get("TMPDIR", "/tmp")

##############################################################################
# 3) Read the metadata (calling_metadata.tsv or user-defined)
##############################################################################
metadata_dict = {}
print(f"DEBUG: Loading metadata from {METADATA_FILE}")

with open(METADATA_FILE, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        # e.g. 'analysis_key' = individual_analysis
        analysis_key = f"{row['individual1']}_{row['analysis']}"
        metadata_dict[analysis_key] = row

print(f"DEBUG: Number of entries in metadata: {len(metadata_dict)}")

# If scattering by chromosome is turned on, define full list of contigs
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

##############################################################################
# 4) Define helper functions for resource usage and expansions
##############################################################################
def get_mem_from_threads(wildcards, threads):
    """Default: 4400 MB per thread."""
    return threads * 4400

# For referencing final scattered or single-chrom outputs
def get_chrom_list():
    return CHROMOSOMES if SCATTER_BY_CHROM else ["all"]

##############################################################################
# 5) Make sure directories exist
##############################################################################
VARIANT_DIR        = os.path.join(OUTPUT_FOLDER, "variant_calls")
MERGED_DIR         = os.path.join(OUTPUT_FOLDER, "variant_merge")
CONTAM_DIR         = os.path.join(OUTPUT_FOLDER, "calculate_contamination")
FILTERED_DIR       = os.path.join(OUTPUT_FOLDER, "filtered_vcfs")
LOG_DIR            = os.path.join(OUTPUT_FOLDER, LOG_SUBFOLDER)

# Create them if they do not exist
for d in [VARIANT_DIR, MERGED_DIR, CONTAM_DIR, FILTERED_DIR, LOG_DIR]:
    os.makedirs(d, exist_ok=True)

##############################################################################
# 6) Final outputs & "rule all"
#    Our final deliverable is the "filtered" VCF for each analysis_key
##############################################################################
rule all:
    """
    Produce final Filtered Mutect2 VCFs for each (individual, analysis) pair.
    """
    input:
        # a) Ensure scattered calls exist
        expand(
            f"{VARIANT_DIR}/{{analysis_key}}_{{chrom}}.vcf.gz",
            analysis_key=metadata_dict.keys(),
            chrom=get_chrom_list()
        ),
        # b) Ensure merged calls (VCF, TBI, stats, f1r2 tar)
        expand(
            [
                f"{MERGED_DIR}/{{analysis_key}}.vcf.gz",
                f"{MERGED_DIR}/{{analysis_key}}.vcf.gz.tbi",
                f"{MERGED_DIR}/{{analysis_key}}.vcf.gz.stats",
                f"{MERGED_DIR}/{{analysis_key}}_read-orientation-model.tar.gz"
            ],
            analysis_key=metadata_dict.keys()
        ),
        # c) Ensure contamination tables
        expand(
            [
                f"{CONTAM_DIR}/{{analysis_key}}.contamination.table",
                f"{CONTAM_DIR}/{{analysis_key}}.segments.table"
            ],
            analysis_key=metadata_dict.keys()
        ),
        # d) Final filtered VCF
        expand(
            f"{FILTERED_DIR}/{{analysis_key}}.filtered.vcf.gz",
            analysis_key=metadata_dict.keys()
        )

##############################################################################
# 7) Call Variants (Mutect2) scattered by chromosome
##############################################################################
rule call_variants:
    """
    Call variants with Mutect2 for each analysis_key, either scattered by chromosome 
    or once with 'all' if SCATTER_BY_CHROM=False.
    """
    input:
        bam1=lambda wc: os.path.join(
            FINAL_BAM_FOLDER,
            metadata_dict[wc.analysis_key]["bam1_file_basename"] + FINAL_BAM_EXTENSION
        ),
        bam2=lambda wc: (
            os.path.join(
                FINAL_BAM_FOLDER,
                metadata_dict[wc.analysis_key]["bam2_file_basename"] + FINAL_BAM_EXTENSION
            )
            if metadata_dict[wc.analysis_key].get("bam2_file_basename")
            else ""
        )
    output:
        variant_file = f"{VARIANT_DIR}/{{analysis_key}}_{{chromosome}}.vcf.gz"
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time   = "72:00:00",
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    params:
        reference     = REFERENCE_UNPACKED,
        pon           = PANEL_OF_NORMALS,
        af_gnomad     = AF_ONLY_GNOMAD,
        normal_sample = lambda wc: metadata_dict[wc.analysis_key].get("sample2", ""),
        individual    = lambda wc: metadata_dict[wc.analysis_key]["individual1"],
        analysis      = lambda wc: metadata_dict[wc.analysis_key]["analysis"]
    log:
        mutect2_log = f"{LOG_DIR}/{{analysis_key}}_{{chromosome}}.mutect2.log"
    shell:
        r"""
        echo "DEBUG: Starting Mutect2 call for {wildcards.analysis_key} {wildcards.chromosome}" >&2

        # Optionally use normal bam if present
        bam2_option=""
        normal_option=""
        if [ -n "{input.bam2}" ] && [ "{input.bam2}" != "" ]; then
            bam2_option="-I {input.bam2}"
            normal_option="-normal {params.normal_sample}"
        fi

        # Scatter by chromosome or do 'all'
        scatter_option=""
        if [ "{SCATTER_BY_CHROM}" = "True" ] && [ "{wildcards.chromosome}" != "all" ]; then
            scatter_option="-L {wildcards.chromosome}"
        fi

        gatk --java-options '-Xms4000m -Xmx10g -Djava.io.tmpdir={resources.tmpdir}' Mutect2 \
            -R "{params.reference}" \
            -I "{input.bam1}" \
            $bam2_option \
            $normal_option \
            --germline-resource "{params.af_gnomad}" \
            --panel-of-normals "{params.pon}" \
            --genotype-germline-sites true \
            --genotype-pon-sites true \
            --f1r2-tar-gz "{VARIANT_DIR}/{params.individual}_{params.analysis}_{wildcards.chromosome}.f1r2.tar.gz" \
            $scatter_option \
            -O "{output.variant_file}" 2> "{log.mutect2_log}"

        echo "DEBUG: Finished Mutect2 call for {wildcards.analysis_key} {wildcards.chromosome}" >&2
        """

##############################################################################
# 8) Merge scattered calls (VCF, stats, f1r2)
##############################################################################
rule merge_vcfs:
    """
    Merge scattered VCF files across all chosen chromosomes into one final VCF per analysis_key.
    """
    input:
        lambda wc: [f"{VARIANT_DIR}/{wc.analysis_key}_{chrom}.vcf.gz" for chrom in get_chrom_list()]
    output:
        vcf  = f"{MERGED_DIR}/{{analysis_key}}.vcf.gz",
        tbi  = f"{MERGED_DIR}/{{analysis_key}}.vcf.gz.tbi"
    params:
        input_files = lambda wc: ' '.join(
            f"-I {VARIANT_DIR}/{wc.analysis_key}_{chrom}.vcf.gz" for chrom in get_chrom_list()
        )
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time   = "72:00:00",
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        merge_log = f"{LOG_DIR}/{{analysis_key}}.merge_vcfs.log"
    shell:
        r"""
        echo "DEBUG: Merging VCFs for {wildcards.analysis_key}" >&2
        gatk --java-options '-Xms4000m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={resources.tmpdir}' GatherVcfs \
            {params.input_files} \
            -O "{output.vcf}" 2> "{log.merge_log}"

        tabix -p vcf "{output.vcf}"
        echo "DEBUG: Finished merging VCFs for {wildcards.analysis_key}" >&2
        """

rule merge_stats:
    """
    Merge .stats files from scattered calls.
    """
    input:
        lambda wc: [f"{VARIANT_DIR}/{wc.analysis_key}_{chrom}.vcf.gz.stats" for chrom in get_chrom_list()]
    output:
        stats_merged = f"{MERGED_DIR}/{{analysis_key}}.vcf.gz.stats"
    params:
        input_files = lambda wc: ' '.join(
            f"-stats {VARIANT_DIR}/{wc.analysis_key}_{chrom}.vcf.gz.stats" for chrom in get_chrom_list()
        )
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time   = "72:00:00",
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        merge_stats_log = f"{LOG_DIR}/{{analysis_key}}.merge_stats.log"
    shell:
        r"""
        echo "DEBUG: Merging stats for {wildcards.analysis_key}" >&2
        gatk --java-options '-Xms4000m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={resources.tmpdir}' MergeMutectStats \
            {params.input_files} \
            -O "{output.stats_merged}" 2> "{log.merge_stats_log}"
        echo "DEBUG: Finished merging stats for {wildcards.analysis_key}" >&2
        """

rule merge_f1r2:
    """
    Merge scattered f1r2 tar.gz files across chromosomes.
    """
    input:
        lambda wc: [f"{VARIANT_DIR}/{wc.analysis_key}_{chrom}.f1r2.tar.gz" for chrom in get_chrom_list()]
    output:
        tar_merged = f"{MERGED_DIR}/{{analysis_key}}_read-orientation-model.tar.gz"
    params:
        input_files = lambda wc: ' '.join(
            f"-I {VARIANT_DIR}/{wc.analysis_key}_{chrom}.f1r2.tar.gz" for chrom in get_chrom_list()
        )
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time   = "72:00:00",
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        merge_f1r2_log = f"{LOG_DIR}/{{analysis_key}}.merge_f1r2.log"
    shell:
        r"""
        echo "DEBUG: Merging f1r2 files for {wildcards.analysis_key}" >&2
        gatk --java-options '-Xms4000m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={resources.tmpdir}' LearnReadOrientationModel \
            {params.input_files} \
            -O "{output.tar_merged}" 2> "{log.merge_f1r2_log}"
        echo "DEBUG: Finished merging f1r2 for {wildcards.analysis_key}" >&2
        """

##############################################################################
# 9) GetPileupSummaries for contamination (for each unique BAM)
##############################################################################
all_bam_basenames = set()
for row in metadata_dict.values():
    if row.get("bam1_file_basename"):
        all_bam_basenames.add(row["bam1_file_basename"])
    if row.get("bam2_file_basename"):
        all_bam_basenames.add(row["bam2_file_basename"])

rule get_pileup_summaries:
    """
    Generate pileup summaries for each unique BAM found in metadata.
    """
    input:
        bam=lambda wc: os.path.join(FINAL_BAM_FOLDER, wc.bam_basename + FINAL_BAM_EXTENSION),
        ref=REFERENCE_UNPACKED,
        gnomad=COMMON_BIALLELIC_GNOMAD
    output:
        pileup_table=f"{CONTAM_DIR}/{{bam_basename}}.getpileupsummaries.table"
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time   = "72:00:00",
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        getps_log = f"{LOG_DIR}/{{bam_basename}}.getpileupsummaries.log"
    shell:
        r"""
        echo "DEBUG: Running GetPileupSummaries on {wildcards.bam_basename}" >&2
        gatk --java-options '-Xms4000m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={resources.tmpdir}' GetPileupSummaries \
            -I "{input.bam}" \
            -R "{input.ref}" \
            -V "{input.gnomad}" \
            -L "{input.gnomad}" \
            -O "{output.pileup_table}" 2> "{log.getps_log}"
        echo "DEBUG: Finished GetPileupSummaries on {wildcards.bam_basename}" >&2
        """

##############################################################################
# 10) CalculateContamination for each analysis_key
##############################################################################
rule calculate_contamination:
    """
    Calculate contamination for each tumor sample (and matched normal, if present).
    """
    input:
        tumor_pileup=lambda wc: os.path.join(
            CONTAM_DIR, metadata_dict[wc.analysis_key]["bam1_file_basename"] + ".getpileupsummaries.table"
        ),
        normal_pileup=lambda wc: (
            os.path.join(
                CONTAM_DIR, metadata_dict[wc.analysis_key]["bam2_file_basename"] + ".getpileupsummaries.table"
            )
            if metadata_dict[wc.analysis_key].get("bam2_file_basename")
            else []
        )
    output:
        contamination_table=f"{CONTAM_DIR}/{{analysis_key}}.contamination.table",
        segments_table=f"{CONTAM_DIR}/{{analysis_key}}.segments.table"
    params:
        matched_option=lambda wc: (
            f"--matched {os.path.join(CONTAM_DIR, metadata_dict[wc.analysis_key]['bam2_file_basename'] + '.getpileupsummaries.table')}"
            if metadata_dict[wc.analysis_key].get("bam2_file_basename")
            else ""
        )
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time   = "72:00:00",
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        calccont_log = f"{LOG_DIR}/{{analysis_key}}.calculatecontamination.log"
    shell:
        r"""
        echo "DEBUG: Calculating contamination for {wildcards.analysis_key}" >&2
        gatk --java-options '-Xms4000m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={resources.tmpdir}' CalculateContamination \
            -I "{input.tumor_pileup}" \
            {params.matched_option} \
            --tumor-segmentation "{output.segments_table}" \
            -O "{output.contamination_table}" 2> "{log.calccont_log}"
        echo "DEBUG: Finished contamination for {wildcards.analysis_key}" >&2
        """

##############################################################################
# 11) Filter Mutect Calls
##############################################################################
rule filter_mutect_calls:
    """
    Filter the merged raw Mutect2 calls using FilterMutectCalls.
    """
    input:
        vcf=lambda wc: os.path.join(MERGED_DIR, f"{wc.analysis_key}.vcf.gz"),
        stats=lambda wc: os.path.join(MERGED_DIR, f"{wc.analysis_key}.vcf.gz.stats"),
        contamination_table=lambda wc: os.path.join(CONTAM_DIR, f"{wc.analysis_key}.contamination.table"),
        tumor_segmentation=lambda wc: os.path.join(CONTAM_DIR, f"{wc.analysis_key}.segments.table"),
        ob_priors=lambda wc: os.path.join(MERGED_DIR, f"{wc.analysis_key}_read-orientation-model.tar.gz")
    output:
        filtered_vcf=f"{FILTERED_DIR}/{{analysis_key}}.filtered.vcf.gz"
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time   = "72:00:00",
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        filter_log = f"{LOG_DIR}/{{analysis_key}}.filter_mutect_calls.log"
    shell:
        r"""
        echo "DEBUG: Filtering VCF for {wildcards.analysis_key}" >&2
        gatk --java-options '-Xms4000m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={resources.tmpdir}' FilterMutectCalls \
            -R "{REFERENCE_UNPACKED}" \
            -V "{input.vcf}" \
            --stats "{input.stats}" \
            --contamination-table "{input.contamination_table}" \
            --tumor-segmentation "{input.tumor_segmentation}" \
            --ob-priors "{input.ob_priors}" \
            -O "{output.filtered_vcf}" 2> "{log.filter_log}"

        echo "DEBUG: Finished filtering VCF for {wildcards.analysis_key}" >&2
        """
