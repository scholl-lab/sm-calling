import os
import functools
import csv

# ----------------------------------------------------------------------------------- #
# Load configuration file containing user-defined settings
configfile: "config.yaml"

# Define temporary directory using an environment variable (usually set by the cluster scheduler)
SCRATCH_DIR = os.environ.get('TMPDIR')

# Helper function to return memory based on the number of threads
def get_mem_from_threads(wildcards, threads):
    """Calculate the amount of memory to allocate based on the number of threads."""
    return threads * 4400
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Extract user-defined input and output directories and reference (unpacked) file from the configuration file
INPUT_DIR = config["final_bam_folder"]
FINAL_BAM_FILE_EXTENSION = config.get("final_bam_file_extension", ".bam")
REFERENCE_FILE = config["reference_unpacked"]
COMMON_BIALLELIC_GNOMAD = config["common_biallelic_gnomad"]

# Read metadata table into a dictionary from the TSV file
metadata_dict = {}
with open('calling_metadata.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        analysis_key = f"{row['individual1']}_{row['analysis']}"
        metadata_dict[analysis_key] = row
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define result directories using functools.partial to join paths with the output folder
prefix_results = functools.partial(os.path.join, config['output_folder'])
VARIANT_DIR = prefix_results('variant_calls')
CONTAMINATION_DIR = prefix_results('calculate_contamination')
LOG_DIR = prefix_results('logs')

# Set of all unique BAM basenames
all_bams = set(row['bam1_file_basename'] for row in metadata_dict.values()).union(
           set(row['bam2_file_basename'] for row in metadata_dict.values()))

# Ensure log directory exists
os.makedirs(LOG_DIR, exist_ok=True)
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Snakemake Rules

# Main rule to define the expected final outputs
rule all:
    input:
        expand(
            f"{CONTAMINATION_DIR}/{{analysis_key}}.contamination.table",
            analysis_key=metadata_dict.keys()
        )

# GetPileupSummaries rule
rule get_pileup_summaries:
    input:
        bam = lambda wildcards: os.path.join(INPUT_DIR, wildcards.bam_basename + FINAL_BAM_FILE_EXTENSION),
        ref = REFERENCE_FILE,
        gnomad = COMMON_BIALLELIC_GNOMAD
    output:
        pileup_table = os.path.join(CONTAMINATION_DIR, "{bam_basename}.getpileupsummaries.table")
    threads: 4  # Adjust as necessary
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        os.path.join(LOG_DIR, "{bam_basename}.getpileupsummaries.log")
    shell:
        """
        gatk --java-options '-Xms4000m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={resources.tmpdir}' GetPileupSummaries \
            -I {input.bam} \
            -R {input.ref} \
            -V {input.gnomad} \
            -L {input.gnomad} \
            -O {output.pileup_table} 2> {log}
        """

# CalculateContamination rule
rule calculate_contamination:
    input:
        tumor_pileup = lambda wildcards: os.path.join(CONTAMINATION_DIR, metadata_dict[wildcards.analysis_key]['bam1_file_basename'] + ".getpileupsummaries.table"),
        normal_pileup = lambda wildcards: os.path.join(CONTAMINATION_DIR, metadata_dict[wildcards.analysis_key]['bam2_file_basename'] + ".getpileupsummaries.table") if metadata_dict[wildcards.analysis_key].get('bam2_file_basename') else []
    output:
        contamination_table = os.path.join(CONTAMINATION_DIR, "{analysis_key}.contamination.table"),
        segments_table = os.path.join(CONTAMINATION_DIR, "{analysis_key}.segments.table")
    params:
        matched_option = lambda wildcards: f"--matched {os.path.join(CONTAMINATION_DIR, metadata_dict[wildcards.analysis_key]['bam2_file_basename'] + '.getpileupsummaries.table')}" if metadata_dict[wildcards.analysis_key].get('bam2_file_basename') else ''
    threads: 4  # Adjust as necessary
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        os.path.join(LOG_DIR, "{analysis_key}.calculatecontamination.log")
    shell:
        """
        gatk --java-options '-Xms4000m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={resources.tmpdir}' CalculateContamination \
            -I {input.tumor_pileup} \
            {params.matched_option} \
            --tumor-segmentation {output.segments_table} \
            -O {output.contamination_table} 2> {log}
        """
# ----------------------------------------------------------------------------------- #
