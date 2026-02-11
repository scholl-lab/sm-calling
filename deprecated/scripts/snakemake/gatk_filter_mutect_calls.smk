import os
import functools
import csv
import yaml

# ----------------------------------------------------------------------------------- #
# Script Description:
# This script uses Snakemake to manage a workflow that filters merged raw MuTect2 
# variant calls for each individual and analysis type. It reads a metadata file to get 
# the necessary details about the files to process. The GATK tool is used for filtering 
# and indexing the resulting VCF files.
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Load configuration file containing user-defined settings
configfile: "config.yaml"

# Define temporary directory using an environment variable
SCRATCH_DIR = os.environ.get('TMPDIR')
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Helper Functions:
# In this section, we define helper functions that assist in various tasks throughout
# the script. 
#
# - get_mem_from_threads: Calculates the amount of memory to allocate based on the 
#                         number of threads.

def get_mem_from_threads(wildcards, threads):
    """Calculate the amount of memory to allocate based on the number of threads."""
    return threads * 4400  # Adjust as necessary
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Extract Input, Output Directories and Reference File:
# Here, we extract user-defined settings from the configuration file.

with open("config.yaml", "r") as stream:
    config = yaml.safe_load(stream)

OUTPUT_DIR = config["output_folder"]
REFERENCE_FILE = config["reference_unpacked"]
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Metadata Loading:
# We load the metadata file which contains details about the individuals and analyses 
# to be processed. The metadata file is expected to contain columns 'individual1' and 
# 'analysis', outlining the details required to process each file.

metadata_dict = {}
with open('calling_metadata.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        analysis_key = f"{row['individual1']}_{row['analysis']}"
        metadata_dict[analysis_key] = row
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Directory Setup:
# Here we define the directories where the results will be stored and ensure that 
# these directories exist.

prefix_results = functools.partial(os.path.join, OUTPUT_DIR)
FILTERED_DIR = prefix_results('filtered_vcfs')
LOG_DIR = prefix_results('logs')

# Ensure log directory exists
os.makedirs(LOG_DIR, exist_ok=True)
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Snakemake Rules:
# Below are the Snakemake rules that define the workflow. The "all" rule specifies the
# desired final outputs, effectively "driving" the workflow. The "filter_mutect_calls" 
# rule performs the core task of filtering VCF files.

# Rule to specify the final output files
rule all:
    input:
        expand(
            os.path.join(FILTERED_DIR, "{analysis_key}.filtered.vcf.gz"),
            analysis_key=metadata_dict.keys()
        )

# Rule to filter merged raw MuTect2 calls
rule filter_mutect_calls:
    input:
        vcf = lambda wildcards: os.path.join(OUTPUT_DIR, 'variant_merge', f"{wildcards.analysis_key}.vcf.gz"),
        stats = lambda wildcards: os.path.join(OUTPUT_DIR, 'variant_merge', f"{wildcards.analysis_key}.vcf.gz.stats"),
        contamination_table = lambda wildcards: os.path.join(OUTPUT_DIR, 'calculate_contamination', f"{wildcards.analysis_key}.contamination.table"),
        tumor_segmentation = lambda wildcards: os.path.join(OUTPUT_DIR, 'calculate_contamination', f"{wildcards.analysis_key}.segments.table"),
        ob_priors = lambda wildcards: os.path.join(OUTPUT_DIR, 'variant_merge', f"{wildcards.analysis_key}_read-orientation-model.tar.gz")
    output:
        filtered_vcf = os.path.join(FILTERED_DIR, "{analysis_key}.filtered.vcf.gz")
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        os.path.join(LOG_DIR, "{analysis_key}.filter_mutect_calls.log")
    shell:
        """
        gatk --java-options '-Xms4000m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={resources.tmpdir}' FilterMutectCalls \
        -R {REFERENCE_FILE} \
        -V {input.vcf} \
        --stats {input.stats} \
        --contamination-table {input.contamination_table} \
        --tumor-segmentation {input.tumor_segmentation} \
        --ob-priors {input.ob_priors} \
        -O {output.filtered_vcf} 2> {log}
        """
# ----------------------------------------------------------------------------------- #
