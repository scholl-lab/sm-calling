import os
import functools
import csv

# ----------------------------------------------------------------------------------- #
# Script Description:
# This script uses Snakemake to manage a workflow that concatenates scattered VCF files
# into single VCF files per individual and analysis type. The script reads a metadata file
# to get the necessary details about the files to process. It leverages the bcftools
# package to perform the concatenation and index the resulting files. 

# The script now supports a flexible definition of chromosomal intervals, allowing for 
# both numbered and named intervals (e.g., 'X', 'Y') as well as defined ranges (e.g., '1-22').
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
# - expand_interval_names: This function effectively handles both individual 
#                          chromosome/interval names and ranges, translating them into 
#                          a list of interval names used in the workflow.

def get_mem_from_threads(wildcards, threads):
    """Calculate the amount of memory to allocate based on the number of threads."""
    return threads * 4400  # Adjust as necessary

def expand_interval_names(interval_ranges, prefix=""):
    """
    Expand a list of interval ranges into a list of interval names.

    Parameters:
    - interval_ranges (list): A list of interval ranges specified as strings.
                              For example: ["1-3", "5", "X-Y"].
    - prefix (str): A prefix to add to each interval name. For example, "chr".

    Returns:
    - list: A list of interval names with the prefix added to each name.
    """
    interval_names = []
    for interval_range in interval_ranges:
        if '-' in str(interval_range):
            start, end = map(int, str(interval_range).split('-'))
            interval_names.extend([f"{prefix}{i}" for i in range(start, end+1)])
        else:
            interval_names.append(f"{prefix}{interval_range}")
    return interval_names
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Configuration Settings and Metadata Loading:
# Here, we extract user-defined settings from the configuration file and load the metadata
# file which contains details about the individuals and analyses to be processed.

# The metadata file is expected to contain columns 'individual1' and 'analysis', 
# outlining the details required to process each file.

# Extract paths and settings from the configuration file
INPUT_DIR = config["scattered_vcf_folder"]
OUTPUT_DIR = config["output_folder"]
METADATA_FILE = config["metadata_file"]
VCF_FILE_EXTENSION = config.get("vcf_file_extension", ".vcf.gz")

# Get the prefix for scatter names; default is empty
SCATTER_NAME_PREFIX = config.get("scatter_name_prefix", "")

# Define the scattering interval names based on the ranges and prefixes specified in the
# configuration file. This allows for a flexible definition of intervals, accommodating 
# various naming conventions and range specifications.
SCATTERING_INTERVAL_RANGES = config.get("scattering_interval_ranges", [str(i) for i in range(1, 23)] + ['X', 'Y'])
SCATTERING_INTERVAL_NAMES = expand_interval_names(SCATTERING_INTERVAL_RANGES, SCATTER_NAME_PREFIX)

# Get the delimiter for scatter names; default is "."
SCATTER_NAME_DELIMITER = config.get("scatter_name_delimiter", ".")

# Read the metadata file to build a dictionary with analysis details
metadata_dict = {}
with open(METADATA_FILE, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        analysis_key = f"{row['individual1']}_{row['analysis']}"
        metadata_dict[analysis_key] = row
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Directory Setup:
# Here we define the directories where the results will be stored, and ensure that these
# directories exist.

# Note: Given the flexibility in defining scatter name prefixes and interval names, 
# it is important to ensure that the input file names follow the convention defined 
# in the configuration file.

# Define result directories
prefix_results = functools.partial(os.path.join, config['output_folder'])
CONCATENATED_VCF_DIR = prefix_results('concatenated_vcfs')
LOG_DIR = prefix_results('logs')

# Ensure output directories exist
os.makedirs(CONCATENATED_VCF_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Snakemake Rules:
# Below are the Snakemake rules that define the workflow. The "all" rule specifies the
# desired final outputs, effectively "driving" the workflow. The "concatenate_vcfs" rule
# performs the core task of concatenating VCF files and indexing the results.

# Rule to specify the final output files
rule all:
    input:
        expand(
            f"{CONCATENATED_VCF_DIR}/{{individual1}}_{{analysis}}{VCF_FILE_EXTENSION}", 
            individual1=[row['individual1'] for row in metadata_dict.values()], 
            analysis=[row['analysis'] for row in metadata_dict.values()]
        )

# Rule to concatenate VCF files and index the resulting file
rule concatenate_vcfs:
    input:
        lambda wildcards: expand(f"{INPUT_DIR}/{wildcards.individual1}_{wildcards.analysis}{SCATTER_NAME_DELIMITER}{{interval}}{VCF_FILE_EXTENSION}", interval=SCATTERING_INTERVAL_NAMES)
    output:
        concatenated_vcf = f"{CONCATENATED_VCF_DIR}/{{individual1}}_{{analysis}}{VCF_FILE_EXTENSION}"
    threads: 2  # Adjust as necessary
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        "bcftools"
    log:
        f"{LOG_DIR}/{{individual1}}_{{analysis}}_concat.log"
    shell:
        """
        echo "Starting concatenate_vcfs at: $(date)" >> {log}
        bcftools concat {input} -o {output.concatenated_vcf} -Oz 2>> {log}
        bcftools index --threads {threads} -t {output.concatenated_vcf} 2>> {log}
        echo "Finished concatenate_vcfs at: $(date)" >> {log}
        """
# ----------------------------------------------------------------------------------- #
