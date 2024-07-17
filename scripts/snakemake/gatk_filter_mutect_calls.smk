import os
import functools
import csv

# ----------------------------------------------------------------------------------- #
# Configuration and Environment Setup

# Load configuration file containing user-defined settings
configfile: "config.yaml"

# Define temporary directory using an environment variable (usually set by the cluster scheduler)
SCRATCH_DIR = os.environ.get('TMPDIR')

# Helper function to return memory based on the number of threads
def get_mem_from_threads(wildcards, threads):
    return threads * 4400
# ----------------------------------------------------------------------------------- #

# Extract user-defined output directory from the configuration file
OUTPUT_DIR = config["output_folder"]
REFERENCE_FILE = config["reference_unpacked"]

# Load metadata from the provided TSV file and store in a dictionary
metadata_dict = {}
with open('calling_metadata.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        analysis_key = f"{row['individual1']}_{row['analysis']}"
        metadata_dict[analysis_key] = row

# Define output directories for various results using functools
prefix_results = functools.partial(os.path.join, OUTPUT_DIR)
FILTERED_DIR = prefix_results('filtered_vcfs')
LOG_DIR = prefix_results('logs')

# ----------------------------------------------------------------------------------- #

# Main rule to define the expected final outputs
rule all:
    input:
        expand(
            os.path.join(FILTERED_DIR, "{individual1}_{analysis}.filtered.vcf.gz"),
            individual1=[row['individual1'] for row in metadata_dict.values()],
            analysis=[row['analysis'] for row in metadata_dict.values()]
        )

# Rule to filter merged raw MuTect2 calls
rule filter_mutect_calls:
    input:
        vcf = lambda wildcards: os.path.join(OUTPUT_DIR, 'variant_merge', f"{wildcards.individual1}_{wildcards.analysis}.vcf.gz"),
        stats = lambda wildcards: os.path.join(OUTPUT_DIR, 'variant_merge', f"{wildcards.individual1}_{wildcards.analysis}.vcf.gz.stats"),
        contamination_table = lambda wildcards: os.path.join(OUTPUT_DIR, 'calculate_contamination', f"{wildcards.individual1}_{wildcards.analysis}.contamination.table"),
        tumor_segmentation = lambda wildcards: os.path.join(OUTPUT_DIR, 'calculate_contamination', f"{wildcards.individual1}_{wildcards.analysis}.segments.table"),
        ob_priors = lambda wildcards: os.path.join(OUTPUT_DIR, 'variant_merge', f"{wildcards.individual1}_{wildcards.analysis}_read-orientation-model.tar.gz")
    output:
        filtered_vcf = os.path.join(FILTERED_DIR, "{individual1}_{analysis}.filtered.vcf.gz")
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    log:
        os.path.join(LOG_DIR, "{individual1}_{analysis}.filter_mutect_calls.log")
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
