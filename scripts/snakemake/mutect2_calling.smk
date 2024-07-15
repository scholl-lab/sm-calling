import os
import functools
import csv
import yaml

# ----------------------------------------------------------------------------------- #
# Load configuration file containing user-defined settings
configfile: "config.yaml"

# Define temporary directory using an environment variable (usually set by the cluster scheduler)
SCRATCH_DIR = os.environ.get('TMPDIR')

# ----------------------------------------------------------------------------------- #
# Extract user-defined input and output directories and reference file from the configuration file
with open("config.yaml", "r") as stream:
    config = yaml.safe_load(stream)

INPUT_DIR = config["final_bam_folder"]
FINAL_BAM_FILE_EXTENSION = config.get("final_bam_file_extension", ".bam")
OUTPUT_DIR = config["output_folder"]
REFERENCE_FILE = config["reference_unpacked"]
PANEL_OF_NORMALS = config["panel_of_normals"]
AF_ONLY_GNOMAD = config["af_only_gnomad"]

# ----------------------------------------------------------------------------------- #
# Read metadata table into a dictionary from the TSV file
metadata_dict = {}
with open('calling_metadata.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        analysis_key = f"{row['individual1']}_{row['analysis']}"
        metadata_dict[analysis_key] = row

# ----------------------------------------------------------------------------------- #
# Define result directories using functools.partial to join paths with the output folder
prefix_results = functools.partial(os.path.join, OUTPUT_DIR)
VARIANT_DIR = prefix_results('variant_calls')
LOG_DIR = prefix_results('logs')

# List of chromosomes to loop through
chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

# ----------------------------------------------------------------------------------- #
# Helper functions

# Function to return memory based on the number of threads
def get_mem_from_threads(wildcards, threads):
    return threads * 4400

# ----------------------------------------------------------------------------------- #
# Define the rules

rule all:
    input:
        expand(
            f"{VARIANT_DIR}/{{analysis_key}}_{{chromosome}}.vcf.gz", 
            analysis_key=metadata_dict.keys(),
            chromosome=chromosomes if config.get('mutect_scatter_by_chromosome', False) else ['all']
        )

rule call_variants:
    input:
        bam1 = lambda wildcards: os.path.join(INPUT_DIR, metadata_dict[wildcards.analysis_key]['bam1_file_basename'] + FINAL_BAM_FILE_EXTENSION),
        bam2 = lambda wildcards: os.path.join(INPUT_DIR, metadata_dict[wildcards.analysis_key].get('bam2_file_basename', '') + (FINAL_BAM_FILE_EXTENSION if metadata_dict[wildcards.analysis_key].get('bam2_file_basename') else '')),
    output:
        variant_file = f"{VARIANT_DIR}/{{analysis_key}}_{{chromosome}}.vcf.gz",
    threads: 4  # Adjust as necessary
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        "gatk"
    params:
        reference = REFERENCE_FILE,
        normal_sample = lambda wildcards: metadata_dict[wildcards.analysis_key].get('sample2', ''),
        individual = lambda wildcards: metadata_dict[wildcards.analysis_key]['individual1'],
        analysis = lambda wildcards: metadata_dict[wildcards.analysis_key]['analysis'],
        panel_of_normals = PANEL_OF_NORMALS,
        af_only_gnomad = AF_ONLY_GNOMAD
    log:
        mutect2 = f"{LOG_DIR}/{{analysis_key}}_{{chromosome}}.mutect2.log"
    shell:
        """
        # Determine if bam2 file is provided and set the respective options
        bam2_option=""
        normal_sample_option=""
        if [ ! -z "{input.bam2}" ] && [ "{input.bam2}" != "{INPUT_DIR}" ]; then
            bam2_option="-I {input.bam2}"
            normal_sample_option="-normal {params.normal_sample}"
        fi
        
        # Determine if scattering by chromosome is enabled and set the respective option
        scatter_option=""
        if [ "{config[mutect_scatter_by_chromosome]}" = "True" ] && [ "{wildcards.chromosome}" != "all" ]; then
            scatter_option="-L {wildcards.chromosome}"
        fi
        
        gatk --java-options '-Xms4000m -Xmx10g -Djava.io.tmpdir={resources.tmpdir}' Mutect2 \
            -R {params.reference} \
            -I {input.bam1} \
            $bam2_option \
            $normal_sample_option \
            --germline-resource {params.af_only_gnomad} \
            --panel-of-normals {params.panel_of_normals} \
            --f1r2-tar-gz {VARIANT_DIR}/{params.individual}_{params.analysis}_{wildcards.chromosome}.f1r2.tar.gz \
            $scatter_option \
            -O {output.variant_file} 2> {log.mutect2}
        """

# ----------------------------------------------------------------------------------- #
