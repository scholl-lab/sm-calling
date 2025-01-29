##############################################################################
# mutect2_calling.smk
##############################################################################
import os
import csv
import yaml
import functools

##############################################################################
# 1) Load config and environment
##############################################################################
configfile: "config.yaml"

# Required config keys (existing from your pipeline):
FINAL_BAM_FOLDER       = config["final_bam_folder"]
METADATA_FILE          = config["metadata_file"]
OUTPUT_DIR             = config["output_folder"]
REFERENCE_FILE         = config["reference_unpacked"]
PANEL_OF_NORMALS       = config["panel_of_normals"]
AF_ONLY_GNOMAD         = config["af_only_gnomad"]
FINAL_BAM_EXTENSION    = config.get("final_bam_file_extension", ".bam")
LOG_SUBFOLDER          = config["log_dir_sub"]

# Additional config keys for advanced scattering:
SCATTER_MODE           = config.get("scatter_mode", "chromosome")  # "chromosome", "interval", or "none"
SCATTER_COUNT          = config.get("scatter_count", 24)           # e.g. number of intervals
INTERVALS_DIR          = config.get("intervals_dir", os.path.join(OUTPUT_DIR, "intervals"))

# Make sure directories exist
VARIANT_DIR = os.path.join(OUTPUT_DIR, "variant_calls")
LOG_DIR     = os.path.join(OUTPUT_DIR, LOG_SUBFOLDER)
for d in [VARIANT_DIR, LOG_DIR, INTERVALS_DIR]:
    os.makedirs(d, exist_ok=True)

##############################################################################
# 2) Read metadata
##############################################################################
metadata_dict = {}
print(f"DEBUG: Loading metadata from {METADATA_FILE}")
with open(METADATA_FILE, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        analysis_key = f"{row['individual1']}_{row['analysis']}"
        metadata_dict[analysis_key] = row

print(f"DEBUG: Found {len(metadata_dict)} entries in metadata.")

##############################################################################
# 3) Define scattering logic: by chromosome or interval
##############################################################################
# By default, chromosomes for GRCh38:
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

def get_scatter_units():
    """
    Return either the list of chromosomes or
    a list like ["0000-scattered", "0001-scattered", ...].
    """
    if SCATTER_MODE == "chromosome":
        return CHROMOSOMES
    elif SCATTER_MODE == "interval":
        return [f"{i:04d}-scattered" for i in range(SCATTER_COUNT)]
    else:
        # If none, just do a single region "all"
        return ["all"]

SCATTER_UNITS = get_scatter_units()

##############################################################################
# 4) Helper function for memory usage
##############################################################################
def get_mem_from_threads(wildcards, threads):
    """
    Basic formula: 8800 MB per thread
    """
    return threads * 8800

##############################################################################
# 5) Prepare scatter intervals rules (only if SCATTER_MODE == 'interval')
##############################################################################
rule scatter_intervals_by_ns:
    """
    Scatter intervals across the reference genome’s Ns using GATK ScatterIntervalsByNs.
    Produces one big interval_list if we're doing 'interval' mode.
    Otherwise creates an empty file to satisfy input.
    """
    input:
        ref=REFERENCE_FILE
    output:
        intervals_list=os.path.join(INTERVALS_DIR, "scattered.interval_list")
    log:
        os.path.join(LOG_DIR, "scatter_intervals_by_ns.log")
    conda:
        "gatk"
    shell:
        r"""
        if [ "{SCATTER_MODE}" = "interval" ]; then
            echo "DEBUG: Starting ScatterIntervalsByNs" > {log}
            gatk ScatterIntervalsByNs \
                -R {input.ref} \
                -O {output.intervals_list} \
                &>> {log}
            echo "DEBUG: Finished ScatterIntervalsByNs" >> {log}
        else
            # If not interval mode, produce empty file
            touch {output.intervals_list}
        fi
        """

rule split_intervals:
    """
    Split intervals into SCATTER_COUNT parts with GATK SplitIntervals.
    Then produce multiple .interval_list in INTERVALS_DIR named e.g. 0000-scattered.interval_list.
    Mark them as temp so we can clean them up after merging.
    Also produce a .bed for each scattered interval_list (if needed).
    """
    input:
        intervals_list=rules.scatter_intervals_by_ns.output.intervals_list
    output:
        interval_lists=temp(expand(
            os.path.join(INTERVALS_DIR, "{i:04d}-scattered.interval_list"),
            i=range(SCATTER_COUNT)
        ))
    params:
        scatter_count=SCATTER_COUNT
    log:
        os.path.join(LOG_DIR, "split_intervals.log")
    conda:
        "gatk"
    shell:
        r"""
        if [ "{SCATTER_MODE}" = "interval" ]; then
            echo "DEBUG: Starting SplitIntervals" > {log}
            gatk SplitIntervals \
                -R {input.intervals_list} \
                -L {input.intervals_list} \
                --scatter-count {params.scatter_count} \
                -O {INTERVALS_DIR} &>> {log}
            echo "DEBUG: Finished SplitIntervals" >> {log}
        else
            # Not scattering by intervals => just create placeholders
            for i in $(seq 0 $(( {params.scatter_count} - 1 ))); do
                touch {INTERVALS_DIR}/$(printf "%04d-scattered.interval_list" $i)
            done
        fi
        """

##############################################################################
# 6) Write a per-analysis-key BAM list (if you want multi-bam input)
##############################################################################
rule write_bam_list:
    """
    For each analysis_key, create a small .txt file listing the tumor & normal BAMs
    (or just tumor if no normal is present).
    """
    output:
        bam_list=lambda wc: os.path.join(VARIANT_DIR, f"{wc.analysis_key}.bam_list.txt")
    run:
        # We fetch from metadata_dict, writing lines to the .txt
        analysis_info = metadata_dict[wildcards.analysis_key]
        bam1 = os.path.join(FINAL_BAM_FOLDER, analysis_info["bam1_file_basename"] + FINAL_BAM_EXTENSION)
        bam2 = ""
        if analysis_info.get("bam2_file_basename"):
            bam2 = os.path.join(FINAL_BAM_FOLDER, analysis_info["bam2_file_basename"] + FINAL_BAM_EXTENSION)
        with open(output.bam_list, "w") as f:
            f.write(bam1 + "\n")
            if bam2:
                f.write(bam2 + "\n")

##############################################################################
# 7) Call variants (Mutect2) for each scatter chunk
##############################################################################
rule call_variants_scatter:
    """
    Run Mutect2 per analysis_key and per scatter unit (chromosome or interval).
    Mark the scattered output VCF as temp so it can be removed post-merge.
    """
    input:
        bam_list=rules.write_bam_list.output,
        interval_list=lambda wc: os.path.join(INTERVALS_DIR, f"{wc.scatter_id}.interval_list"),
        ref=REFERENCE_FILE
    output:
        vcf=temp(f"{VARIANT_DIR}/{{analysis_key}}_{{scatter_id}}.vcf.gz")
    threads: 2
    resources:
        mem_mb=get_mem_from_threads,
        time="72:00:00",
        tmpdir=SCRATCH_DIR
    params:
        panel_of_normals=PANEL_OF_NORMALS,
        af_only_gnomad=AF_ONLY_GNOMAD
    log:
        mutect2_log=os.path.join(LOG_DIR, "{analysis_key}.{scatter_id}.mutect2.log")
    conda:
        "gatk"
    shell:
        r"""
        echo "DEBUG: Starting Mutect2 for {wildcards.analysis_key}, scatter {wildcards.scatter_id}" >&2

        # Build up multiple -I arguments from the bam_list
        INPUT_ARGS=""
        while read -r B; do
            INPUT_ARGS="$INPUT_ARGS -I $B"
        done < {input.bam_list}

        # Possibly use -L "chrX" or the interval_list, depending on SCATTER_MODE
        SCATTER_ARG=""
        if [ "{SCATTER_MODE}" = "chromosome" ]; then
            # scatter_id is e.g. 'chr1', 'chr2', etc. => -L {scatter_id}
            SCATTER_ARG="-L {wildcards.scatter_id}"
        elif [ "{SCATTER_MODE}" = "interval" ]; then
            SCATTER_ARG="-L {input.interval_list}"
        else
            # If 'none' => do entire genome
            SCATTER_ARG=""
        fi

        gatk --java-options '-Xms4000m -Xmx10g -Djava.io.tmpdir={resources.tmpdir}' Mutect2 \
            -R {input.ref} \
            $INPUT_ARGS \
            --panel-of-normals {params.panel_of_normals} \
            --germline-resource {params.af_only_gnomad} \
            --genotype-germline-sites true \
            --genotype-pon-sites true \
            $SCATTER_ARG \
            -O {output.vcf} \
            2> {log.mutect2_log}

        echo "DEBUG: Finished Mutect2 for {wildcards.analysis_key}, scatter {wildcards.scatter_id}" >&2
        """

##############################################################################
# 8) Merge scattered VCFs with GatherVcfs
##############################################################################
rule gather_vcfs:
    """
    Gather scattered VCFs into one final VCF for each analysis_key.
    """
    input:
        lambda wc: [
            f"{VARIANT_DIR}/{wc.analysis_key}_{scatter_unit}.vcf.gz"
            for scatter_unit in SCATTER_UNITS
        ]
    output:
        vcf=f"{VARIANT_DIR}/{{analysis_key}}.vcf.gz",
        tbi=f"{VARIANT_DIR}/{{analysis_key}}.vcf.gz.tbi"
    threads: 2
    resources:
        mem_mb=get_mem_from_threads,
        time="72:00:00",
        tmpdir=SCRATCH_DIR
    log:
        gather_log=os.path.join(LOG_DIR, "{analysis_key}.gather_vcfs.log")
    conda:
        "gatk"
    shell:
        r"""
        echo "DEBUG: Gathering VCFs for {wildcards.analysis_key}" >&2
        GATHER_INPUTS=""
        for chunk_vcf in {input}; do
            GATHER_INPUTS="$GATHER_INPUTS -I $chunk_vcf"
        done

        gatk --java-options '-Xms4000m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={resources.tmpdir}' GatherVcfs \
            $GATHER_INPUTS \
            -O {output.vcf} 2> {log.gather_log}

        tabix -p vcf {output.vcf}
        echo "DEBUG: Finished gather for {wildcards.analysis_key}" >&2
        """

##############################################################################
# 9) Final rule all
##############################################################################
rule all:
    """
    Produce final scattered-and-gathered VCFs for each analysis_key.
    (Does not include filtering, contamination steps, etc. here;
     those can be separate rules or in separate pipeline sections.)
    """
    input:
        # Ensure intervals are scattered if needed
        rules.split_intervals.output.interval_lists,
        # Then expand per analysis_key * scatter chunk
        expand(
            f"{VARIANT_DIR}/{{analysis_key}}_{{scatter_id}}.vcf.gz",
            analysis_key=metadata_dict.keys(),
            scatter_id=SCATTER_UNITS
        ),
        # Finally ensure the merged VCF
        expand(
            [
                f"{VARIANT_DIR}/{{analysis_key}}.vcf.gz",
                f"{VARIANT_DIR}/{{analysis_key}}.vcf.gz.tbi"
            ],
            analysis_key=metadata_dict.keys()
        )
