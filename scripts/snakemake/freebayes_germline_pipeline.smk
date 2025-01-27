##############################################################################
# freebayes_pipeline.smk
##############################################################################
import os

# 1) Load the config file
configfile: "config.yaml"

# Retrieve config values
BAM_LIST_FILE     = config["bam_list_file"]
REFERENCE_GENOME  = config["reference_genome"]
SCATTER_MODE      = config["scatter_mode"]         # "chromosome" or "interval"
CHROMS            = config.get("chromosomes", [])
SCATTER_COUNT     = config.get("scatter_count", 25)
FREEBAYES_PARAMS  = config.get("freebayes_params", {})
RESULTS_DIR       = config["results_dir"]
LOGS_DIR          = config["logs_dir"]
INTERVALS_DIR     = config["intervals_dir"]
FINAL_MERGED_VCF  = config["final_merged_vcf"]

FREEBAYES_ENV     = config["freebayes_env"]  # e.g. "freebayes"
GATK_ENV          = config["gatk_env"]       # e.g. "gatk"

# Make sure directories exist
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(LOGS_DIR, exist_ok=True)
os.makedirs(INTERVALS_DIR, exist_ok=True)

##############################################################################
# 2) Read the BAM list file so Snakemake knows the inputs
##############################################################################
with open(BAM_LIST_FILE, "r") as f:
    BAM_FILES = [line.strip() for line in f if line.strip()]

print(f"DEBUG: Found {len(BAM_FILES)} BAM(s) in {BAM_LIST_FILE}")

##############################################################################
# 3) Function to convert freebayes_params dict to a command-line string
##############################################################################
def freebayes_param_str(params_dict):
    """Construct a command string from the dictionary in config."""
    parts = []
    for k, v in params_dict.items():
        parts.append(f"{k} {v}")
    return " ".join(parts)

FREEBAYES_CMD_PARAMS = freebayes_param_str(FREEBAYES_PARAMS)

##############################################################################
# 4) Define naming scheme for GATK's scattered intervals
#
# By default, GATK SplitIntervals can produce "0000-scattered.interval_list",
# "0001-scattered.interval_list", etc. We'll expect exactly that pattern.
##############################################################################
def scattered_interval_name(i):
    """
    Return e.g. "0000-scattered.interval_list" for i=0,
    "0001-scattered.interval_list" for i=1, etc.
    """
    return f"{i:04d}-scattered.interval_list"

def scattered_bed_name(i):
    """
    Return the corresponding bed name, e.g. "0000-scattered.bed".
    """
    return f"{i:04d}-scattered.bed"


# Create Python lists of expected scattered intervals & BED files:
SCATTERED_INTERVAL_LISTS = [
    os.path.join(INTERVALS_DIR, scattered_interval_name(i))
    for i in range(SCATTER_COUNT)
]
SCATTERED_BED_FILES = [
    os.path.join(INTERVALS_DIR, scattered_bed_name(i))
    for i in range(SCATTER_COUNT)
]


##############################################################################
# 5) The final "all" rule - produce one merged VCF
##############################################################################
rule all:
    input:
        FINAL_MERGED_VCF


##############################################################################
# 6) Scatter intervals by Ns (if SCATTER_MODE == 'interval')
##############################################################################
rule scatter_intervals_by_ns:
    """
    Scatter intervals across the reference genome's Ns using GATK ScatterIntervalsByNs.
    If SCATTER_MODE != 'interval', just create an empty file.
    """
    input:
        reference=REFERENCE_GENOME
    output:
        intervals_list=os.path.join(INTERVALS_DIR, "scattered.interval_list")
    log:
        os.path.join(LOGS_DIR, "scatter_intervals_by_ns.log")
    conda:
        GATK_ENV
    shell:
        r"""
        if [ "{SCATTER_MODE}" = "interval" ]; then
            set -e
            echo "Starting ScatterIntervalsByNs" > {log}
            gatk ScatterIntervalsByNs \
                -R {input.reference} \
                -O {output.intervals_list} &>> {log}
            echo "Finished ScatterIntervalsByNs" >> {log}
        else
            # If not interval mode, produce an empty intervals file
            touch {output.intervals_list}
        fi
        """


##############################################################################
# 7) Split intervals => produce "0000-scattered.interval_list" etc. + BED
##############################################################################
rule split_intervals:
    """
    Split intervals into SCATTER_COUNT parts using GATK SplitIntervals,
    then convert each .interval_list => .bed

    GATK by default names them: 0000-scattered.interval_list, 0001-scattered.interval_list, etc.
    We rely on that pattern, so we also specify --output-prefix "" below.
    """
    input:
        intervals_list=os.path.join(INTERVALS_DIR, "scattered.interval_list"),
        reference=REFERENCE_GENOME
    output:
        interval_lists=SCATTERED_INTERVAL_LISTS,
        bed_files=SCATTERED_BED_FILES
    params:
        scatter_count=SCATTER_COUNT
    log:
        os.path.join(LOGS_DIR, "split_intervals.log")
    conda:
        GATK_ENV
    shell:
        r"""
        if [ "{SCATTER_MODE}" = "interval" ]; then
            set -e
            echo "Starting SplitIntervals" > {log}
            # Force GATK to produce exactly "0000-scattered.interval_list", etc.
            gatk SplitIntervals \
                -R {input.reference} \
                -L {input.intervals_list} \
                --scatter-count {params.scatter_count} \
                -O {INTERVALS_DIR} &>> {log}
            echo "Finished SplitIntervals" >> {log}

            # Convert each .interval_list => .bed
            for f in {INTERVALS_DIR}/*-scattered.interval_list; do
                bed="${{f%.interval_list}}.bed"
                sed 's/:\|-/\t/g' "$f" \
                    | grep -v "^@" \
                    | cut -f1-3 \
                    > "$bed"
            done

        else
            # Not scattering => create empty placeholders
            for i in $(seq 0 $(( {{params.scatter_count}} - 1 ))); do
                touch {INTERVALS_DIR}/$(printf "%04d-scattered.interval_list" $i)
                touch {INTERVALS_DIR}/$(printf "%04d-scattered.bed" $i)
            done
        fi
        """


##############################################################################
# 8) Determine scatter units for subsequent rules
#    e.g. "0000-scattered", "0001-scattered", ...
##############################################################################
def get_scatter_units():
    if SCATTER_MODE == "chromosome":
        return CHROMS
    elif SCATTER_MODE == "interval":
        # Each scattered file is "####-scattered.interval_list"
        # We'll use "####-scattered" as the scatter_id
        return [f"{i:04d}-scattered" for i in range(SCATTER_COUNT)]
    else:
        return []

SCATTER_UNITS = get_scatter_units()


##############################################################################
# 9) FreeBayes scatter rule
##############################################################################
rule freebayes_scatter:
    input:
        reference=REFERENCE_GENOME,
        bam_files=BAM_FILES,
        bed_file=lambda wc: os.path.join(INTERVALS_DIR, f"{wc.scatter_id}.bed")
    output:
        vcf = f"{RESULTS_DIR}/freebayes_{{scatter_id}}.vcf",
        log = f"{LOGS_DIR}/freebayes_{{scatter_id}}.log"
    threads: 2
    resources:
        mem_mb = 8192
    conda:
        FREEBAYES_ENV
    shell:
        r"""
        set -e
        echo "Starting FreeBayes on {wildcards.scatter_id} (threads={threads}, mem={resources.mem_mb}MB)" >&2

        REGION_ARG=""
        if [ "{SCATTER_MODE}" = "chromosome" ]; then
            REGION_ARG="-r chr{wildcards.scatter_id}"
        elif [ "{SCATTER_MODE}" = "interval" ]; then
            REGION_ARG="--targets {input.bed_file}"
        fi

        BAM_LIST_TXT="{RESULTS_DIR}/bam_list_{wildcards.scatter_id}.txt"
        rm -f "$BAM_LIST_TXT"
        for b in {input.bam_files}; do
            echo "$b" >> "$BAM_LIST_TXT"
        done

        freebayes \
            -f {input.reference} \
            $REGION_ARG \
            -L "$BAM_LIST_TXT" \
            {FREEBAYES_CMD_PARAMS} \
            > {output.vcf} 2> {output.log}

        echo "Finished FreeBayes on {wildcards.scatter_id}" >&2
        """


##############################################################################
# 10) Merge scattered VCFs into one final VCF
##############################################################################
rule merge_vcfs:
    """
    Merge scattered FreeBayes VCFs into final merged VCF (bgzipped + indexed).
    """
    input:
        expand(
            os.path.join(RESULTS_DIR, "freebayes_{scatter_id}.vcf"),
            scatter_id=SCATTER_UNITS
        )
    output:
        merged_vcf=FINAL_MERGED_VCF
    log:
        os.path.join(LOGS_DIR, "merge_freebayes.log")
    shell:
        r"""
        set -e
        echo "Merging scattered FreeBayes VCFs" > {log}
        bcftools concat -O z -o {output.merged_vcf} {input} &>> {log}
        tabix -p vcf {output.merged_vcf}
        echo "Done merging final VCF." >> {log}
        """
