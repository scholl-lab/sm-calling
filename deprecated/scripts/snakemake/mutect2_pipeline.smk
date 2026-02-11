##############################################################################
# mutect2_pipeline.smk
##############################################################################
import os
import csv

##############################################################################
# 1) Load the config file
##############################################################################
configfile: "config_mutect2.yaml"

# Extract user-defined paths and references from config
FINAL_BAM_FOLDER         = config["final_bam_folder"]
OUTPUT_FOLDER            = config["output_folder"]
REFERENCE_UNPACKED       = config["reference_unpacked"]
PANEL_OF_NORMALS         = config["panel_of_normals"]
AF_ONLY_GNOMAD           = config["af_only_gnomad"]
COMMON_BIALLELIC_GNOMAD  = config["common_biallelic_gnomad"]
FINAL_BAM_EXTENSION      = config.get("final_bam_file_extension", ".bam")
SCATTER_MODE             = config.get("scatter_mode", "none")  # "chromosome", "interval", or "none"
SCATTER_COUNT            = config.get("scatter_count", 400)
INTERVALS_DIR            = config.get("intervals_dir", "analysis/intervals")
LOG_SUBFOLDER            = config["log_dir_sub"]
METADATA_FILE            = config["metadata_file"]

# Define chromosomes if scattering by chromosome
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

##############################################################################
# 2) Optional environment-based scratch directory (e.g., for cluster jobs)
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
        # e.g., 'analysis_key' = individual_analysis
        analysis_key = f"{row['individual1']}_{row['analysis']}"
        metadata_dict[analysis_key] = row

print(f"DEBUG: Number of entries in metadata: {len(metadata_dict)}")

##############################################################################
# 4) Define helper functions for resource usage and expansions
##############################################################################
def get_mem_from_threads(wildcards, threads):
    """Default: 4400 MB per thread."""
    return threads * 4400

def get_scatter_units():
    """
    Return the list of units by which we scatter.
    If 'chromosome', return CHROMOSOMES.
    If 'interval', return e.g. ['0000-scattered', '0001-scattered', ... ].
    If 'none', return ['all'].
    """
    if SCATTER_MODE == "chromosome":
        return CHROMOSOMES
    elif SCATTER_MODE == "interval":
        return [f"{i:04d}-scattered" for i in range(SCATTER_COUNT)]
    else:
        return ["all"]

SCATTER_UNITS = get_scatter_units()

##############################################################################
# 5) Make sure directories exist
##############################################################################
VARIANT_DIR          = os.path.join(OUTPUT_FOLDER, "variant_calls")
MERGED_DIR           = os.path.join(OUTPUT_FOLDER, "variant_merge")
CONTAM_DIR           = os.path.join(OUTPUT_FOLDER, "calculate_contamination")
FILTERED_DIR         = os.path.join(OUTPUT_FOLDER, "filtered_vcfs")
LOG_DIR              = os.path.join(OUTPUT_FOLDER, LOG_SUBFOLDER)
INTERVALS_OUTPUT_DIR = INTERVALS_DIR

for d in [VARIANT_DIR, MERGED_DIR, CONTAM_DIR, FILTERED_DIR, LOG_DIR, INTERVALS_OUTPUT_DIR]:
    os.makedirs(d, exist_ok=True)

##############################################################################
# 6) Final outputs & "rule all"
##############################################################################
rule all:
    """
    Produce final Filtered Mutect2 VCFs for each (individual, analysis) pair.
    """
    input:
        # a) Ensure scattered calls exist
        expand(
            f"{VARIANT_DIR}/{{analysis_key}}_{{chromosome}}.vcf.gz",
            analysis_key=metadata_dict.keys(),
            chromosome=SCATTER_UNITS
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
# 7) Scatter intervals if SCATTER_MODE is 'interval'
##############################################################################
if SCATTER_MODE == "interval":
    rule scatter_intervals_by_ns:
        """
        Scatter intervals across the reference genome's Ns using GATK ScatterIntervalsByNs.
        """
        input:
            reference=REFERENCE_UNPACKED
        output:
            scattered_interval_list=os.path.join(INTERVALS_OUTPUT_DIR, "scattered.interval_list")
        log:
            os.path.join(LOG_DIR, "scatter_intervals_by_ns.log")
        conda:
            "gatk"
        shell:
            r"""
            echo "DEBUG: Starting ScatterIntervalsByNs" > {log}
            gatk ScatterIntervalsByNs \
                -R {input.reference} \
                -O {output.scattered_interval_list} &>> {log}
            echo "DEBUG: Finished ScatterIntervalsByNs" >> {log}
            """

    rule split_intervals:
        """
        Split intervals into SCATTER_COUNT parts using GATK SplitIntervals,
        then convert each .interval_list => .bed
        """
        input:
            scattered_interval_list=os.path.join(INTERVALS_OUTPUT_DIR, "scattered.interval_list"),
            reference=REFERENCE_UNPACKED
        output:
            interval_lists=temp(expand(os.path.join(INTERVALS_OUTPUT_DIR, "{scatter_id}-scattered.interval_list"),
                                       scatter_id=[f"{i:04d}" for i in range(SCATTER_COUNT)])),
            bed_files=temp(expand(os.path.join(INTERVALS_OUTPUT_DIR, "{scatter_id}-scattered.bed"),
                                  scatter_id=[f"{i:04d}" for i in range(SCATTER_COUNT)]))
        log:
            os.path.join(LOG_DIR, "split_intervals.log")
        conda:
            "gatk"
        shell:
            r"""
            echo "DEBUG: Starting SplitIntervals" > {log}
            gatk SplitIntervals \
                -R {input.reference} \
                -L {input.scattered_interval_list} \
                --scatter-count {SCATTER_COUNT} \
                -O {INTERVALS_OUTPUT_DIR} &>> {log}
            echo "DEBUG: Finished SplitIntervals" >> {log}

            # Convert each .interval_list => .bed
            for interval in {INTERVALS_OUTPUT_DIR}/*-scattered.interval_list; do
                bed="${{interval%.interval_list}}.bed"
                sed 's/[:\-]/\t/g' "$interval" \
                    | grep -v "^@" \
                    | cut -f1-3 \
                    > "$bed"
            done
            """

##############################################################################
# 8) Call Variants (Mutect2) scattered by chromosome or interval
##############################################################################
rule call_variants:
    """
    Call variants with Mutect2, either scattered by chromosome or interval.
    """
    input:
        bam1=lambda wc: os.path.join(
            FINAL_BAM_FOLDER,
            metadata_dict[wc.analysis_key]["bam1_file_basename"] + FINAL_BAM_EXTENSION
        ),
        # Return a 1-element list if bam2_file_basename exists, otherwise an empty list
        bam2=lambda wc: (
            [os.path.join(
                FINAL_BAM_FOLDER,
                metadata_dict[wc.analysis_key]["bam2_file_basename"] + FINAL_BAM_EXTENSION
            )]
            if metadata_dict[wc.analysis_key].get("bam2_file_basename")
            else []
        ),
        # Return a 1-element list if scatter_mode=interval and chromosome!=all, otherwise empty list
        intervals=lambda wc: (
            [os.path.join(INTERVALS_OUTPUT_DIR, f"{wc.chromosome}-scattered.interval_list")]
            if SCATTER_MODE == "interval" and wc.chromosome != "all"
            else []
        )
    output:
        variant_file = f"{VARIANT_DIR}/{{analysis_key}}_{{chromosome}}.vcf.gz",
        stats = f"{VARIANT_DIR}/{{analysis_key}}_{{chromosome}}.vcf.gz.stats",
        f1r2 = f"{VARIANT_DIR}/{{analysis_key}}_{{chromosome}}.f1r2.tar.gz"
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
        
        # We now expect input.bam2 to be a list (length 0 or 1)
        # We can build an array of -I options from both bam1 and bam2
        BAMS=()
        
        # Always add bam1 as -I
        BAMS+=("-I {input.bam1}")
        
        # If bam2 is present in the list, add it
        for b2 in {input.bam2}; do
            if [ -s "$b2" ]; then
                echo "DEBUG: Found second BAM: $b2" >&2
                BAMS+=("-I $b2")
                # Also pass the normal sample name if desired
                if [ -n "{params.normal_sample}" ]; then
                    BAMS+=("-normal {params.normal_sample}")
                fi
            fi
        done
        
        # Intervals (list of length 0 or 1)
        SCATTER_ARGS=""
        for intervals_file in {input.intervals}; do
            if [ -s "$intervals_file" ]; then
                echo "DEBUG: Found intervals file: $intervals_file" >&2
                SCATTER_ARGS="-L $intervals_file"
            fi
        done
        
        if [ "{SCATTER_MODE}" = "chromosome" ] && [ "{wildcards.chromosome}" != "all" ]; then
            echo "DEBUG: Using chromosome scattering for {wildcards.chromosome}" >&2
            SCATTER_ARGS="-L {wildcards.chromosome}"
        fi
        
        # Construct the final BAMS string
        BAMS_STR="${{BAMS[*]}}"
        echo "DEBUG: Final BAMS_STR: $BAMS_STR" >&2
        
        gatk --java-options '-Xms4000m -Xmx10g -Djava.io.tmpdir={resources.tmpdir}' Mutect2 \
            -R "{params.reference}" \
            $BAMS_STR \
            --germline-resource "{params.af_gnomad}" \
            --panel-of-normals "{params.pon}" \
            --genotype-germline-sites true \
            --genotype-pon-sites true \
            --f1r2-tar-gz "{VARIANT_DIR}/{params.individual}_{params.analysis}_{wildcards.chromosome}.f1r2.tar.gz" \
            $SCATTER_ARGS \
            -O "{output.variant_file}" 2> "{log.mutect2_log}"

        echo "DEBUG: Finished Mutect2 call for {wildcards.analysis_key} {wildcards.chromosome}" >&2
        """

##############################################################################
# 9) Merge scattered calls (VCF, stats, f1r2)
##############################################################################
rule merge_vcfs:
    """
    Merge scattered VCF files across all chosen chromosomes/intervals.
    """
    input:
        vcf_files=lambda wc: [
            f"{VARIANT_DIR}/{wc.analysis_key}_{chrom}.vcf.gz"
            for chrom in SCATTER_UNITS
            if chrom != "all"
        ]
    output:
        vcf  = f"{MERGED_DIR}/{{analysis_key}}.vcf.gz",
        tbi  = f"{MERGED_DIR}/{{analysis_key}}.vcf.gz.tbi"
    params:
        input_files = lambda wc: ' '.join(
            f"-I {VARIANT_DIR}/{wc.analysis_key}_{chrom}.vcf.gz"
            for chrom in SCATTER_UNITS if chrom != "all"
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
        stats_files=lambda wc: [
            f"{VARIANT_DIR}/{wc.analysis_key}_{chrom}.vcf.gz.stats"
            for chrom in SCATTER_UNITS if chrom != "all"
        ]
    output:
        stats_merged = f"{MERGED_DIR}/{{analysis_key}}.vcf.gz.stats"
    params:
        input_files = lambda wc: ' '.join(
            f"-stats {VARIANT_DIR}/{wc.analysis_key}_{chrom}.vcf.gz.stats"
            for chrom in SCATTER_UNITS if chrom != "all"
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
    Merge scattered f1r2 tar.gz files across scatter units.
    """
    input:
        f1r2_files=lambda wc: [
            f"{VARIANT_DIR}/{wc.analysis_key}_{chrom}.f1r2.tar.gz"
            for chrom in SCATTER_UNITS if chrom != "all"
        ]
    output:
        tar_merged = f"{MERGED_DIR}/{{analysis_key}}_read-orientation-model.tar.gz"
    params:
        input_files = lambda wc: ' '.join(
            f"-I {VARIANT_DIR}/{wc.analysis_key}_{chrom}.f1r2.tar.gz"
            for chrom in SCATTER_UNITS if chrom != "all"
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
# 10) GetPileupSummaries for contamination (for each unique BAM)
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
# 11) CalculateContamination for each analysis_key
##############################################################################
rule calculate_contamination:
    """
    Calculate contamination for each tumor sample (and matched normal, if present).
    """
    input:
        # Return a 1-element list if bam1_file_basename is defined, otherwise empty
        tumor_pileup=lambda wc: (
            [
                os.path.join(
                    CONTAM_DIR,
                    metadata_dict[wc.analysis_key]["bam1_file_basename"] + ".getpileupsummaries.table"
                )
            ]
            if metadata_dict[wc.analysis_key].get("bam1_file_basename")
            else []
        ),
        # Return a 1-element list if bam2_file_basename is defined, otherwise empty
        normal_pileup=lambda wc: (
            [
                os.path.join(
                    CONTAM_DIR,
                    metadata_dict[wc.analysis_key]["bam2_file_basename"] + ".getpileupsummaries.table"
                )
            ]
            if metadata_dict[wc.analysis_key].get("bam2_file_basename")
            else []
        )
    output:
        contamination_table=f"{CONTAM_DIR}/{{analysis_key}}.contamination.table",
        segments_table=f"{CONTAM_DIR}/{{analysis_key}}.segments.table"
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

        # The lists {input.tumor_pileup} and {input.normal_pileup} might have 0 or 1 file each.
        # We'll extract them into shell variables if they exist:
        TUMOR_FILE=""
        for tfile in {input.tumor_pileup}; do
            if [ -s "$tfile" ]; then
                TUMOR_FILE="$tfile"
                echo "DEBUG: Found tumor pileup: $TUMOR_FILE" >&2
            fi
        done

        MATCHED_OPTION=""
        for nfile in {input.normal_pileup}; do
            if [ -s "$nfile" ]; then
                MATCHED_OPTION="--matched $nfile"
                echo "DEBUG: Found normal pileup: $nfile" >&2
            fi
        done

        if [ -z "$TUMOR_FILE" ]; then
            echo "WARNING: No tumor pileup table found for {wildcards.analysis_key}. Skipping." >&2
            # You could either skip or fail here.
            # If skipping, you might just 'touch' the output or raise an error:
            touch {output.contamination_table}
            touch {output.segments_table}
            exit 0
        fi

        echo "DEBUG: matched_option=$MATCHED_OPTION" >&2

        # Run GATK CalculateContamination
        gatk --java-options '-Xms4000m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={resources.tmpdir}' CalculateContamination \
            -I "$TUMOR_FILE" \
            $MATCHED_OPTION \
            --tumor-segmentation "{output.segments_table}" \
            -O "{output.contamination_table}" 2> "{log.calccont_log}"

        echo "DEBUG: Finished contamination for {wildcards.analysis_key}" >&2
        """

##############################################################################
# 12) Filter Mutect Calls
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
