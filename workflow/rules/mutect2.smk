"""Mutect2 somatic variant calling pipeline rules.

Implements the full GATK Mutect2 best-practice workflow:
scatter -> call -> gather -> contamination -> filter.
"""


wildcard_constraints:
    sample="[A-Za-z0-9_.-]+",
    scatter_unit="(chr[0-9XYMT]+|[0-9]{4}-scattered|all)",
    bam_basename="[A-Za-z0-9_.-]+",


rule pon_mutect2:
    """Per-normal tumor-only Mutect2 calling for cohort Panel of Normals."""
    input:
        bam=lambda wc: os.path.join(BAM_FOLDER, wc.bam_basename + BAM_EXT),
        reference=REF,
        target_bed=[TARGET_BED] if TARGET_BED else [],
    output:
        vcf=os.path.join(MUTECT2_PON_DIR, "{bam_basename}.pon.vcf.gz"),
        tbi=os.path.join(MUTECT2_PON_DIR, "{bam_basename}.pon.vcf.gz.tbi"),
        stats=os.path.join(MUTECT2_PON_DIR, "{bam_basename}.pon.vcf.gz.stats"),
    params:
        java_opts=get_java_opts,
        interval_arg=f"-L {TARGET_BED}" if TARGET_BED else "",
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 17600 * attempt,
        runtime=lambda wildcards, attempt: 1440 * attempt,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "pon_mutect2/{bam_basename}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/pon_mutect2/{bam_basename}.tsv")
    shell:
        r"""
        gatk --java-options '{params.java_opts}' Mutect2 \
            -R {input.reference} \
            -I {input.bam} \
            --max-mnp-distance 0 \
            --native-pair-hmm-threads {threads} \
            {params.interval_arg} \
            -O {output.vcf} \
            2> {log}
        """


rule pon_genomicsdb:
    """Combine per-normal PoN VCFs into a GenomicsDB workspace."""
    input:
        vcfs=lambda wc: [
            os.path.join(MUTECT2_PON_DIR, f"{bn}.pon.vcf.gz")
            for bn in get_normal_bams()
        ],
        reference=REF,
        target_bed=[TARGET_BED] if TARGET_BED else [],
    output:
        db=temp(directory(os.path.join(MUTECT2_PON_DIR, "pon_db"))),
    params:
        java_opts=get_java_opts,
        vcf_args=lambda wc: " ".join(
            f"-V {os.path.join(MUTECT2_PON_DIR, f'{bn}.pon.vcf.gz')}"
            for bn in get_normal_bams()
        ),
        interval_arg=(
            f"-L {TARGET_BED}"
            if TARGET_BED
            else " ".join(f"-L {c}" for c in CHROMOSOMES)
        ),
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 17600 * attempt,
        runtime=lambda wildcards, attempt: 1440 * attempt,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "pon_genomicsdb/genomicsdb.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/pon_genomicsdb/genomicsdb.tsv")
    shell:
        r"""
        rm -rf {output.db}
        gatk --java-options '{params.java_opts}' GenomicsDBImport \
            -R {input.reference} \
            --genomicsdb-workspace-path {output.db} \
            --merge-input-intervals true \
            {params.interval_arg} \
            {params.vcf_args} \
            2> {log}
        """


rule pon_create:
    """Create cohort Somatic Panel of Normals from GenomicsDB."""
    input:
        db=os.path.join(MUTECT2_PON_DIR, "pon_db"),
        reference=REF,
        gnomad=AF_GNOMAD,
    output:
        vcf=os.path.join(MUTECT2_PON_DIR, "pon.vcf.gz"),
        tbi=os.path.join(MUTECT2_PON_DIR, "pon.vcf.gz.tbi"),
    params:
        java_opts=get_java_opts,
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 17600 * attempt,
        runtime=lambda wildcards, attempt: 1440 * attempt,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "pon_create/pon.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/pon_create/pon.tsv")
    shell:
        r"""
        gatk --java-options '{params.java_opts}' CreateSomaticPanelOfNormals \
            -R {input.reference} \
            --germline-resource {input.gnomad} \
            -V gendb://{input.db} \
            -O {output.vcf} \
            2> {log}
        """


rule pon_mask_hotspots:
    """Mask recurrent somatic hotspot regions out of the cohort Panel of Normals."""
    input:
        vcf=os.path.join(MUTECT2_PON_DIR, "pon.vcf.gz"),
        tbi=os.path.join(MUTECT2_PON_DIR, "pon.vcf.gz.tbi"),
        mask_bed=lambda wc: PON_MASK_BED,
    output:
        vcf=os.path.join(MUTECT2_PON_DIR, "pon.masked.vcf.gz"),
        tbi=os.path.join(MUTECT2_PON_DIR, "pon.masked.vcf.gz.tbi"),
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
        runtime=lambda wildcards, attempt: 120 * attempt,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "pon_mask_hotspots/pon_masked.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/pon_mask_hotspots/pon_masked.tsv")
    shell:
        r"""
        bcftools view \
            -T ^{input.mask_bed} \
            -Oz \
            -W=tbi \
            -o {output.vcf} \
            {input.vcf} \
            2> {log}
        """


rule mutect2_call:
    """Scatter-parallel Mutect2 somatic variant calling."""
    input:
        tumor_bam=get_tumor_bam,
        normal_bam=get_normal_bam,
        reference=REF,
        pon=get_pon,
        gnomad=AF_GNOMAD,
        intervals=get_interval_input,
    output:
        vcf=temp(
            ensure(
                os.path.join(MUTECT2_SCATTER_DIR, "{sample}_{scatter_unit}.vcf.gz"),
                non_empty=True,
            )
        ),
        stats=temp(os.path.join(MUTECT2_SCATTER_DIR, "{sample}_{scatter_unit}.vcf.gz.stats")),
        f1r2=temp(os.path.join(MUTECT2_SCATTER_DIR, "{sample}_{scatter_unit}.f1r2.tar.gz")),
    params:
        tumor_input_arg=get_tumor_input_arg,
        normal_input_arg=get_normal_input_arg,
        interval_arg=get_interval_arg,
        java_opts=get_java_opts,
        extra=MUTECT2_ALL_ARGS,
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 17600 * attempt,
        runtime=lambda wildcards, attempt: 1440 * attempt,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "mutect2_call/{sample}_{scatter_unit}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/mutect2_call/{sample}_{scatter_unit}.tsv")
    shell:
        r"""
        extract_sm() {{
            local sm
            sm=$(samtools view -H "$1" | awk -F'\t' '/^@RG/ {{s=""; for (i=1; i<=NF; i++) if ($i ~ /^SM:/) s=substr($i, 4); print s}}' | sort -u)
            if [ -z "$sm" ] || [ "$(printf '%s\n' "$sm" | wc -l)" -ne 1 ]; then
                echo "ERROR: $1 must contain a single non-empty @RG SM tag" >&2
                return 1
            fi
            printf '%s\n' "$sm"
        }}

        : > {log}
        tumor_sms=()
        for tbam in {input.tumor_bam}; do
            tsm=$(extract_sm "$tbam") 2>> {log}
            tumor_sms+=("$tsm")
        done
        normal_arg=()
        if [ -n "{input.normal_bam}" ]; then
            normal_sm=$(extract_sm "{input.normal_bam}") 2>> {log}
            for tsm in "${{tumor_sms[@]}}"; do
                if [ "$tsm" = "$normal_sm" ]; then
                    echo "ERROR: tumor_sm ($tsm) equals normal_sm ($normal_sm)" >> {log}
                    exit 1
                fi
            done
            normal_arg=(-normal "$normal_sm")
        fi

        gatk --java-options '{params.java_opts}' Mutect2 \
            -R {input.reference} \
            {params.tumor_input_arg} \
            {params.normal_input_arg} \
            "${{normal_arg[@]}}" \
            --germline-resource {input.gnomad} \
            --panel-of-normals {input.pon} \
            --f1r2-tar-gz {output.f1r2} \
            --native-pair-hmm-threads {threads} \
            {params.extra} \
            {params.interval_arg} \
            -O {output.vcf} \
            2>> {log}
        """


rule gather_mutect2_vcfs:
    """Merge scattered VCF files using GATK GatherVcfs."""
    input:
        vcfs=lambda wc: [
            os.path.join(MUTECT2_SCATTER_DIR, f"{wc.sample}_{unit}.vcf.gz")
            for unit in SCATTER_UNITS
        ],
    output:
        vcf=os.path.join(MUTECT2_MERGED_DIR, "{sample}.vcf.gz"),
        tbi=os.path.join(MUTECT2_MERGED_DIR, "{sample}.vcf.gz.tbi"),
    params:
        input_args=lambda wc: " ".join(
            f"-I {os.path.join(MUTECT2_SCATTER_DIR, f'{wc.sample}_{unit}.vcf.gz')}"
            for unit in SCATTER_UNITS
        ),
        java_opts=get_java_opts,
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 17600 * attempt,
        runtime=lambda wildcards, attempt: 1440 * attempt,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "gather_vcfs/{sample}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/gather_vcfs/{sample}.tsv")
    shell:
        r"""
        gatk --java-options '{params.java_opts}' GatherVcfs \
            {params.input_args} \
            -O {output.vcf} \
            2> {log}

        tabix -p vcf {output.vcf}
        """


rule merge_mutect2_stats:
    """Merge scattered .stats files using GATK MergeMutectStats."""
    input:
        stats=lambda wc: [
            os.path.join(
                MUTECT2_SCATTER_DIR,
                f"{wc.sample}_{unit}.vcf.gz.stats",
            )
            for unit in SCATTER_UNITS
        ],
    output:
        stats=os.path.join(MUTECT2_MERGED_DIR, "{sample}.vcf.gz.stats"),
    params:
        input_args=lambda wc: " ".join(
            f"-stats {os.path.join(MUTECT2_SCATTER_DIR, f'{wc.sample}_{unit}.vcf.gz.stats')}"
            for unit in SCATTER_UNITS
        ),
        java_opts=get_java_opts,
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
        runtime=lambda wildcards, attempt: 720 * attempt,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "merge_stats/{sample}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/merge_stats/{sample}.tsv")
    shell:
        r"""
        gatk --java-options '{params.java_opts}' MergeMutectStats \
            {params.input_args} \
            -O {output.stats} \
            2> {log}
        """


rule learn_read_orientation:
    """Merge scattered f1r2 tar.gz files using GATK LearnReadOrientationModel."""
    input:
        f1r2=lambda wc: [
            os.path.join(
                MUTECT2_SCATTER_DIR,
                f"{wc.sample}_{unit}.f1r2.tar.gz",
            )
            for unit in SCATTER_UNITS
        ],
    output:
        model=os.path.join(MUTECT2_MERGED_DIR, "{sample}_read-orientation-model.tar.gz"),
    params:
        input_args=lambda wc: " ".join(
            f"-I {os.path.join(MUTECT2_SCATTER_DIR, f'{wc.sample}_{unit}.f1r2.tar.gz')}"
            for unit in SCATTER_UNITS
        ),
        java_opts=get_java_opts,
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 17600 * attempt,
        runtime=lambda wildcards, attempt: 1440 * attempt,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "learn_read_orientation/{sample}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/learn_read_orientation/{sample}.tsv")
    shell:
        r"""
        gatk --java-options '{params.java_opts}' LearnReadOrientationModel \
            {params.input_args} \
            -O {output.model} \
            2> {log}
        """


rule get_pileup_summaries:
    """Generate pileup summaries for each unique BAM (for contamination estimation)."""
    input:
        bam=lambda wc: os.path.join(BAM_FOLDER, wc.bam_basename + BAM_EXT),
        reference=REF,
        gnomad=COMMON_GNOMAD,
        target_bed=[TARGET_BED] if TARGET_BED else [],
    output:
        table=os.path.join(MUTECT2_CONTAM_DIR, "{bam_basename}.getpileupsummaries.table"),
    params:
        java_opts=get_java_opts,
        interval_arg=(
            f"-L {TARGET_BED} --interval-set-rule INTERSECTION" if TARGET_BED else ""
        ),
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 17600 * attempt,
        runtime=lambda wildcards, attempt: 360 * attempt,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "get_pileup_summaries/{bam_basename}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/get_pileup_summaries/{bam_basename}.tsv")
    shell:
        r"""
        gatk --java-options '{params.java_opts}' GetPileupSummaries \
            -I {input.bam} \
            -R {input.reference} \
            -V {input.gnomad} \
            -L {input.gnomad} \
            {params.interval_arg} \
            -O {output.table} \
            2> {log}
        """


rule calculate_contamination:
    """Calculate contamination per (sample, tumor_bam) pair (with optional matched normal)."""
    input:
        tumor_pileup=get_tumor_pileup,
        normal_pileup=get_normal_pileup,
    output:
        contamination=os.path.join(
            MUTECT2_CONTAM_DIR, "{sample}_{bam_basename}.contamination.table"
        ),
        segments=os.path.join(
            MUTECT2_CONTAM_DIR, "{sample}_{bam_basename}.segments.table"
        ),
    params:
        java_opts=get_java_opts,
        matched_arg=lambda wc: (
            f"--matched {get_normal_pileup(wc)[0]}" if get_normal_pileup(wc) else ""
        ),
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 17600 * attempt,
        runtime=lambda wildcards, attempt: 1440 * attempt,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "calculate_contamination/{sample}_{bam_basename}.log"),
    benchmark:
        os.path.join(
            LOG_DIR, "benchmarks/calculate_contamination/{sample}_{bam_basename}.tsv"
        )
    shell:
        r"""
        gatk --java-options '{params.java_opts}' CalculateContamination \
            -I {input.tumor_pileup} \
            {params.matched_arg} \
            --tumor-segmentation {output.segments} \
            -O {output.contamination} \
            2> {log}
        """


rule filter_mutect_calls:
    """Apply Mutect2 filters using contamination, segmentation, and orientation data."""
    input:
        vcf=os.path.join(MUTECT2_MERGED_DIR, "{sample}.vcf.gz"),
        reference=REF,
        stats=os.path.join(MUTECT2_MERGED_DIR, "{sample}.vcf.gz.stats"),
        contamination=get_contamination_tables,
        segmentation=get_segmentation_tables,
        ob_priors=os.path.join(MUTECT2_MERGED_DIR, "{sample}_read-orientation-model.tar.gz"),
    output:
        vcf=ensure(
            protected(os.path.join(MUTECT2_FILTERED_DIR, "{sample}.filtered.vcf.gz")),
            non_empty=True,
        ),
        tbi=protected(os.path.join(MUTECT2_FILTERED_DIR, "{sample}.filtered.vcf.gz.tbi")),
    params:
        java_opts=get_java_opts,
        contam_args=get_contamination_args,
        extra=FILTER_MUTECT_CALLS_EXTRA,
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 17600 * attempt,
        runtime=lambda wildcards, attempt: 1440 * attempt,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "filter_mutect/{sample}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/filter_mutect/{sample}.tsv")
    shell:
        r"""
        tmp_vcf="{output.vcf}.raw.vcf.gz"
        trap 'rm -f "$tmp_vcf" "$tmp_vcf.tbi" "$tmp_vcf.filteringStats.tsv"' EXIT

        gatk --java-options '{params.java_opts}' FilterMutectCalls \
            -R {input.reference} \
            -V {input.vcf} \
            --stats {input.stats} \
            {params.contam_args} \
            --ob-priors {input.ob_priors} \
            {params.extra} \
            -O "$tmp_vcf" \
            2> {log}

        bcftools norm \
            -m-any \
            -f {input.reference} \
            -Oz \
            -W=tbi \
            -o {output.vcf} \
            "$tmp_vcf" \
            2>> {log}
        """
