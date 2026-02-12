"""Mutect2 somatic variant calling pipeline rules.

Implements the full GATK Mutect2 best-practice workflow:
scatter -> call -> gather -> contamination -> filter.
"""


wildcard_constraints:
    sample="[A-Za-z0-9_.-]+",
    scatter_unit="(chr[0-9XYMT]+|[0-9]{4}-scattered|all)",
    bam_basename="[A-Za-z0-9_.-]+",


rule mutect2_call:
    """Scatter-parallel Mutect2 somatic variant calling."""
    input:
        tumor_bam=get_tumor_bam,
        normal_bam=get_normal_bam,
        reference=REF,
        pon=PON,
        gnomad=AF_GNOMAD,
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
        normal_arg=get_normal_name,
        interval_arg=get_interval_arg,
        java_opts=get_java_opts,
        extra=MUTECT2_ALL_ARGS,
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 17600 * attempt,
        runtime=lambda wildcards, attempt: 4320 * attempt,
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(LOG_DIR, "mutect2_call/{sample}_{scatter_unit}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/mutect2_call/{sample}_{scatter_unit}.tsv")
    shell:
        r"""
        gatk --java-options '{params.java_opts}' Mutect2 \
            -R {input.reference} \
            -I {input.tumor_bam} \
            {params.normal_arg} \
            --germline-resource {input.gnomad} \
            --panel-of-normals {input.pon} \
            --f1r2-tar-gz {output.f1r2} \
            {params.extra} \
            {params.interval_arg} \
            -O {output.vcf} \
            2> {log}
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
    output:
        table=os.path.join(MUTECT2_CONTAM_DIR, "{bam_basename}.getpileupsummaries.table"),
    params:
        java_opts=get_java_opts,
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 17600 * attempt,
        runtime=lambda wildcards, attempt: 4320 * attempt,
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
            -O {output.table} \
            2> {log}
        """


rule calculate_contamination:
    """Calculate contamination for each sample (with optional matched normal)."""
    input:
        tumor_pileup=get_tumor_pileup,
        normal_pileup=get_normal_pileup,
    output:
        contamination=os.path.join(MUTECT2_CONTAM_DIR, "{sample}.contamination.table"),
        segments=os.path.join(MUTECT2_CONTAM_DIR, "{sample}.segments.table"),
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
        os.path.join(LOG_DIR, "calculate_contamination/{sample}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/calculate_contamination/{sample}.tsv")
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
        contamination=os.path.join(MUTECT2_CONTAM_DIR, "{sample}.contamination.table"),
        segmentation=os.path.join(MUTECT2_CONTAM_DIR, "{sample}.segments.table"),
        ob_priors=os.path.join(MUTECT2_MERGED_DIR, "{sample}_read-orientation-model.tar.gz"),
    output:
        vcf=ensure(
            protected(os.path.join(MUTECT2_FILTERED_DIR, "{sample}.filtered.vcf.gz")),
            non_empty=True,
        ),
    params:
        java_opts=get_java_opts,
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
        gatk --java-options '{params.java_opts}' FilterMutectCalls \
            -R {input.reference} \
            -V {input.vcf} \
            --stats {input.stats} \
            --contamination-table {input.contamination} \
            --tumor-segmentation {input.segmentation} \
            --ob-priors {input.ob_priors} \
            -O {output.vcf} \
            2> {log}
        """
