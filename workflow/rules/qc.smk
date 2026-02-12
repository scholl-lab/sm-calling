"""Variant calling quality control rules.

Generates per-sample bcftools stats and an aggregated MultiQC report.
"""


rule bcftools_stats_mutect2:
    """Per-sample variant statistics for Mutect2 filtered VCFs."""
    input:
        vcf=os.path.join(MUTECT2_FILTERED_DIR, "{sample}.filtered.vcf.gz"),
        ref=REF,
    output:
        stats=os.path.join(QC_DIR, "bcftools_stats", "mutect2", "{sample}.stats.txt"),
    params:
        extra=BCFTOOLS_STATS_EXTRA,
    conda:
        "../envs/bcftools.yaml"
    log:
        os.path.join(LOG_DIR, "bcftools_stats/mutect2/{sample}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/bcftools_stats/mutect2/{sample}.tsv")
    shell:
        r"""
        bcftools stats \
            -s- \
            -F {input.ref} \
            {params.extra} \
            {input.vcf} \
            > {output.stats} \
            2> {log}
        """


rule bcftools_stats_freebayes:
    """Variant statistics for FreeBayes merged VCF."""
    input:
        vcf=os.path.join(FREEBAYES_DIR, "final_merged.vcf.gz"),
        ref=REF,
    output:
        stats=os.path.join(
            QC_DIR, "bcftools_stats", "freebayes", "all_samples.stats.txt"
        ),
    params:
        extra=BCFTOOLS_STATS_EXTRA,
    conda:
        "../envs/bcftools.yaml"
    log:
        os.path.join(LOG_DIR, "bcftools_stats/freebayes/all_samples.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/bcftools_stats/freebayes/all_samples.tsv")
    shell:
        r"""
        bcftools stats \
            -s- \
            -F {input.ref} \
            {params.extra} \
            {input.vcf} \
            > {output.stats} \
            2> {log}
        """


rule multiqc:
    """Aggregate all QC stats into a single MultiQC HTML report."""
    input:
        stats=get_qc_inputs(),
    output:
        html=os.path.join(QC_DIR, "multiqc_report.html"),
        data=directory(os.path.join(QC_DIR, "multiqc_data")),
    params:
        use_input_files_only=True,
    conda:
        "../envs/multiqc.yaml"
    log:
        os.path.join(LOG_DIR, "multiqc/multiqc.log"),
    shell:
        r"""
        multiqc \
            --force \
            --no-data-dir \
            --outdir $(dirname {output.html}) \
            --filename $(basename {output.html}) \
            {input.stats} \
            2> {log}
        """
