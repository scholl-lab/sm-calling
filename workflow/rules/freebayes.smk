"""FreeBayes germline variant calling pipeline rules.

Implements scatter-parallel FreeBayes calling with merge and normalization.
Shares scatter rules from scatter.smk.
"""


def _get_germline_bam_paths():
    """Return full paths for all germline BAMs."""
    bams = []
    for _, row in samples_df[samples_df["analysis_type"] == "germline"].iterrows():
        bams.append(os.path.join(BAM_FOLDER, row["tumor_bam"] + BAM_EXT))
    return sorted(set(bams))


rule freebayes_call:
    """Scatter-parallel FreeBayes germline variant calling."""
    input:
        reference=REF,
        bam_files=_get_germline_bam_paths(),
    output:
        vcf=temp(os.path.join(FREEBAYES_SCATTER_DIR, "freebayes_{scatter_unit}.vcf.gz")),
    params:
        region_arg=get_freebayes_region_arg,
        extra=FREEBAYES_EXTRA,
    retries: 3
    resources:
        mem_mb=lambda wildcards, attempt: 16384 * (2 ** (attempt - 1)),
        runtime=lambda wildcards, attempt: 4320 * attempt,
    conda:
        "../envs/freebayes.yaml"
    log:
        os.path.join(LOG_DIR, "freebayes_call/{scatter_unit}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/freebayes_call/{scatter_unit}.tsv")
    shell:
        r"""
        # Write BAM list to temp file
        BAM_LIST=$(mktemp)
        trap 'rm -f "$BAM_LIST"' EXIT
        for b in {input.bam_files}; do
            echo "$b" >> "$BAM_LIST"
        done

        freebayes \
            -f {input.reference} \
            {params.region_arg} \
            -L "$BAM_LIST" \
            {params.extra} \
            2> {log} \
        | bgzip -c > {output.vcf}

        tabix -p vcf {output.vcf}
        """


rule merge_freebayes_vcfs:
    """Merge and normalize scattered FreeBayes VCFs into the final output."""
    input:
        vcfs=expand(
            os.path.join(FREEBAYES_SCATTER_DIR, "freebayes_{scatter_unit}.vcf.gz"),
            scatter_unit=SCATTER_UNITS,
        ),
    output:
        vcf=protected(os.path.join(FREEBAYES_DIR, "final_merged.vcf.gz")),
    params:
        norm_extra=BCFTOOLS_NORM_EXTRA,
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: 16384 * attempt,
        runtime=lambda wildcards, attempt: 1440 * attempt,
    conda:
        "../envs/bcftools.yaml"
    log:
        os.path.join(LOG_DIR, "merge_freebayes/merge.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/merge_freebayes/merge.tsv")
    shell:
        r"""
        bcftools concat \
            -n \
            -Oz \
            --threads 4 \
            {input.vcfs} 2> {log} \
        | bcftools norm \
            {params.norm_extra} \
            -f {REF} \
            --threads 4 \
            - 2>> {log} \
        | bcftools +fill-tags \
            -Oz \
            -o {output.vcf} \
            --threads 4 \
            -W=tbi - -- \
            2>> {log}
        """
