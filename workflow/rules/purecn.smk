"""PureCN copy number analysis rules.

Implements the PureCN pipeline for tumor purity, ploidy, and allele-specific
copy number estimation. Only active when purecn.enabled is true in config.
Requires Mutect2 caller (tumor_only or tumor_normal samples).
"""


if PURECN_ENABLED:

    rule purecn_intervals:
        """Generate optimised interval file for PureCN (one-time)."""
        localrule: True
        input:
            bed=PURECN_INTERVALS_BED,
            ref=REF,
        output:
            intervals=os.path.join(PURECN_REF_DIR, "intervals.txt"),
        params:
            genome=PURECN_GENOME,
        conda:
            "../envs/purecn.yaml"
        log:
            os.path.join(LOG_DIR, "purecn_intervals/intervals.log"),
        shell:
            r"""
            PURECN=$(Rscript -e 'cat(system.file("extdata", package = "PureCN"))')

            Rscript "$PURECN"/IntervalFile.R \
                --infile {input.bed} \
                --fasta {input.ref} \
                --outfile {output.intervals} \
                --offtarget \
                --genome {params.genome} \
                2> {log}
            """

    rule purecn_coverage:
        """Compute GC-normalised coverage for a BAM (per unique BAM)."""
        input:
            bam=lambda wc: os.path.join(BAM_FOLDER, wc.bam_basename + BAM_EXT),
            intervals=os.path.join(PURECN_REF_DIR, "intervals.txt"),
        output:
            coverage=os.path.join(
                PURECN_COV_DIR, "{bam_basename}_coverage_loess.txt.gz"
            ),
        params:
            outdir=PURECN_COV_DIR,
        retries: 2
        resources:
            mem_mb=lambda wildcards, attempt: 8000 * attempt,
            runtime=lambda wildcards, attempt: 720 * attempt,
        conda:
            "../envs/purecn.yaml"
        log:
            os.path.join(LOG_DIR, "purecn_coverage/{bam_basename}.log"),
        benchmark:
            os.path.join(LOG_DIR, "benchmarks/purecn_coverage/{bam_basename}.tsv")
        shell:
            r"""
            PURECN=$(Rscript -e 'cat(system.file("extdata", package = "PureCN"))')

            Rscript "$PURECN"/Coverage.R \
                --bam {input.bam} \
                --intervals {input.intervals} \
                --outdir {params.outdir} \
                2> {log}
            """

    rule purecn_normaldb:
        """Build normal database from normal BAM coverages (one-time)."""
        input:
            coverages=lambda wc: [
                os.path.join(PURECN_COV_DIR, f"{bn}_coverage_loess.txt.gz")
                for bn in get_normal_bams()
            ],
        output:
            normaldb=os.path.join(PURECN_REF_DIR, "normalDB.rds"),
            bias=os.path.join(PURECN_REF_DIR, "mapping_bias.rds"),
        params:
            outdir=PURECN_REF_DIR,
            genome=PURECN_GENOME,
        retries: 2
        resources:
            mem_mb=lambda wildcards, attempt: 16384 * attempt,
            runtime=lambda wildcards, attempt: 1440 * attempt,
        conda:
            "../envs/purecn.yaml"
        log:
            os.path.join(LOG_DIR, "purecn_normaldb/normaldb.log"),
        benchmark:
            os.path.join(LOG_DIR, "benchmarks/purecn_normaldb/normaldb.tsv")
        shell:
            r"""
            PURECN=$(Rscript -e 'cat(system.file("extdata", package = "PureCN"))')

            Rscript "$PURECN"/NormalDB.R \
                --coveragefiles {input.coverages} \
                --outdir {params.outdir} \
                --genome {params.genome} \
                2> {log}
            """

    def _get_purecn_normaldb(wildcards):
        """Return pre-built normalDB or pipeline-generated one."""
        if PURECN_NORMALDB:
            return PURECN_NORMALDB
        return os.path.join(PURECN_REF_DIR, "normalDB.rds")

    def _get_purecn_mapping_bias(wildcards):
        """Return pre-built mapping bias file or pipeline-generated one."""
        if PURECN_MAPPING_BIAS:
            return PURECN_MAPPING_BIAS
        return os.path.join(PURECN_REF_DIR, "mapping_bias.rds")

    def _get_purecn_snp_blacklist_arg(wildcards):
        """Return --snpblacklist arg or empty string."""
        if PURECN_SNP_BLACKLIST:
            return f"--snpblacklist {PURECN_SNP_BLACKLIST}"
        return ""

    rule purecn_run:
        """Run PureCN purity/ploidy estimation for a Mutect2 sample."""
        input:
            coverage=lambda wc: os.path.join(
                PURECN_COV_DIR,
                samples_df.loc[wc.sample, "tumor_bam"] + "_coverage_loess.txt.gz",
            ),
            vcf=os.path.join(MUTECT2_FILTERED_DIR, "{sample}.filtered.vcf.gz"),
            stats=os.path.join(MUTECT2_MERGED_DIR, "{sample}.vcf.gz.stats"),
            intervals=os.path.join(PURECN_REF_DIR, "intervals.txt"),
            normaldb=_get_purecn_normaldb,
            mapping_bias=_get_purecn_mapping_bias,
        output:
            rds=os.path.join(PURECN_RESULTS_DIR, "{sample}", "{sample}.rds"),
            csv=os.path.join(PURECN_RESULTS_DIR, "{sample}", "{sample}.csv"),
            loh=os.path.join(
                PURECN_RESULTS_DIR, "{sample}", "{sample}_loh.csv"
            ),
        params:
            outdir=lambda wc: os.path.join(PURECN_RESULTS_DIR, wc.sample),
            genome=PURECN_GENOME,
            seed=PURECN_SEED,
            postoptimize="--postoptimize" if PURECN_POSTOPTIMIZE else "",
            snp_blacklist=_get_purecn_snp_blacklist_arg,
            extra=PURECN_EXTRA,
        retries: 2
        resources:
            mem_mb=lambda wildcards, attempt: 16384 * attempt,
            runtime=lambda wildcards, attempt: 2880 * attempt,
        conda:
            "../envs/purecn.yaml"
        log:
            os.path.join(LOG_DIR, "purecn_run/{sample}.log"),
        benchmark:
            os.path.join(LOG_DIR, "benchmarks/purecn_run/{sample}.tsv")
        shell:
            r"""
            PURECN=$(Rscript -e 'cat(system.file("extdata", package = "PureCN"))')

            Rscript "$PURECN"/PureCN.R \
                --out {params.outdir}/{wildcards.sample} \
                --tumor {input.coverage} \
                --sampleid {wildcards.sample} \
                --vcf {input.vcf} \
                --statsfile {input.stats} \
                --normaldb {input.normaldb} \
                --mappingbiasfile {input.mapping_bias} \
                --intervals {input.intervals} \
                --genome {params.genome} \
                --force \
                --seed {params.seed} \
                {params.postoptimize} \
                {params.snp_blacklist} \
                {params.extra} \
                2> {log}
            """
