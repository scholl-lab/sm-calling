"""Interval scattering rules shared by all callers.

Only activated when scatter mode is "interval". Chromosome-based scattering
needs no preprocessing -- chromosomes are passed directly as -L/-r arguments.
"""

if SCATTER_MODE == "interval":

    rule scatter_intervals_by_ns:
        """Create an interval list from reference genome Ns using GATK ScatterIntervalsByNs."""
        input:
            reference=REF,
        output:
            interval_list=os.path.join(INTERVALS_DIR, "scattered.interval_list"),
        log:
            os.path.join(LOG_DIR, "scatter/scatter_intervals_by_ns.log"),
        benchmark:
            os.path.join(LOG_DIR, "benchmarks/scatter/scatter_intervals_by_ns.tsv")
        conda:
            "../envs/gatk.yaml"
        shell:
            r"""
            gatk ScatterIntervalsByNs \
                -R {input.reference} \
                -O {output.interval_list} \
                2> {log}
            """

    rule split_intervals:
        """Split the scattered interval list into SCATTER_COUNT parts."""
        input:
            interval_list=os.path.join(INTERVALS_DIR, "scattered.interval_list"),
            reference=REF,
        output:
            interval_lists=temp(
                expand(
                    os.path.join(INTERVALS_DIR, "{scatter_id}-scattered.interval_list"),
                    scatter_id=[f"{i:04d}" for i in range(SCATTER_COUNT)],
                )
            ),
            bed_files=temp(
                expand(
                    os.path.join(INTERVALS_DIR, "{scatter_id}-scattered.bed"),
                    scatter_id=[f"{i:04d}" for i in range(SCATTER_COUNT)],
                )
            ),
        params:
            scatter_count=SCATTER_COUNT,
            intervals_dir=INTERVALS_DIR,
        log:
            os.path.join(LOG_DIR, "scatter/split_intervals.log"),
        benchmark:
            os.path.join(LOG_DIR, "benchmarks/scatter/split_intervals.tsv")
        conda:
            "../envs/gatk.yaml"
        shell:
            r"""
            gatk SplitIntervals \
                -R {input.reference} \
                -L {input.interval_list} \
                --scatter-count {params.scatter_count} \
                -O {params.intervals_dir} \
                2> {log}

            # Convert each .interval_list to .bed
            for interval in {params.intervals_dir}/*-scattered.interval_list; do
                bed="${{interval%.interval_list}}.bed"
                sed 's/[:\-]/\t/g' "$interval" \
                    | grep -v "^@" \
                    | cut -f1-3 \
                    > "$bed"
            done
            """
