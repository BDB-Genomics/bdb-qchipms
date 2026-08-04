rule locus_annotation:
    input:
        enriched_tsv=config["locus_annotation"]["enriched_tsv"],
        peaks_bed=config["locus_annotation"]["peaks_bed"],
        chrom_sizes=config["locus_annotation"]["chrom_sizes"]
    output:
        annotated_bed=config["locus_annotation"]["annotated_bed"]

    log: "results/logs/locus_annotation.log"
    benchmark: "results/benchmarks/locus_annotation.benchmark.txt"
    threads: config["locus_annotation"]["threads"]
    resources: mem_mb=get_mem_mb
    conda: "../envs/python_env.yaml"
    message: "Annotating enriched proteins against MACS2 genomic peaks..."

    shell:
        """
        python scripts/intersect_peaks.py {input.enriched_tsv} {input.peaks_bed} {output.annotated_bed} > {log} 2>&1
        """
