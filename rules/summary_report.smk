rule summary_report:
    input:
        enriched_tsv=config["summary_report"]["enriched_tsv"],
        annotated_bed=config["summary_report"]["annotated_bed"],
        volcano_png=config["lfq_enrichment"]["volcano_png"]
    output:
        html_report=config["summary_report"]["html_report"]

    log: "results/logs/summary_report.log"
    benchmark: "results/benchmarks/summary_report.benchmark.txt"
    threads: config["summary_report"]["threads"]
    resources: mem_mb=get_mem_mb
    conda: "../envs/python_env.yaml"
    message: "Compiling final qChIP-MS HTML summary report..."

    shell:
        """
        python scripts/summary_report.py {input.enriched_tsv} {input.annotated_bed} {output.html_report} > {log} 2>&1
        """
