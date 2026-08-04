rule lfq_enrichment_and_plot:
    input:
        protein_table=config["lfq_enrichment"]["protein_table"],
        valid="results/logs/inputs_validated.done"
    output:
        enriched_tsv=config["lfq_enrichment"]["enriched_tsv"],
        volcano_png=config["lfq_enrichment"]["volcano_png"]

    log: "results/logs/lfq_enrichment.log"
    benchmark: "results/benchmarks/lfq_enrichment.benchmark.txt"
    threads: config["lfq_enrichment"]["threads"]
    resources: mem_mb=get_mem_mb
    conda: "../envs/r_env.yaml"
    message: "Performing LFQ statistical enrichment and generating volcano plot..."

    shell:
        """
        Rscript rules/scripts/qchip_ms_enrichment.R {input.protein_table} {output.enriched_tsv} {output.volcano_png} > {log} 2>&1
        """
