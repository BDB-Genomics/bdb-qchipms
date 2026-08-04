rule summary_report:
    input:
        enriched_tsv=config["summary_report"]["enriched_tsv"],
        annotated_bed=config["summary_report"]["annotated_bed"],
        volcano_png=config["lfq_enrichment"]["volcano_png"]
    output:
        html_report=config["summary_report"]["html_report"]
    resources:
        mem_mb=get_mem_mb,
        threads=1
    shell:
        """
        python scripts/summary_report.py {input.enriched_tsv} {input.annotated_bed} {output.html_report}
        """
