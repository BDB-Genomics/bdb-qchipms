rule lfq_enrichment_and_plot:
    input:
        protein_table=config["lfq_enrichment"]["protein_table"],
        valid="results/logs/inputs_validated.done"
    output:
        enriched_tsv=config["lfq_enrichment"]["enriched_tsv"],
        volcano_png=config["lfq_enrichment"]["volcano_png"]
    resources:
        mem_mb=get_mem_mb,
        threads=config["resources"]["threads"]
    shell:
        """
        Rscript rules/scripts/qchip_ms_enrichment.R {input.protein_table} {output.enriched_tsv} {output.volcano_png}
        """
