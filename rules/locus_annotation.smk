rule locus_annotation:
    input:
        enriched_tsv=config["locus_annotation"]["enriched_tsv"],
        peaks_bed=config["locus_annotation"]["peaks_bed"],
        chrom_sizes=config["locus_annotation"]["chrom_sizes"]
    output:
        annotated_bed=config["locus_annotation"]["annotated_bed"]
    resources:
        mem_mb=get_mem_mb,
        threads=config["resources"]["threads"]
    shell:
        """
        python scripts/intersect_peaks.py {input.enriched_tsv} {input.peaks_bed} {output.annotated_bed}
        """
