configfile: "config/config.yaml"

def get_mem_mb(wildcards, attempt):
    """Dynamically allocate memory based on attempt number"""
    return min(attempt * config["resources"]["base_mem_mb"], config["resources"]["max_mem_mb"])

include: "rules/validate_inputs.smk"
include: "rules/lfq_enrichment_and_plot.smk"
include: "rules/locus_annotation.smk"
include: "rules/summary_report.smk"

rule all:
    input:
        config["lfq_enrichment"]["enriched_tsv"],
        config["lfq_enrichment"]["volcano_png"],
        config["locus_annotation"]["annotated_bed"],
        config["summary_report"]["html_report"]

onstart:
    print("--------------------------------------------------")
    print("🚀 [qChIP-MS Pipeline] Execution Started")
    print("--------------------------------------------------")

onsuccess:
    print("--------------------------------------------------")
    print("✅ [qChIP-MS Pipeline] Completed Successfully!")
    print("--------------------------------------------------")

onerror:
    print("--------------------------------------------------")
    print("❌ [qChIP-MS Pipeline] Execution Failed!")
    print("Please check rule log files in results/logs/ for details.")
    print("--------------------------------------------------")
