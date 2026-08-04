rule validate_inputs:
    input:
        samples=config["validate_inputs"]["samples_tsv"],
        protein_table=config["validate_inputs"]["protein_table"]
    output:
        touch("results/logs/inputs_validated.done")
    shell:
        """
        python scripts/validate_config.py config/config.yaml
        touch {output}
        """
