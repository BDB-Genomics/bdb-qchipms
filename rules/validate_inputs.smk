rule validate_inputs:
    input:
        samples=config["validate_inputs"]["samples_tsv"],
        protein_table=config["validate_inputs"]["protein_table"]
    output:
        touch("results/logs/inputs_validated.done")

    log: "results/logs/validate_inputs.log"
    benchmark: "results/benchmarks/validate_inputs.benchmark.txt"
    threads: config["validate_inputs"]["threads"]
    resources: mem_mb=get_mem_mb
    conda: "../envs/python_env.yaml"
    message: "Validating pipeline configuration and input matrices..."

    shell:
        """
        python scripts/validate_config.py config/config.yaml > {log} 2>&1
        """
