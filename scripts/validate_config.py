#!/usr/bin/env python3
"""Validate the qChIP-MS config and input files before Snakemake builds the DAG."""

import sys
from pathlib import Path
import yaml
import csv


class Colors:
    HEADER = "\033[95m"
    BLUE = "\033[94m"
    CYAN = "\033[96m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"


from typing import NoReturn


def fail(errors: list[str]) -> NoReturn:
    print(
        f"\n{Colors.BOLD}{Colors.RED}┏━ CONFIGURATION VALIDATION FAILED ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓{Colors.ENDC}"
    )
    for m in errors:
        print(f"    • {m}")
    print(
        f"\n{Colors.BOLD}{Colors.RED}┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛{Colors.ENDC}"
    )
    sys.exit(1)


def load_config(config_path: Path) -> dict:
    if not config_path.exists():
        fail([f"Config file not found: {config_path}"])
        return {}
    try:
        with config_path.open("r", encoding="utf-8") as handle:
            return yaml.safe_load(handle) or {}
    except yaml.YAMLError as exc:
        fail([f"Could not parse YAML config '{config_path}': {exc}"])


def validate_inputs(config: dict, root: Path):
    errors = []

    # Check if validate_inputs block exists
    if "validate_inputs" not in config:
        errors.append("Missing 'validate_inputs' block in config.")
        return errors

    val_block = config["validate_inputs"]

    # Validate samples_tsv
    if "samples_tsv" not in val_block:
        errors.append("Missing 'samples_tsv' in 'validate_inputs' config.")
    else:
        samples_path = root / val_block["samples_tsv"]
        if not samples_path.exists():
            errors.append(f"Samples TSV not found at: {samples_path}")
        else:
            with samples_path.open("r", encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t")
                header = next(reader, None)
                if not header:
                    errors.append("Samples TSV is empty.")

    # Validate protein_table
    if "protein_table" not in val_block:
        errors.append("Missing 'protein_table' in 'validate_inputs' config.")
    else:
        protein_path = root / val_block["protein_table"]
        if not protein_path.exists():
            errors.append(f"Protein intensities table not found at: {protein_path}")

    return errors


def main():
    root = Path.cwd()
    config_arg = Path(sys.argv[1]) if len(sys.argv) > 1 else root / "config/config.yaml"

    print(f"{Colors.BLUE}Validating qChIP-MS Configuration...{Colors.ENDC}")
    config = load_config(config_arg)

    errors = validate_inputs(config, root)
    if errors:
        fail(errors)

    print(f"{Colors.GREEN}[CONFIG VALIDATION] OK{Colors.ENDC}")


if __name__ == "__main__":
    main()
