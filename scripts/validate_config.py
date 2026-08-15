#!/usr/bin/env python3
"""Validate the qChIP-MS config and input files before Snakemake builds the DAG.

Applies HPC streaming checks (O(1) memory, fast single-line header inspection)
to ensure all input tables and references exist and are valid before computation.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Any, NoReturn
import yaml


class Colors:
    HEADER = "\033[95m"
    BLUE = "\033[94m"
    CYAN = "\033[96m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"


def fail(errors: list[str]) -> NoReturn:
    """Print error block and exit with non-zero status."""
    print(
        f"\n{Colors.BOLD}{Colors.RED}┏━ CONFIGURATION VALIDATION FAILED ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓{Colors.ENDC}"
    )
    for m in errors:
        print(f"    • {m}")
    print(
        f"{Colors.BOLD}{Colors.RED}┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛{Colors.ENDC}\n"
    )
    sys.exit(1)


def load_config(config_path: Path) -> dict[str, Any]:
    """Load and parse the YAML configuration file safely."""
    if not config_path.exists():
        fail([f"Config file not found: {config_path}"])
        return {}
    try:
        with config_path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
            return data if isinstance(data, dict) else {}
    except yaml.YAMLError as exc:
        fail([f"Could not parse YAML config '{config_path}': {exc}"])
        return {}


def check_file_non_empty(file_path: Path, label: str, errors: list[str]) -> bool:
    """O(1) file existence and size check using OS stat."""
    if not file_path.exists():
        errors.append(f"{label} not found at: {file_path}")
        return False
    if file_path.stat().st_size == 0:
        errors.append(f"{label} at {file_path} is empty (0 bytes).")
        return False
    return True


def validate_samples_tsv(samples_path: Path, errors: list[str]) -> None:
    """Fast header stream verification for samples manifest."""
    if not check_file_non_empty(samples_path, "Samples TSV", errors):
        return

    with samples_path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            errors.append(f"Samples TSV ({samples_path}) has no header line.")
            return

        expected = {"sample_id", "group"}
        found = {h.strip().lower() for h in header}
        if not expected.issubset(found):
            errors.append(
                f"Samples TSV ({samples_path}) missing required columns: {expected - found}"
            )


def validate_protein_table(protein_path: Path, errors: list[str]) -> None:
    """O(1) memory inspection of protein intensities header."""
    if not check_file_non_empty(protein_path, "Protein table", errors):
        return

    # Read only the first line into memory to verify intensity matrix structure
    with protein_path.open("r", encoding="utf-8") as f:
        header_line = f.readline()
        if not header_line:
            errors.append(f"Protein table ({protein_path}) is empty.")
            return

        columns = [c.strip() for c in header_line.split("\t")]
        has_protein_id = any("protein" in c.lower() for c in columns)
        has_ctrl = any("igg" in c.lower() or "control" in c.lower() for c in columns)
        has_treat = any(
            "trf2" in c.lower() or "terf2" in c.lower() or "treat" in c.lower()
            for c in columns
        )

        if not has_protein_id:
            errors.append(
                f"Protein table ({protein_path}) missing Protein ID identifier column."
            )
        if not (has_ctrl and has_treat):
            errors.append(
                f"Protein table ({protein_path}) must contain both treatment and control intensity columns."
            )


def validate_inputs(config: dict[str, Any], root: Path) -> list[str]:
    """Validate all required config sections and reference files."""
    errors: list[str] = []

    if "validate_inputs" not in config:
        errors.append("Missing 'validate_inputs' section in config.yaml.")
        return errors

    val_block = config["validate_inputs"]

    # 1. Validate samples TSV
    if "samples_tsv" not in val_block:
        errors.append("Missing 'samples_tsv' in 'validate_inputs' config.")
    else:
        validate_samples_tsv(root / val_block["samples_tsv"], errors)

    # 2. Validate protein table
    if "protein_table" not in val_block:
        errors.append("Missing 'protein_table' in 'validate_inputs' config.")
    else:
        validate_protein_table(root / val_block["protein_table"], errors)

    # 3. Validate reference files if declared
    if "references" in config and isinstance(config["references"], dict):
        refs = config["references"]
        if "chrom_sizes" in refs:
            chrom_path = root / refs["chrom_sizes"]
            check_file_non_empty(chrom_path, "Chromosome sizes reference", errors)

    return errors


def main() -> None:
    root = Path.cwd()
    config_arg = Path(sys.argv[1]) if len(sys.argv) > 1 else root / "config/config.yaml"

    print(f"{Colors.BLUE}Validating qChIP-MS Configuration...{Colors.ENDC}")
    config = load_config(config_arg)

    errors = validate_inputs(config, root)
    if errors:
        fail(errors)

    print(f"{Colors.GREEN}[CONFIG VALIDATION] OK — All inputs and matrices verified.{Colors.ENDC}")


if __name__ == "__main__":
    main()
