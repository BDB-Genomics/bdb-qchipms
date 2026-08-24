#!/usr/bin/env python3
"""
Unit and Integration Test Suite for bdb-qchipms Pipeline.
Tests configuration schema, sample metadata, peak intersection, and Snakemake DAG.
"""

import os
import shutil
import subprocess
import pytest
import pandas as pd
import yaml
from pathlib import Path


@pytest.fixture
def repo_root():
    return Path(__file__).resolve().parent.parent


def test_config_yaml_exists_and_valid(repo_root):
    config_path = repo_root / "config" / "config.yaml"
    assert config_path.exists(), "config/config.yaml is missing"

    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    assert "references" in cfg, "Missing 'references' section in config"
    assert "validate_inputs" in cfg, "Missing 'validate_inputs' section"
    assert "lfq_enrichment" in cfg, "Missing 'lfq_enrichment' section"
    assert "locus_annotation" in cfg, "Missing 'locus_annotation' section"
    assert "summary_report" in cfg, "Missing 'summary_report' section"


def test_samples_tsv_structure(repo_root):
    samples_path = repo_root / "config" / "samples.tsv"
    assert samples_path.exists(), "config/samples.tsv is missing"

    df = pd.read_csv(samples_path, sep="\t")
    required_cols = {"sample_id", "group", "bait_target"}
    assert required_cols.issubset(set(df.columns)), f"Missing required columns in samples.tsv: {required_cols - set(df.columns)}"
    assert len(df) >= 2, "samples.tsv must have at least 2 samples"
    assert "treatment" in df["group"].values, "Must contain at least one treatment group"
    assert "control" in df["group"].values, "Must contain at least one control group"


def test_references_exist(repo_root):
    chrom_sizes_path = repo_root / "data" / "references" / "hg38.chrom.sizes"
    assert chrom_sizes_path.exists(), "data/references/hg38.chrom.sizes is missing"
    assert chrom_sizes_path.stat().st_size > 0, "hg38.chrom.sizes is empty"

    uniprot_path = repo_root / "data" / "references" / "uniprot_human.fasta"
    if uniprot_path.exists():
        assert uniprot_path.stat().st_size > 0, "uniprot_human.fasta is empty"


def test_macs2_peaks_bed_valid(repo_root):
    peaks_path = repo_root / "data" / "macs2_peaks.bed"
    assert peaks_path.exists(), "data/macs2_peaks.bed is missing"
    assert peaks_path.stat().st_size > 0, "data/macs2_peaks.bed is empty"

    with open(peaks_path, "r") as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]

    assert len(lines) > 0, "No peak records found in macs2_peaks.bed"
    for line in lines:
        parts = line.split("\t")
        assert len(parts) >= 3, f"Invalid BED line: {line}"
        chrom, start, end = parts[0], int(parts[1]), int(parts[2])
        assert start >= 0, f"Negative start coordinate in {line}"
        assert end > start, f"End coordinate <= start in {line}"


def test_protein_intensities_structure(repo_root):
    intensities_path = repo_root / "data" / "protein_intensities.tsv"
    if intensities_path.exists():
        df = pd.read_csv(intensities_path, sep="\t")
        assert len(df.columns) >= 3, "protein_intensities.tsv must have protein ID and sample columns"
        assert len(df) > 0, "protein_intensities.tsv must contain data rows"


def test_snakefile_rules_and_smk_files(repo_root):
    snakefile = repo_root / "Snakefile"
    assert snakefile.exists(), "Snakefile is missing"
    
    rules_dir = repo_root / "rules"
    assert rules_dir.exists(), "rules/ directory is missing"
    
    required_smk = [
        "validate_inputs.smk",
        "lfq_enrichment_and_plot.smk",
        "locus_annotation.smk",
        "summary_report.smk"
    ]
    for smk in required_smk:
        p = rules_dir / smk
        assert p.exists(), f"Rule file {smk} is missing in rules/"
        with open(p, "r") as f:
            content = f.read()
            assert "rule " in content, f"No rule definition in {smk}"


def test_snakemake_dry_run(repo_root):
    # Check if snakemake is installed
    snakemake_bin = shutil.which("snakemake")
    if snakemake_bin:
        cmd = [snakemake_bin, "-n", "--quiet"]
        res = subprocess.run(cmd, cwd=repo_root, capture_output=True, text=True)
        assert res.returncode == 0, f"Snakemake dry-run failed:\n{res.stderr}"
    else:
        # If running in environment without snakemake binary, verify snakefile structure
        with open(repo_root / "Snakefile", "r") as f:
            snake_content = f.read()
        assert "rule all:" in snake_content
