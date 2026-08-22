#!/usr/bin/env python3
"""
Unit and Integration Test Suite for bdb-qchipms Pipeline.
Tests configuration schema, sample metadata, peak intersection, and Snakemake DAG.
"""

import os
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
    uniprot_path = repo_root / "data" / "references" / "uniprot_human.fasta"
    chrom_sizes_path = repo_root / "data" / "references" / "hg38.chrom.sizes"

    assert uniprot_path.exists(), "data/references/uniprot_human.fasta is missing"
    assert chrom_sizes_path.exists(), "data/references/hg38.chrom.sizes is missing"
    assert uniprot_path.stat().st_size > 0, "uniprot_human.fasta is empty"
    assert chrom_sizes_path.stat().st_size > 0, "hg38.chrom.sizes is empty"


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


def test_snakemake_dry_run(repo_root):
    # Try generic snakemake first, then fallback to conda runner if available
    cmd = ["snakemake", "-n", "--quiet"]
    try:
        res = subprocess.run(cmd, cwd=repo_root, capture_output=True, text=True)
    except FileNotFoundError:
        cmd = ["conda", "run", "-n", "snakemake_runner", "snakemake", "-n", "--quiet"]
        res = subprocess.run(cmd, cwd=repo_root, capture_output=True, text=True)
    assert res.returncode == 0, f"Snakemake dry-run failed:\n{res.stderr}"


def test_output_artifacts_exist(repo_root):
    enriched_tsv = repo_root / "results" / "enrichment" / "qchip_ms_enriched_proteins.tsv"
    volcano_png = repo_root / "results" / "plots" / "qchip_ms_volcano.png"
    annotated_bed = repo_root / "results" / "annotation" / "qchip_ms_locus_annotated.bed"
    summary_html = repo_root / "results" / "reports" / "qchip_ms_summary.html"

    assert enriched_tsv.exists(), "Enriched proteins TSV missing"
    assert volcano_png.exists(), "Volcano plot PNG missing"
    assert annotated_bed.exists(), "Annotated BED missing"
    assert summary_html.exists(), "Summary report HTML missing"
