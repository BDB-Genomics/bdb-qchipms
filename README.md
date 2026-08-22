# bdb-qchipms: Quantitative Locus-Specific Chromatin Proteomics & Epigenomics Pipeline

[![CI](https://github.com/BDB-Genomics/bdb-qchipms/actions/workflows/ci.yml/badge.svg)](https://github.com/BDB-Genomics/bdb-qchipms/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/Snakemake-%E2%89%A57.0-brightgreen.svg)](https://snakemake.github.io)
[![Python: 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)

**Author & Maintainer:** Himanshu Bhandary (BDB-Genomics Research Consortium)  
**Correspondence:** `2032ushimanshu@gmail.com` | [GitHub](https://github.com/BDB-Genomics)

---

## 📌 Overview

`bdb-qchipms` is an automated, reproducible **Snakemake workflow** designed to characterize the **locus-specific protein composition** and transcription factor interactomes bound at designated genomic loci.

By integrating chromatin immunoprecipitation, label-free quantitative mass spectrometry (LFQ), and post-GWAS functional genomics, the pipeline automates raw data ingestion, statistical normalization, quality control, differential enrichment modeling, and interactive reporting.

```
                  ┌───────────────────────────────┐
                  │ Raw Mass-Spec LFQ & BAM Data  │
                  └──────────────┬────────────────┘
                                 │
                                 ▼
                  ┌───────────────────────────────┐
                  │ Data Ingestion & Sample QC    │
                  │ (FDR < 0.01, Missing Impute)  │
                  └──────────────┬────────────────┘
                                 │
                                 ▼
                  ┌───────────────────────────────┐
                  │ Differential Protein Complex  │
                  │ Enrichment Modeling (limma/LFQ│
                  └──────────────┬────────────────┘
                                 │
                                 ▼
                  ┌───────────────────────────────┐
                  │ Multi-Omic Integration &      │
                  │ Interactive Summary Reports   │
                  └───────────────────────────────┘
```

---

## 🔬 Core Features

* **Automated Workflow Orchestration:** Full DAG-based Snakemake pipeline with rule-level modularity and dependency tracking.
* **Rigorous Statistical Quality Control:** Filtering by peptide count, false discovery rate (FDR < 0.01), and log2 fold-change cutoffs against non-targeting controls (IgG / input).
* **Volcano & Interactome Visualization:** Automated generation of high-resolution volcano plots, correlation heatmaps, and protein-protein interaction subgraphs.
* **Containerized & Reproducible:** Fully specified Conda environments and Docker container specifications for zero-drift execution across workstations and HPC clusters.

---

## 🚀 Quick Start

### 1. Prerequisites
* [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/)
* [Snakemake](https://snakemake.readthedocs.io/) (>= 7.0)

### 2. Installation
```bash
git clone https://github.com/BDB-Genomics/bdb-qchipms.git
cd bdb-qchipms
mamba env create -f environment.yml
conda activate bdb-qchipms
```

### 3. Dry-Run Verification
Verify workflow DAG and rule dependencies without executing computation:
```bash
snakemake -n --cores 4
```

### 4. Pipeline Execution
Run the complete workflow using local cores or cluster submission (SLURM):
```bash
# Local execution
snakemake --use-conda --cores 8

# Cluster execution (SLURM profile)
snakemake --profile config/slurm --use-conda
```

---

## 📂 Repository Structure

```
bdb-qchipms/
├── config/                  # Pipeline parameters, sample manifests, and paths
│   └── config.yaml
├── workflow/ / rules/       # Modular Snakemake rules
│   ├── lfq_enrichment_and_plot.smk
│   └── summary_report.smk
├── scripts/                 # Core Python and R statistical analysis modules
│   ├── intersect_peaks.py
│   └── summary_report.py
├── tests/                   # Automated pytest and workflow validation suite
│   └── test_pipeline.py
├── Dockerfile               # Production container definition
├── environment.yml          # Top-level Conda dependency specification
└── Snakefile                # Master pipeline execution entrypoint
```

---

## 🧪 Testing

Run the automated test suite with pytest:
```bash
pytest tests/ -v
```

---

## 📜 Citation & License

This project is licensed under the **MIT License**.

If you use `bdb-qchipms` in your research, please cite:
```bibtex
@software{bhandary2026bdbqchipms,
  author = {Himanshu Bhandary},
  title = {bdb-qchipms: Quantitative Locus-Specific Chromatin Proteomics Pipeline},
  year = {2026},
  url = {https://github.com/BDB-Genomics/bdb-qchipms}
}
```
