# BDB-Genomics Blind Benchmark Suite

This directory is strictly isolated for running pre-registered, untainted blind evaluations of the BDB-Genomics biophysical engine.

## Directory Structure

- `data/`: Raw input datasets (e.g. held-out CRISPRi enhancer-gene screens, full human gene TSS annotations). True target labels must be hidden during prediction runs.
- `scripts/`: Benchmark runners, blind inference pipelines, and scoring harnesses.
- `results/`: Frozen prediction CSVs, scoring matrices, and statistical evaluation reports.

## Blind Protocol Rules

1. **Pre-Registration:** Lock and commit the model code and scoring weights before running inference.
2. **Hidden Labels:** Input to the scoring engine must only contain genomic coordinates (enhancer position, candidate window, promoter coordinates). Ground truth target gene identities must be withheld.
3. **Genome-Wide Candidate Universe:** Candidate genes must be drawn from the complete human genome annotation (all ~20,000 genes within the genomic search window), avoiding artificially restricted panels.
4. **Single-Shot Execution:** Predictions are generated once, timestamped, committed to git, and evaluated without parameter post-tuning.
