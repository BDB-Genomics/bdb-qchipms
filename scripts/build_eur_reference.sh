#!/usr/bin/env bash
# Build EUR-only PLINK reference from 1000G high-coverage hg38 panel.
# Source: PLINK2 resources mirror of Byrska-Bishop et al. 2022 callset (3202 samples).
# Output: data/atlas/ref_panel/{eur_autosomes_wgs.bed,bim,fam} + EUR allele frequencies.
set -euo pipefail
cd "$(dirname "$0")/.."

REF=data/atlas/ref_panel
TOOLS=tools/plink2

echo "[1/4] Decompressing pvar..."
zstd -d -f "$REF/all_hg38.pvar.zst" -o "$REF/all_hg38.pvar"

echo "[2/4] Extracting EUR sample list..."
awk 'NR>1 && $5=="EUR" {print $1}' "$REF/hg38_corrected.psam" > "$REF/eur_samples.txt"
wc -l "$REF/eur_samples.txt"

echo "[3/4] Building EUR-only autosome BED (single pass over pgen)..."
"$TOOLS" \
  --pfile "$REF/all_hg38" \
  --keep "$REF/eur_samples.txt" \
  --chr 1-22 \
  --snps-only just-acgt \
  --max-alleles 2 \
  --memory 2000 \
  --threads 8 \
  --make-bed \
  --out "$REF/eur_autosomes_wgs"

echo "[4/4] Computing EUR allele frequencies..."
"$TOOLS" \
  --bfile "$REF/eur_autosomes_wgs" \
  --memory 2000 \
  --threads 8 \
  --freq \
  --out "$REF/eur_af"

echo "DONE: $REF/eur_autosomes_wgs.{bed,bim,fam}"
