#!/bin/bash
# ==============================================================================
# BDB-Genomics: Master 2026 Multi-Omic & Epigenomic Suite Execution Engine
# ==============================================================================
set -euo pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$DIR"

echo "========================================================================"
echo "BDB-GENOMICS: EXECUTING FULL 2026 MULTI-OMIC & EPIGENOMICS FRONTIER SUITE"
echo "Date: $(date -u)"
echo "Directory: $DIR"
echo "========================================================================"

# 1. Run Nature 2026 Expanded TF Codebook Scanner
echo -e "\n>>> [1/5] Running Nature 2026 Expanded Human TF Codebook Scanner..."
python3 tf_codebook_scanner.py
python3 scan_t2d_loci.py

# 2. Run ENCODE-rE2G Enhancer-to-Gene Linker (Engreitz Nature 2026)
echo -e "\n>>> [2/5] Running ENCODE-rE2G Enhancer-to-Gene Linking Engine..."
python3 encode_re2g_linker.py

# 3. Run SpliceImpactR Alternative Splicing & PPI Rewirer (Fiszbein NAR 2026)
echo -e "\n>>> [3/5] Running SpliceImpactR Isoform & PPI Rewiring Module..."
python3 splice_impact_annotator.py

# 4. Run Dark Proteome & Alternate RNA Decoding Database Generator (Slavov Nature 2026)
echo -e "\n>>> [4/5] Running Dark Proteome Search Database Generator..."
python3 dark_proteome_search_generator.py

# 5. Run CellOracle In Silico GRN Perturbation Engine (Morris Nature 2023)
echo -e "\n>>> [5/5] Running CellOracle In Silico GRN Perturbation Engine..."
python3 celloracle_in_silico_perturb.py

# 6. Run UpSet Multi-Omic Intersection Matrix Visualizer
echo -e "\n>>> [6/6] Generating Multi-Omic UpSet Intersection Matrix..."
python3 multiomic_upset_visualizer.py

# 7. Generate Cryptographic Checksums
echo -e "\n>>> Cryptographically hashing all artifacts with SHA-256..."
sha256sum *.py *.tsv *.fasta *.md > checksums.sha256 2>/dev/null || sha256sum *.py *.tsv *.fasta > checksums.sha256

echo "========================================================================"
echo "[SUCCESS] Complete 2026 Multi-Omic Frontier Suite executed successfully!"
echo "========================================================================"
cat checksums.sha256
