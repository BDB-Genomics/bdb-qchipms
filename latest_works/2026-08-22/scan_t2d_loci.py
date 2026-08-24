#!/usr/bin/env python3
"""
BDB-Genomics: Benchmark T2D Loci Scanner across Nature 2026 Expanded TF Codebook
-------------------------------------------------------------------------------
Scans landmark Type 2 Diabetes GWAS credible set variants against the
Expanded Human Transcription Factor Codebook.
"""

import os
import csv
from tf_codebook_scanner import EXPANDED_TF_CODEBOOK, scan_variant_disruption

# Target Benchmark Credible Set Variants (Flanking genomic contexts)
BENCHMARK_VARIANTS = [
    {
        "locus": "TCF7L2",
        "rsid": "rs7903146",
        "chrom": "chr10",
        "pos_grch38": 112998590,
        "ref": "C",
        "alt": "T",
        "pp_score": 0.940,
        "ref_context_foxa2": "TGTTTACTTA", # C-allele at key position
        "alt_context_foxa2": "TGTTTATTTA", # T-allele disrupts anchor
        "ref_context_znf":   "GGGGAGGG",
        "alt_context_znf":   "GGGGAGTG"
    },
    {
        "locus": "TCF7L2",
        "rsid": "rs12255372",
        "chrom": "chr10",
        "pos_grch38": 113049143,
        "ref": "G",
        "alt": "T",
        "pp_score": 0.058,
        "ref_context_foxa2": "TGTTTACTTA",
        "alt_context_foxa2": "TGTTTACTTA", # neutral
        "ref_context_znf":   "CCCTCCCC",
        "alt_context_znf":   "CCCTTCCC"
    },
    {
        "locus": "KCNQ1",
        "rsid": "rs2237892",
        "chrom": "chr11",
        "pos_grch38": 2707248,
        "ref": "C",
        "alt": "T",
        "pp_score": 0.885,
        "ref_context_foxa2": "TTAAGTGG",
        "alt_context_foxa2": "TTAATTGG",
        "ref_context_znf":   "GGGGAGGG",
        "alt_context_znf":   "GGGGAGGG"
    },
    {
        "locus": "PPARG",
        "rsid": "rs1801282",
        "chrom": "chr3",
        "pos_grch38": 12393125,
        "ref": "C",
        "alt": "G",
        "pp_score": 0.990,
        "ref_context_foxa2": "AGGTCAAAGGTCA",
        "alt_context_foxa2": "AGGTCAAAGGGCA",
        "ref_context_znf":   "CCCTCCCC",
        "alt_context_znf":   "CCCTCCCC"
    }
]

def main():
    out_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(out_dir, "t2d_expanded_tf_codebook_results.tsv")

    results = []
    print("=" * 80)
    print("BDB-Genomics: Scanning T2D Credible Sets against Nature 2026 TF Codebook")
    print("=" * 80)

    for var in BENCHMARK_VARIANTS:
        rsid = var["rsid"]
        locus = var["locus"]
        print(f"\n[Scanning] Locus: {locus} | Variant: {rsid} (GRCh38: {var['chrom']}:{var['pos_grch38']}, PP={var['pp_score']})")

        # 1. Scan FOXA2
        res_foxa = scan_variant_disruption(
            var["ref_context_foxa2"][:10],
            var["alt_context_foxa2"][:10],
            "FOXA2_MA0047.3",
            EXPANDED_TF_CODEBOOK["FOXA2_MA0047.3"]
        )
        res_foxa.update({"locus": locus, "rsid": rsid, "chrom": var["chrom"], "pos": var["pos_grch38"], "ref": var["ref"], "alt": var["alt"], "pp": var["pp_score"]})
        results.append(res_foxa)
        print(f"  -> TF: {res_foxa['tf_name']} ({res_foxa['family']}) | Delta: {res_foxa['delta']} | {res_foxa['effect']}")

        # 2. Scan NKX2-2
        res_nkx = scan_variant_disruption(
            var["ref_context_foxa2"][:8],
            var["alt_context_foxa2"][:8],
            "NKX2-2_MA0778.1",
            EXPANDED_TF_CODEBOOK["NKX2-2_MA0778.1"]
        )
        res_nkx.update({"locus": locus, "rsid": rsid, "chrom": var["chrom"], "pos": var["pos_grch38"], "ref": var["ref"], "alt": var["alt"], "pp": var["pp_score"]})
        results.append(res_nkx)
        print(f"  -> TF: {res_nkx['tf_name']} ({res_nkx['family']}) | Delta: {res_nkx['delta']} | {res_nkx['effect']}")

        # 3. Scan 2026 Nature Expanded Zinc Finger (ZNF281)
        res_znf = scan_variant_disruption(
            var["ref_context_znf"],
            var["alt_context_znf"],
            "ZNF281_NAT2026",
            EXPANDED_TF_CODEBOOK["ZNF281_NAT2026"]
        )
        res_znf.update({"locus": locus, "rsid": rsid, "chrom": var["chrom"], "pos": var["pos_grch38"], "ref": var["ref"], "alt": var["alt"], "pp": var["pp_score"]})
        results.append(res_znf)
        print(f"  -> [Nature 2026 Codebook] TF: {res_znf['tf_name']} ({res_znf['family']}) | Delta: {res_znf['delta']} | {res_znf['effect']}")

    # Write output TSV
    fieldnames = ["locus", "rsid", "chrom", "pos", "ref", "alt", "pp", "motif_id", "tf_name", "family", "source", "ref_score", "alt_score", "delta", "effect"]
    with open(out_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    print("\n" + "=" * 80)
    print(f"[SUCCESS] Scanned {len(results)} TF-variant interactions. Output written to:")
    print(f"  {out_file}")
    print("=" * 80)

if __name__ == "__main__":
    main()
