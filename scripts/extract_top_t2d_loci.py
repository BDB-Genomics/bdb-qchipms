#!/usr/bin/env python3
"""
BDB-Genomics: Top T2D Loci Extractor & Clumper (Zero 3rd-party dependencies)
----------------------------------------------------------------------------
Streams through 5.05M DIAGRAM summary statistics, filters genome-wide
significant hits (P < 5e-8), and clumps them into top non-overlapping 
regional loci (1 Mb windows) ranked by statistical significance.
Uses pure standard library (gzip, csv, math).
"""

import gzip
import csv
import os
import sys

def clump_loci(gwas_path, p_cutoff=5e-8, window_bp=500000, top_n=60):
    print("=" * 80)
    print(f"BDB-Genomics: Extracting Top {top_n} T2D Loci from DIAGRAM GWAS (5.05M vars)")
    print("=" * 80)
    
    significant_vars = []
    
    print(f"[1/3] Streaming through {gwas_path} for P < {p_cutoff}...")
    with gzip.open(gwas_path, "rt") as f:
        header = f.readline().strip().split("\t")
        col_map = {name: i for i, name in enumerate(header)}
        
        chrom_idx = col_map.get("hm_chrom", col_map.get("chromosome"))
        pos_idx = col_map.get("hm_pos", col_map.get("base_pair_location"))
        rsid_idx = col_map.get("hm_rsid", col_map.get("variant_id"))
        p_idx = col_map.get("p_value")
        beta_idx = col_map.get("beta", col_map.get("hm_beta"))
        se_idx = col_map.get("standard_error")
        ea_idx = col_map.get("hm_effect_allele", col_map.get("effect_allele"))
        oa_idx = col_map.get("hm_other_allele", col_map.get("other_allele"))
        eaf_idx = col_map.get("hm_effect_allele_frequency", col_map.get("effect_allele_frequency"))

        count = 0
        for line in f:
            count += 1
            if count % 1000000 == 0:
                print(f"  Processed {count:,} variants... found {len(significant_vars):,} significant hits")
                
            parts = line.strip().split("\t")
            try:
                pval_str = parts[p_idx]
                if not pval_str or pval_str == "NA":
                    continue
                pval = float(pval_str)
                if pval <= p_cutoff:
                    chrom = str(parts[chrom_idx]).replace("chr", "")
                    pos = int(parts[pos_idx])
                    rsid = parts[rsid_idx]
                    beta = parts[beta_idx]
                    se = parts[se_idx]
                    ea = parts[ea_idx]
                    oa = parts[oa_idx]
                    eaf = parts[eaf_idx]
                    
                    significant_vars.append({
                        "rsid": rsid,
                        "chrom": chrom,
                        "pos": pos,
                        "pval": pval,
                        "beta": beta,
                        "se": se,
                        "effect_allele": ea,
                        "other_allele": oa,
                        "eaf": eaf
                    })
            except (ValueError, IndexError):
                continue

    print(f"\n[2/3] Total genome-wide significant variants found: {len(significant_vars):,}")
    
    # Sort by p-value ascending
    significant_vars.sort(key=lambda x: x["pval"])
    
    # Distance-based clumping
    print(f"\n[3/3] Performing distance clumping (+/- {window_bp/1000:.0f} kb window)...")
    clumped_loci = []
    claimed_regions = [] # (chrom, start, end)
    
    for row in significant_vars:
        c = row["chrom"]
        p = row["pos"]
        
        # Check if already covered by an existing top locus
        is_claimed = False
        for (cl_chrom, cl_start, cl_end) in claimed_regions:
            if c == cl_chrom and (cl_start <= p <= cl_end):
                is_claimed = True
                break
                
        if not is_claimed:
            win_start = max(1, p - window_bp)
            win_end = p + window_bp
            claimed_regions.append((c, win_start, win_end))
            
            locus_num = len(clumped_loci) + 1
            clumped_loci.append({
                "locus_id": f"T2D_LOCUS_{locus_num:03d}",
                "lead_rsid": row["rsid"],
                "chrom": f"chr{c}",
                "lead_pos": p,
                "window_start": win_start,
                "window_end": win_end,
                "lead_pval": row["pval"],
                "lead_beta": row["beta"],
                "lead_se": row["se"],
                "lead_effect_allele": row["effect_allele"],
                "lead_other_allele": row["other_allele"],
                "lead_eaf": row["eaf"]
            })
            
            if len(clumped_loci) >= top_n:
                break
                
    out_dir = "data/atlas"
    os.makedirs(out_dir, exist_ok=True)
    out_tsv = os.path.join(out_dir, f"t2d_top{top_n}_loci_manifest.tsv")
    
    fieldnames = [
        "locus_id", "lead_rsid", "chrom", "lead_pos", 
        "window_start", "window_end", "lead_pval", "lead_beta", 
        "lead_se", "lead_effect_allele", "lead_other_allele", "lead_eaf"
    ]
    
    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for locus in clumped_loci:
            writer.writerow(locus)
            
    print("\n" + "=" * 80)
    print(f"[SUCCESS] Clumped Top {len(clumped_loci)} Loci across the human genome!")
    print(f"Manifest written to: {out_tsv}")
    print("=" * 80)
    print("\nTop 10 Lead Loci Preview:")
    for locus in clumped_loci[:10]:
        print(f"  {locus['locus_id']}: {locus['lead_rsid']} | {locus['chrom']}:{locus['lead_pos']} (P={locus['lead_pval']:.2e}, beta={locus['lead_beta']})")
    
    return out_tsv

if __name__ == "__main__":
    gwas_file = "data/incoming/diagram_t2d_eur_gwas_catalog_raw.tsv.gz"
    if not os.path.exists(gwas_file):
        print(f"Error: {gwas_file} not found!", file=sys.stderr)
        sys.exit(1)
    clump_loci(gwas_file, p_cutoff=5e-8, window_bp=500000, top_n=60)
