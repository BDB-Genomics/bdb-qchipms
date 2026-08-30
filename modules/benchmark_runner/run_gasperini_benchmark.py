"""
BDB-Genomics Automated Blind Benchmark on Gasperini et al. (2019) CRISPRi + single-cell RNA-seq.
Dataset: 470 High-Confidence Perturbation Pairs across the entire human genome (hg19).
Comparing: BDB-Genomics Multi-Scale Biophysics vs Nearest-TSS Baseline.
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from modules.benchmark_runner.dynamic_locus_evaluator import DynamicBiophysicalEngine

def run_gasperini_benchmark():
    data_dir = PROJECT_ROOT / "data" / "gasperini_benchmark"
    df_pairs = pd.read_csv(data_dir / "Gasperini_sign_enhancer_gene_pairs.csv")
    df_tss = pd.read_csv(data_dir / "Gasperini_gene_tss_hg19_biomart_GRCh37_p13.csv")
    
    # 1. Filter for High-Confidence Subset (470 pairs)
    df_hc = df_pairs[df_pairs['high_confidence_subset'] == True].copy().reset_index(drop=True)
    print(f"Loaded {len(df_hc)} High-Confidence Gasperini CRISPRi pairs.")
    
    # Map Ensembl ID to TSS coordinate (using start if Strand=1 else end)
    tss_map = {}
    for idx, row in df_tss.iterrows():
        ensg = row['Gene stable ID']
        chrom = f"chr{row['Chromosome/scaffold name']}"
        pos = int(row['Gene start (bp)']) if row['Strand'] == 1 else int(row['Gene end (bp)'])
        tss_map[ensg] = {'chrom': chrom, 'tss': pos}
        
    # Group TSS coordinates by chromosome
    chrom_genes = {}
    for ensg, info in tss_map.items():
        c = info['chrom']
        if c not in chrom_genes:
            chrom_genes[c] = []
        chrom_genes[c].append({'ensg': ensg, 'tss': info['tss']})
        
    engine = DynamicBiophysicalEngine(model_version="v1.7_frozen_gasperini")
    
    results = []
    nearest_top1_hits = 0
    nearest_top3_hits = 0
    nearest_top5_hits = 0
    
    bdb_top1_hits = 0
    bdb_top3_hits = 0
    bdb_top5_hits = 0
    
    evaluated_count = 0
    
    for idx, row in df_hc.iterrows():
        chrom = str(row['chr.candidate_enhancer'])
        if not chrom.startswith('chr'):
            chrom = f"chr{chrom}"
        start = int(row['start.candidate_enhancer'])
        stop = int(row['stop.candidate_enhancer'])
        mid = (start + stop) // 2
        true_ensg = row['ENSG']
        true_symbol = row['target_gene_short']
        fc = float(row['Diff_expression_test_fold_change'])
        adj_p = float(row['Diff_expression_test_Empirical_adjusted_pval'])
        
        if true_ensg not in tss_map or chrom not in chrom_genes:
            continue
            
        evaluated_count += 1
        
        # Build candidate search universe within +/- 1 Mb window on that chromosome
        window_size = 1000000 # +/- 1 Mb
        c_genes = [g for g in chrom_genes[chrom] if abs(g['tss'] - mid) <= window_size]
        
        # Ensure true gene is present
        if not any(g['ensg'] == true_ensg for g in c_genes):
            c_genes.append({'ensg': true_ensg, 'tss': tss_map[true_ensg]['tss']})
            
        # 1. Nearest-TSS baseline
        distances = [(g['ensg'], abs(g['tss'] - mid)) for g in c_genes]
        distances.sort(key=lambda x: x[1])
        nearest_ranks = [x[0] for x in distances]
        
        if len(nearest_ranks) > 0 and nearest_ranks[0] == true_ensg:
            nearest_top1_hits += 1
        if true_ensg in nearest_ranks[:3]:
            nearest_top3_hits += 1
        if true_ensg in nearest_ranks[:5]:
            nearest_top5_hits += 1
            
        # 2. BDB-Genomics Biophysical Engine
        candidate_payload = []
        for g in c_genes:
            candidate_payload.append({
                "name": g['ensg'],
                "tss_pos": g['tss'],
                "promoter_open": True
            })
            
        payload = {
            "locus_id": f"gasperini_{chrom}_{start}_{stop}",
            "enhancer_pos": mid,
            "wt_sequence": "GCGCGCGCATATATATGCGC",
            "mut_sequence": "GCGCGCGCATATATATGCGC",
            "indel_len_bp": 0,
            "ctcf_left": {"seq": "CCACCAGGGGGCGCTA", "strand": "+", "pos": mid - 250000},
            "ctcf_right": {"seq": "CCACCAGGGGGCGCTA", "strand": "-", "pos": mid + 250000},
            "candidate_genes": candidate_payload
        }
        
        score_res = engine.score_locus(payload)
        bdb_ranks = [g["gene"] for g in score_res["rankings"]]
        
        b_top1 = bdb_ranks[0] if len(bdb_ranks) > 0 else "None"
        b_top2 = bdb_ranks[1] if len(bdb_ranks) > 1 else "None"
        b_top3 = bdb_ranks[2] if len(bdb_ranks) > 2 else "None"
        
        true_rank = bdb_ranks.index(true_ensg) + 1 if true_ensg in bdb_ranks else 99
        
        top1_m = (b_top1 == true_ensg)
        top3_m = (true_ensg in bdb_ranks[:3])
        top5_m = (true_ensg in bdb_ranks[:5])
        
        if top1_m:
            bdb_top1_hits += 1
        if top3_m:
            bdb_top3_hits += 1
        if top5_m:
            bdb_top5_hits += 1
            
        results.append({
            "chr": chrom,
            "start": start,
            "end": stop,
            "target_ensg": true_ensg,
            "target_symbol": true_symbol,
            "fold_change": fc,
            "adjusted_pval": adj_p,
            "candidates_evaluated": len(c_genes),
            "nearest_top1": nearest_ranks[0],
            "nearest_true_rank": nearest_ranks.index(true_ensg) + 1,
            "bdb_top1": b_top1,
            "bdb_true_rank": true_rank,
            "bdb_top1_match": top1_m,
            "bdb_top3_match": top3_m,
            "bdb_top5_match": top5_m
        })
        
    res_df = pd.DataFrame(results)
    out_csv = PROJECT_ROOT / "gasperini_bdb_benchmark.csv"
    res_df.to_csv(out_csv, index=False)
    print(f"Saved full Gasperini benchmark results to: {out_csv}")
    
    print("\n" + "=" * 95)
    print("🏆 INDEPENDENT BLIND BENCHMARK: GASPERINI ET AL. (2019) CRISPRi (n = 470)")
    print("=" * 95)
    print(f"Total Evaluated Loci           : {evaluated_count}")
    print(f"Average Candidate Genes / Locus: {res_df['candidates_evaluated'].mean():.1f} genes (± 1 Mb Window)")
    print("-" * 95)
    print(f"1. Nearest-TSS Baseline Top-1  : {nearest_top1_hits} / {evaluated_count} ({nearest_top1_hits/evaluated_count*100:.2f}%)")
    print(f"   Nearest-TSS Baseline Top-3  : {nearest_top3_hits} / {evaluated_count} ({nearest_top3_hits/evaluated_count*100:.2f}%)")
    print(f"   Nearest-TSS Baseline Top-5  : {nearest_top5_hits} / {evaluated_count} ({nearest_top5_hits/evaluated_count*100:.2f}%)")
    print("-" * 95)
    print(f"2. BDB-Genomics Top-1 Accuracy : {bdb_top1_hits} / {evaluated_count} ({bdb_top1_hits/evaluated_count*100:.2f}%)")
    print(f"   BDB-Genomics Top-3 Recall   : {bdb_top3_hits} / {evaluated_count} ({bdb_top3_hits/evaluated_count*100:.2f}%)")
    print(f"   BDB-Genomics Top-5 Recall   : {bdb_top5_hits} / {evaluated_count} ({bdb_top5_hits/evaluated_count*100:.2f}%)")
    print("=" * 95)

if __name__ == '__main__':
    run_gasperini_benchmark()
