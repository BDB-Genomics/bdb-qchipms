from rank_candidates_engine import LOCI_DEFINITIONS, evaluate_locus_candidates

def main():
    print("=" * 145)
    print("TABLE 1: 8q24.21 SCREEN AUDIT (Fulco et al. 2016 Science e1-e4 + rs6983267 Colorectal Enhancer)")
    print("=" * 145)
    header = ["Hit_ID", "Genomic_Coord_GRCh38", "1D_Nearest_TSS", "Dist_Nearest", "BDB_Rank1_Target", "Dist_Target", "Contact_Prior", "BDB_Score", "Comp_Share", "Triage_Decision"]
    print(f"{header[0]:<12} | {header[1]:<24} | {header[2]:<16} | {header[3]:<12} | {header[4]:<16} | {header[5]:<12} | {header[6]:<13} | {header[7]:<10} | {header[8]:<10} | {header[9]}")
    print("-" * 145)

    for h in LOCI_DEFINITIONS[:5]:
        res = evaluate_locus_candidates(h)
        nearest = res["nearest_1d"]
        
        # For rs698 in colorectal context:
        # Candidate CCAT2 is proximal lncRNA; MYC is distal oncogene.
        # If auditing protein-coding targets or multi-gene models:
        top = res["rank1_bdb"]
        
        c_str = f"{h['chrom']}:{h['start']}-{h['stop']}" if h['start'] != h['stop'] else f"{h['chrom']}:{h['start']}"
        dist_near = f"{nearest['dist_tss']:,} bp"
        dist_top = f"{top['dist_tss']:,} bp"
        dec = f"PROMOTE ({top['gene']}) / DEMOTE ({nearest['gene']})" if top["gene"] != nearest["gene"] else f"CONCORDANT ({top['gene']})"
        print(f"{h['id']:<12} | {c_str:<24} | {nearest['gene']:<16} | {dist_near:<12} | {top['gene']:<16} | {dist_top:<12} | {top['p_contact']:<13.4f} | {top['raw_score']:<10.4f} | {top['competitive_share']*100:<9.1f}% | {dec}")

    print("\n" + "=" * 145)
    print("TABLE 2: MULTI-PROMOTER COMPETITION & LONG-RANGE LOCI (GWAS / CRISPR BENCHMARKS)")
    print("=" * 145)
    print(f"{header[0]:<12} | {header[1]:<24} | {header[2]:<16} | {header[3]:<12} | {header[4]:<16} | {header[5]:<12} | {header[6]:<13} | {header[7]:<10} | {header[8]:<10} | {header[9]}")
    print("-" * 145)

    for h in LOCI_DEFINITIONS[5:]:
        res = evaluate_locus_candidates(h)
        nearest = res["nearest_1d"]
        top = res["rank1_bdb"]
        c_str = f"{h['chrom']}:{h['start']}"
        dist_near = f"{nearest['dist_tss']:,} bp"
        dist_top = f"{top['dist_tss']:,} bp"
        dec = f"PROMOTE ({top['gene']}) / DEMOTE ({nearest['gene']})" if top["gene"] != nearest["gene"] else f"CONCORDANT ({top['gene']})"
        print(f"{h['id']:<12} | {c_str:<24} | {nearest['gene']:<16} | {dist_near:<12} | {top['gene']:<16} | {dist_top:<12} | {top['p_contact']:<13.4f} | {top['raw_score']:<10.4f} | {top['competitive_share']*100:<9.1f}% | {dec}")

if __name__ == '__main__':
    main()
