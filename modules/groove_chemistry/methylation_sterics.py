"""
Subject: Epigenetic Chemistry (Scale 1: 0.34 nm)
Core Concept: 5-Methylcytosine (5mC) Steric Hindrance & Hydrophobic Exclusion
Foundational References:
- Klose & Bird, Trends in Biochem Sci 31, 89–97 (2006) [Genomic DNA methylation: the mark and its mediators]
- Yin et al., Science 356, eaaj2239 (2017) [Impact of cytosine methylation on TF-DNA binding specificity]

Methodology:
1. Methylation adds a hydrophobic -CH3 group at cytosine C5 position in the major groove.
2. In non-methyl-binding TFs, this causes severe steric clash (volume addition ~35 A^3) and
   alters the donor/acceptor barcode (HDAA -> MDAA), decreasing binding affinity by >10-fold.
3. In Methyl-CpG-binding proteins (MBDs, MeCP2), it forms a hydrophobic pocket stabilizing binding.
"""

from typing import Dict

def evaluate_cpg_methylation_effect(base_sequence: str, is_methylated: bool, tf_type: str = "STERIC_SENSITIVE") -> Dict[str, any]:
    """
    Evaluates major groove chemical alteration and binding free energy shift from 5mC methylation.
    """
    seq = base_sequence.upper()
    cpg_count = seq.count("CG")
    
    if not is_methylated or cpg_count == 0:
        return {
            'cpg_sites': cpg_count,
            'methylation_status': 'UNMETHYLATED (C5-H)',
            'major_groove_code': 'HDAA (Canonical Cytosine)',
            'delta_G_shift_kcal_mol': 0.0,
            'binding_impact': 'NEUTRAL / UNPERMANENT'
        }
        
    if tf_type == "STERIC_SENSITIVE":
        # Steric clash destabilizes typical pioneer/TF binding (e.g. CTCF, ETS factors)
        delta_g = +2.4 * cpg_count # kcal/mol penalty per methylated CpG
        impact = f'SEVERE_INHIBITION ({round(delta_g, 1)} kcal/mol steric barrier disrupting TF binding)'
    else:
        # Specialized reader (MeCP2 / MBD2)
        delta_g = -1.8 * cpg_count # kcal/mol stabilization
        impact = f'ENHANCED_RECRUITMENT (Hydrophobic interaction stabilizes binding by {round(delta_g, 1)} kcal/mol)'
        
    return {
        'cpg_sites': cpg_count,
        'methylation_status': 'METHYLATED (5-Methylcytosine C5-CH3)',
        'major_groove_code': 'MDAA (Methylated Cytosine Barcode)',
        'delta_G_shift_kcal_mol': round(delta_g, 2),
        'binding_impact': impact
    }

if __name__ == '__main__':
    print("--- Module: 5mC Epigenetic Chemistry ---")
    res_unmeth = evaluate_cpg_methylation_effect("CCCGGATCG", is_methylated=False)
    res_meth = evaluate_cpg_methylation_effect("CCCGGATCG", is_methylated=True, tf_type="STERIC_SENSITIVE")
    print(f"Unmethylated: {res_unmeth['major_groove_code']} | Impact: {res_unmeth['binding_impact']}")
    print(f"Methylated:   {res_meth['major_groove_code']} | Impact: {res_meth['binding_impact']}")
