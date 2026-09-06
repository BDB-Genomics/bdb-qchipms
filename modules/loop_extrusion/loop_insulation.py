"""
Subject: Statistical Genetics & Chromatin 3D Architecture (Scale 4: 1 Mb)
Core Concept: Cohesin Loop Extrusion, Convergent CTCF Anchors & Insulation
Foundational References:
- Sanborn et al., PNAS 112, E6456–E6465 (2015) [Chromatin extrusion model]
- Rao et al., Cell 159, 1665–1680 (2014) [3D map of the human genome at kilobase resolution]
- Fudenberg et al., Cell Reports 15, 2038–2049 (2016) [Formation of chromosomal domains by loop extrusion]

Methodology & Physical Equations:
1. Cohesin processively extrudes chromatin loops symmetrically until it encounters CTCF bound in a CONVERGENT orientation.
   - Convergent pair: Left CTCF (+) pointing Forward ---> <--- Right CTCF (-) pointing Reverse. (Strong Loop Anchor)
   - Tandem pair: ---> ---> or <--- <--- (Weak / Leaky Loop Anchor)
   - Divergent pair: <--- ---> (Non-insulating anchor)
2. Loop contact probability P(s) as a function of genomic distance s:
   P(s) proportional to s^(-gamma), where gamma ~ 1.0 within TADs, and gamma ~ 1.5 across TAD boundaries.
3. Boundary Insulation Score (IS):
   Measures the log2 ratio of intra-domain contacts vs inter-domain contacts crossing the boundary.
"""

from typing import Dict, List, Tuple

# Canonical CTCF 20-bp consensus core motif (JASPAR MA0139.1)
CTCF_CONSENSUS_CORE = "CCACCAGGGGGCGCTA"

def compute_ctcf_orientation_affinity(motif_sequence: str, strand: str) -> Dict[str, any]:
    """
    Evaluates CTCF motif binding presence and orientation vector with robust error handling.
    """
    try:
        if not motif_sequence or not isinstance(motif_sequence, str):
            score = 0.0
            motif = ""
        else:
            motif = motif_sequence.upper().strip()
            core_len = len(CTCF_CONSENSUS_CORE)
            if core_len == 0:
                score = 0.0
            else:
                matches = sum(1 for a, b in zip(motif, CTCF_CONSENSUS_CORE) if a == b)
                score = float(matches) / float(core_len)
    except Exception:
        score = 0.0
        motif = str(motif_sequence) if motif_sequence else ""

    valid_strand = strand if strand in ['+', '-'] else '+'

    return {
        'sequence': motif,
        'strand': valid_strand,
        'consensus_match_pct': round(score * 100.0, 1),
        'is_functional_ctcf': bool(score >= 0.70),
        'vector': '+ (Forward --->)' if valid_strand == '+' else '- (Reverse <---)'
    }

def evaluate_loop_anchor_pair(left_ctcf: Dict[str, any], right_ctcf: Dict[str, any], distance_bp: int) -> Dict[str, any]:
    """
    Determines loop extrusion stalling stability and contact insulation.
    """
    left_functional = left_ctcf['is_functional_ctcf']
    right_functional = right_ctcf['is_functional_ctcf']
    
    if not (left_functional and right_functional):
        return {
            'loop_status': 'NO_STABLE_LOOP (Broken Anchor Motif)',
            'insulation_strength': 'LEAKY / UNINSULATED',
            'contact_probability_boost': 1.0,
            'physical_mechanism': 'Cohesin extrudes past non-functional or deleted CTCF sites.'
        }
        
    # Check orientation rule
    is_convergent = (left_ctcf['strand'] == '+' and right_ctcf['strand'] == '-')
    is_tandem = (left_ctcf['strand'] == right_ctcf['strand'])
    is_divergent = (left_ctcf['strand'] == '-' and right_ctcf['strand'] == '+')
    
    if is_convergent:
        status = 'STABLE_CONVERGENT_LOOP_ANCHOR'
        insulation = 'STRONG_INSULATION (Domain Boundary Formed)'
        boost = 8.5 # High corner-peak contact enrichment
        mech = 'Cohesin halts at opposing CTCF N-termini, forming a high-frequency regulatory contact domain.'
    elif is_tandem:
        status = 'WEAK_TANDEM_ANCHOR'
        insulation = 'MODERATE_INSULATION'
        boost = 2.2
        mech = 'Partial pausing; significant extrusion leakage occurs past tandem motif.'
    else:
        status = 'DIVERGENT_NON_INSULATING'
        insulation = 'MINIMAL_INSULATION'
        boost = 1.1
        mech = 'Extrusion complex cannot stall at divergent motifs.'
        
    return {
        'distance_kb': round(distance_bp / 1000.0, 1),
        'left_anchor': left_ctcf['vector'],
        'right_anchor': right_ctcf['vector'],
        'loop_status': status,
        'insulation_strength': insulation,
        'contact_probability_boost': boost,
        'physical_mechanism': mech
    }

if __name__ == '__main__':
    print("--- Module 04 (3D Genome Architecture): Cohesin Loop Extrusion & CTCF ---")
    
    left_ctcf = compute_ctcf_orientation_affinity("CCACCAGGGGGCGCTA", "+")
    right_ctcf = compute_ctcf_orientation_affinity("CCACCAGGGGGCGCTA", "-")
    right_broken = compute_ctcf_orientation_affinity("CCTTAATTTTTTCCTA", "-") # Disruptive non-coding variant
    
    print("\n1. Convergent CTCF Pair (Intact 450 kb TAD):")
    res_convergent = evaluate_loop_anchor_pair(left_ctcf, right_ctcf, 450000)
    print(f"   • Status:      {res_convergent['loop_status']}")
    print(f"   • Insulation:  {res_convergent['insulation_strength']}")
    print(f"   • Mechanism:   {res_convergent['physical_mechanism']}")
    
    print("\n2. Disrupted Boundary Variant (Non-coding SNP in Right Anchor):")
    res_broken = evaluate_loop_anchor_pair(left_ctcf, right_broken, 450000)
    print(f"   • Status:      {res_broken['loop_status']}")
    print(f"   • Insulation:  {res_broken['insulation_strength']}")
    print(f"   • Mechanism:   {res_broken['physical_mechanism']}")
