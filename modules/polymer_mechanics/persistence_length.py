"""
Subject: General Physics & Polymer Mechanics (Scale 2: 10 nm)
Core Concept: Worm-Like Chain (WLC) Model & Sequence-Dependent Persistence Length
Foundational References:
- Marko & Siggia, Macromolecules 28, 8759–8770 (1995) [Worm-Like Chain elasticity]
- Bustamante et al., Science 265, 1599–1600 (1994) [Entropic elasticity of DNA]
- Nelson et al., Nature 325, 411–414 (1987) [A-tract curvature & persistence length]

Methodology & Physical Equations:
1. Interpolation formula for WLC force-extension relationship (Marko-Siggia):
   F(x) = (k_B * T / L_p) * [ 1 / (4 * (1 - x / L_0)^2) - 1/4 + x / L_0 ]
   where L_p is the persistence length, L_0 is the contour length, and x is extension.
2. Sequence-dependent persistence length:
   - Canonical random B-DNA: L_p ~ 50.0 nm (~150 bp)
   - Poly(dA:dT) rigid tracts (>= 4 bp): L_p increases up to 85.0–110.0 nm due to bifurcated H-bonds and propeller twist.
   - Flexible alternating Pyrimidine-Purine steps (e.g. CA/TG): L_p decreases to 35.0–40.0 nm.
3. Bending free energy for a segment of length L bent to radius R:
   E_bend = (k_B * T * L_p * L) / (2 * R^2)
"""

import math
from typing import Dict

# Physical Constants at 310.15 K (37 C)
KB = 1.380649e-23         # J/K
TEMP = 310.15             # Kelvin
KBT_JOULES = KB * TEMP    # ~4.28e-21 Joules
BP_RISE_NM = 0.34         # nm rise per base pair
JOULES_TO_KCAL_MOL = 6.02214076e23 / 4184.0 # Conversion factor

# Sequence dinucleotide rigidity offsets (nm offset relative to 50 nm baseline)
DINUCLEOTIDE_LP_MODIFIERS = {
    'AA': +12.0, 'TT': +12.0, 'AT': +6.0, 'TA': -4.0, # A-tracts are stiff; TA steps are flexible
    'GG': +2.0,  'CC': +2.0,  'GC': +4.0, 'CG': -6.0, # CG steps are flexible hinges
    'CA': -8.0,  'TG': -8.0,  'AC': -2.0, 'GT': -2.0, # CA/TG steps have high roll angle flexibility
    'GA': +1.0,  'TC': +1.0,  'AG': +1.0, 'CT': +1.0
}

def calculate_sequence_persistence_length(sequence: str) -> Dict[str, float]:
    """
    Computes sequence-dependent persistence length (L_p in nm and bp)
    and mechanical bending energy to wrap around a standard 10 nm radius (R=4.5 nm nucleosome core).
    """
    seq = sequence.upper()
    if len(seq) < 2:
        return {'L_p_nm': 50.0, 'L_p_bp': 147.0, 'bending_energy_kcal_mol': 0.0}
    
    base_lp = 50.0 # nm (canonical B-DNA baseline)
    total_mod = 0.0
    
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        total_mod += DINUCLEOTIDE_LP_MODIFIERS.get(dinuc, 0.0)
    
    avg_mod = total_mod / (len(seq) - 1)
    effective_lp_nm = max(25.0, min(120.0, base_lp + avg_mod))
    effective_lp_bp = effective_lp_nm / BP_RISE_NM
    
    # Bending energy for this sequence segment bent around nucleosome core radius R = 4.5 nm (45 A)
    contour_len_nm = len(seq) * BP_RISE_NM
    r_bend_nm = 4.5
    e_bend_joules = (KBT_JOULES * effective_lp_nm * 1e-9 * contour_len_nm * 1e-9) / (2.0 * (r_bend_nm * 1e-9)**2)
    e_bend_kcal_mol = e_bend_joules * JOULES_TO_KCAL_MOL
    
    return {
        'length_bp': len(seq),
        'contour_length_nm': round(contour_len_nm, 2),
        'persistence_length_nm': round(effective_lp_nm, 2),
        'persistence_length_bp': round(effective_lp_bp, 1),
        'nucleosome_bending_energy_kcal_mol': round(e_bend_kcal_mol, 2)
    }

if __name__ == '__main__':
    print("--- Module 02 (Polymer Mechanics): Worm-Like Chain Persistence Length ---")
    
    test_sequences = {
        "Canonical Random DNA (147 bp)": "ACGT" * 36 + "ACG",
        "Poly(dA:dT) Rigid Tract (147 bp)": "AAAA" * 36 + "AAA",
        "Flexible Hinge DNA (CA/TG repeats)": "CACACACACACACA" * 10 + "CACACACA"
    }
    
    for name, s in test_sequences.items():
        res = calculate_sequence_persistence_length(s)
        print(f"\nSequence: {name}")
        print(f"  • Persistence Length (L_p): {res['persistence_length_nm']} nm ({res['persistence_length_bp']} bp)")
        print(f"  • Mechanical Bending Energy (R=4.5 nm): {res['nucleosome_bending_energy_kcal_mol']} kcal/mol")
