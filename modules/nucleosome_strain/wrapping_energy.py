"""
Subject: Molecular Biology & Chromatin Biophysics (Scale 3: 10 nm)
Core Concept: Nucleosome Wrapping Free Energy & Histone Octamer Barrier
Foundational References:
- Segal et al., Nature 442, 772–778 (2006) [Genomic code for nucleosome positioning]
- Richmond & Davey, Nature 423, 145–150 (2003) [Structure of DNA in nucleosome core particle at 1.9 A]
- Travers et al., PNAS 106, 501–506 (2009) [DNA mechanics and nucleosome organization]

Methodology & Physical Equations:
1. Canonical nucleosome wraps 147 bp of DNA in 1.65 superhelical turns around an octamer core.
2. Total wrapping free energy balance:
   Delta G_wrap = Delta G_bend(sequence) + Delta G_histone_contacts
   where Delta G_histone_contacts ~ -35.0 to -40.0 kcal/mol (14 electrostatic contact points with basic Arg/Lys patches).
3. If Delta G_bend exceeds ~35.0 kcal/mol (e.g. in rigid poly(dA:dT) tracts), Delta G_wrap becomes unfavorable (> 0),
   causing intrinsic nucleosome exclusion and maintaining open, accessible chromatin (ATAC-seq peak formation).
4. Periodic 10.5 bp rotational phasing (AA/TT vs GC steps) minimizes mechanical strain by aligning flexible
   steps facing inward toward the histone core.
"""

from typing import Dict
import math

from modules.polymer_mechanics.persistence_length import calculate_sequence_persistence_length

HISTONE_BINDING_AFFINITY_KCAL = -36.0 # kcal/mol standard stabilization from 14 contact points

def evaluate_nucleosome_occupancy_energy(sequence_147bp: str) -> Dict[str, float]:
    """
    Evaluates the net free energy (Delta G_wrap) for a 147 bp sequence wrapping around a histone octamer.
    Returns net free energy, predicted occupancy state, and ATAC-seq accessibility likelihood.
    """
    seq = sequence_147bp.upper()
    if len(seq) != 147:
        # Pad or slice to standard nucleosome core particle length
        if len(seq) < 147:
            seq = (seq * (147 // len(seq) + 1))[:147]
        else:
            seq = seq[:147]
            
    # Calculate bending strain energy from sequence persistence length
    poly_props = calculate_sequence_persistence_length(seq)
    delta_g_bend = poly_props['nucleosome_bending_energy_kcal_mol']
    
    # Check 10.5 bp helical phasing of AA/TT vs GC steps (rotational positioning bonus)
    phasing_bonus = 0.0
    for i in range(0, 147 - 10, 10):
        dinuc = seq[i:i+2]
        if dinuc in ('AA', 'TT', 'AT', 'TA'):
            phasing_bonus -= 0.6 # Favorable minor-groove-inward bend
        elif dinuc in ('GC', 'CG', 'GG', 'CC'):
            phasing_bonus -= 0.4 # Favorable major-groove-inward bend
            
    net_delta_g = delta_g_bend + HISTONE_BINDING_AFFINITY_KCAL + phasing_bonus
    
    # Equilibrium occupancy probability: P_occ = 1 / (1 + exp(Delta G / (k_B * T)))
    # At 37 C, k_B * T ~ 0.616 kcal/mol
    kbt = 0.616
    try:
        occupancy_prob = 1.0 / (1.0 + math.exp(net_delta_g / kbt))
    except OverflowError:
        occupancy_prob = 0.0 if net_delta_g > 0 else 1.0
        
    is_accessible_open_chromatin = bool(net_delta_g > 0.0 or occupancy_prob < 0.25)
    
    return {
        'sequence_length': len(seq),
        'persistence_length_nm': poly_props['persistence_length_nm'],
        'delta_G_bend_kcal_mol': round(delta_g_bend, 2),
        'histone_contact_affinity_kcal_mol': HISTONE_BINDING_AFFINITY_KCAL,
        'phasing_stabilization_kcal_mol': round(phasing_bonus, 2),
        'net_delta_G_wrap_kcal_mol': round(net_delta_g, 2),
        'predicted_nucleosome_occupancy_prob': round(occupancy_prob, 3),
        'chromatin_state': 'OPEN_ACCESSIBLE (Nucleosome Excluded)' if is_accessible_open_chromatin else 'CLOSED_OCCUPIED (Stable Nucleosome)'
    }

if __name__ == '__main__':
    print("--- Module 03 (Chromatin Biophysics): Nucleosome Wrapping Free Energy ---")
    
    seq_canonical = "ACGT" * 36 + "ACG"
    seq_rigid_tract = "AAAA" * 36 + "AAA"
    
    res_canon = evaluate_nucleosome_occupancy_energy(seq_canonical)
    res_rigid = evaluate_nucleosome_occupancy_energy(seq_rigid_tract)
    
    print(f"\n1. Flexible Canonical DNA (147 bp):")
    print(f"   • Bending Strain: {res_canon['delta_G_bend_kcal_mol']} kcal/mol")
    print(f"   • Net ΔG Wrap:    {res_canon['net_delta_G_wrap_kcal_mol']} kcal/mol")
    print(f"   • Occupancy:      {res_canon['predicted_nucleosome_occupancy_prob']} -> {res_canon['chromatin_state']}")
    
    print(f"\n2. Rigid Poly(dA:dT) Tract (147 bp):")
    print(f"   • Bending Strain: {res_rigid['delta_G_bend_kcal_mol']} kcal/mol")
    print(f"   • Net ΔG Wrap:    {res_rigid['net_delta_G_wrap_kcal_mol']} kcal/mol")
    print(f"   • Occupancy:      {res_rigid['predicted_nucleosome_occupancy_prob']} -> {res_rigid['chromatin_state']}")
