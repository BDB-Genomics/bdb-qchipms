"""
Subject: Biophysics & Physical Chemistry
Core Concept: Poisson-Boltzmann Electrostatics & Minor Groove Electric Wells
Foundational References:
- Rohs et al., Nature 461, 1248–1253 (2009) [PMID: 19865164]
- Honig & Nicholls, Science 268, 1144–1149 (1995) [PMID: 7761829]

Methodology & Physical Equations:
1. Linearized Poisson-Boltzmann Equation in cylindrical dielectric medium:
   nabla^2 Phi(r) = kappa^2 * Phi(r)
   where kappa = sqrt(8 * pi * e^2 * I / (epsilon * k_B * T)) is the Debye-Huckel screening parameter.
2. In narrow minor grooves (width w ~ 5.8 A), the two phosphate rails (charge -e per nucleotide)
   approach each other, focusing the electrostatic potential well:
   Phi_well proportional to - (2 * e) / (epsilon_eff * w)
3. Arginine finger insertion binding free energy:
   Delta G_elec = q_Arg * Phi_well = (+1 e) * Phi_well (exergonic, stabilizes binding by -3 to -6 kcal/mol).
"""

import math
from typing import Dict

# Physical Constants at 310.15 K (37 C)
EPSILON_WATER = 78.5      # Bulk water dielectric
EPSILON_DNA_CORE = 4.0    # Hydrophobic core dielectric
EPSILON_GROOVE = 20.0     # Effective dielectric in structured water groove (Honig et al.)
E_CHARGE = 1.602176634e-19 # Coulombs
KB = 1.380649e-23         # J/K
TEMP = 310.15             # Kelvin
IONIC_STRENGTH = 0.15     # 150 mM physiological salt (Molar)
KBT_EV = (KB * TEMP) / E_CHARGE # ~0.0267 eV = 1 k_B*T/e

def calculate_minor_groove_electrostatics(groove_width_angstrom: float) -> Dict[str, float]:
    """
    Calculates the electrostatic potential well depth (in k_B*T/e and kcal/mol)
    and Arginine insertion binding free energy as a function of minor groove width.
    """
    w_meters = groove_width_angstrom * 1e-10
    
    # Debye screening length lambda_D = 1 / kappa
    # kappa ~ 0.329 * sqrt(I) nm^-1 for water at 37 C
    kappa = 0.329e9 * math.sqrt(IONIC_STRENGTH) # m^-1
    debye_length_angstrom = (1.0 / kappa) * 1e10 # ~7.95 A
    
    # Electrostatic potential at center of groove from two opposing phosphate charges (-e each)
    # Screened Coulomb potential: Phi = 2 * (-e / (4 * pi * epsilon_0 * epsilon_r * (w/2))) * exp(-kappa * w/2)
    eps_0 = 8.8541878128e-12
    r = w_meters / 2.0
    phi_volts = 2.0 * (-E_CHARGE / (4.0 * math.pi * eps_0 * EPSILON_GROOVE * r)) * math.exp(-kappa * r)
    
    # Potential in thermal energy units (k_B*T/e)
    phi_kbt = phi_volts / KBT_EV
    
    # Arginine (+1 e) insertion free energy in kcal/mol (1 eV = 23.06 kcal/mol)
    delta_g_kcal_mol = (phi_volts * 1.0) * 23.0605
    
    return {
        'groove_width_A': round(groove_width_angstrom, 2),
        'debye_screening_length_A': round(debye_length_angstrom, 2),
        'potential_volts': round(phi_volts, 3),
        'potential_kbt_e': round(phi_kbt, 2),
        'arginine_insertion_deltaG_kcal_mol': round(delta_g_kcal_mol, 2)
    }

if __name__ == '__main__':
    print("--- Module 01 (Biophysics): Poisson-Boltzmann Minor Groove Electrostatics ---")
    print("Evaluating canonical A-tract narrow groove (5.6 A) vs G/C wide groove (7.2 A):\n")
    
    for w in [5.2, 5.6, 6.0, 6.8, 7.2, 8.0]:
        res = calculate_minor_groove_electrostatics(w)
        print(f"Width: {res['groove_width_A']} Å | Phi: {res['potential_kbt_e']:>6} k_B*T/e ({res['potential_volts']:>6} V) | Arg Insertion ΔG: {res['arginine_insertion_deltaG_kcal_mol']:>6} kcal/mol")
