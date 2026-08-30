"""
Module 01: Groove Chemistry - Chemical Barcodes (Scale 1: 0.34 nm)
==================================================================
Computes the Seeman-Argos crystallographic chemical codes for major and minor grooves.
- Major Groove: Unique 4-character hydrogen bond donor/acceptor patterns (G-C: AADH, A-T: ADAM).
- Minor Groove: Degenerate 3-character patterns (A-T / T-A -> A-H-A; G-C / C-G -> A-D-A).
"""

from typing import Dict, List

# Chemical group codes:
# A = Hydrogen Bond Acceptor
# D = Hydrogen Bond Donor
# H = Non-polar Hydrogen / Methyl Group (Hydrophobic)
# M = Methyl Group (Thymine C5-methyl)

MAJOR_GROOVE_CODES = {
    'A': 'ADAM',  # N7(A), N6-H(D), O4(A), C5-CH3(M) - A:T
    'T': 'MADA',  # C5-CH3(M), O4(A), N6-H(D), N7(A) - T:A
    'G': 'AADH',  # N7(A), O6(A), N4-H(D), C5-H(H)   - G:C
    'C': 'HDAA',  # C5-H(H), N4-H(D), O6(A), N7(A)   - C:G
}

MINOR_GROOVE_CODES = {
    'A': 'A-H-A',  # N3(A), C2-H(H), O2(A) - A:T
    'T': 'A-H-A',  # O2(A), C2-H(H), N3(A) - T:A (Degenerate)
    'G': 'A-D-A',  # N3(A), N2-H(D), O2(A) - G:C
    'C': 'A-D-A',  # O2(A), N2-H(D), N3(A) - C:G (Degenerate)
}

def analyze_sequence_barcodes(sequence: str) -> List[Dict[str, str]]:
    """
    Computes major and minor groove chemical patterns for a DNA sequence.
    """
    seq = sequence.upper()
    results = []
    for i, base in enumerate(seq):
        if base in MAJOR_GROOVE_CODES:
            results.append({
                'pos': str(i + 1),
                'base': base,
                'major_code': MAJOR_GROOVE_CODES[base],
                'minor_code': MINOR_GROOVE_CODES[base],
                'minor_degeneracy': 'DEGENERATE (A-T == T-A)' if base in ('A', 'T') else 'DEGENERATE (G-C == C-G)'
            })
    return results

if __name__ == '__main__':
    test_seq = "GATCA"
    print(f"--- Chemical Barcode Analysis for: {test_seq} ---")
    barcodes = analyze_sequence_barcodes(test_seq)
    for b in barcodes:
        print(f"Pos {b['pos']} ({b['base']}) | Major Groove: {b['major_code']} (Distinct) | Minor Groove: {b['minor_code']} ({b['minor_degeneracy']})")
