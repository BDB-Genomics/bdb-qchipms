"""
Subject: Rotational Biophysics (Scale 2: 10 nm)
Core Concept: Helical Phasing & 10.5 bp Periodicity under Indel Mutations
Foundational References:
- Wang & Giaever, Science 239, 1374–1377 (1988) [Action at a distance along a DNA: role of DNA supercoiling and phasing]
- Ptashne, A Genetic Switch: Phage Lambda and Higher Organisms (3rd ed.) [Stereospecific cooperative binding]

Methodology:
1. B-DNA has a helical repeat of 10.5 base pairs per 360-degree turn (34.3 degrees per base pair).
2. Insertion or deletion of Delta L base pairs rotates the downstream helical face by:
   theta = (Delta L * 360 / 10.5) mod 360 degrees
3. If two cooperative transcription factors are separated by d bp:
   - Constructive phase (theta ~ 0 or 360): Proteins are on the same face -> Cooperative interaction.
   - Destructive phase (theta ~ 180): Proteins are on opposite faces -> Synergistic binding completely broken.
"""

from typing import Dict

HELICAL_REPEAT_BP = 10.5

def calculate_helical_phasing_rotation(indel_size_bp: int) -> Dict[str, any]:
    """
    Computes angular rotation of the DNA double helix caused by an insertion or deletion.
    """
    rotation_deg = (indel_size_bp * 360.0 / HELICAL_REPEAT_BP) % 360.0
    
    # Measure angular deviation from optimal same-face alignment (0 / 360 deg)
    phase_deviation_deg = min(rotation_deg, 360.0 - rotation_deg)
    
    if phase_deviation_deg < 45.0:
        face_status = "SAME_FACE (Cooperative Interaction Preserved)"
        cooperativity_loss_pct = (phase_deviation_deg / 45.0) * 15.0
    elif phase_deviation_deg > 135.0:
        face_status = "OPPOSITE_FACE (180° Inverted - Cooperativity Destroyed)"
        cooperativity_loss_pct = 90.0 + (180.0 - phase_deviation_deg) * 0.2
    else:
        face_status = "PARTIAL_ROTATION (Sterically Displaced)"
        cooperativity_loss_pct = 40.0 + (phase_deviation_deg - 45.0) * 0.5
        
    return {
        'indel_size_bp': indel_size_bp,
        'helical_rotation_degrees': round(rotation_deg, 1),
        'phase_deviation_degrees': round(phase_deviation_deg, 1),
        'helical_face_alignment': face_status,
        'cooperativity_loss_percentage': round(min(100.0, cooperativity_loss_pct), 1)
    }

if __name__ == '__main__':
    print("--- Module: Helical Phasing & 10.5 bp Periodicity ---")
    for indel in [0, 5, 10, 11, 15]:
        res = calculate_helical_phasing_rotation(indel)
        print(f"Indel: {res['indel_size_bp']:>2} bp -> Rotation: {res['helical_rotation_degrees']:>5}° | Status: {res['helical_face_alignment']} | Loss: {res['cooperativity_loss_percentage']}%")
