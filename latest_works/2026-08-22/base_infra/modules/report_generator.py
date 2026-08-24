#!/usr/bin/env python3
"""
BDB-Genomics: Base Master Report & UpSet Evidence Aggregator
------------------------------------------------------------
Combines all 7 pillars into a master tabular target dossier and generates
multi-tier UpSet priority classifications.
"""

import os
import csv
from typing import Dict, List, Any

def generate_master_dossier_tsv(records: List[Dict[str, Any]], out_path: str):
    """
    Writes complete master multi-omic evidence dossier to TSV.
    """
    if not records:
        return
    fieldnames = list(records[0].keys())
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in records:
            writer.writerow(r)
    print(f"\n[Master Dossier] Successfully wrote {len(records)} multi-omic entries to:\n  {out_path}")
