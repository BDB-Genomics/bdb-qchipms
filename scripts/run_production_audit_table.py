"""
Production Audit Pipeline: Full Candidate Ranking & Genome-Wide Nearest Calls.
Grounded in Ensembl GRCh38 Release 113, dbSNP Build 156, and published literature ground truth.
Outputs raw TSV with all candidate evaluations and explicit filter statuses.
Strictly zero em-dashes.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from rank_candidates_engine import print_audit_report

if __name__ == '__main__':
    print_audit_report()
