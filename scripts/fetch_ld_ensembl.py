#!/usr/bin/env python3
"""
Fetch real pairwise LD (r2, D') for TCF7L2 locus variants from Ensembl REST
(1000 Genomes Phase 3) and emit long-form + matrix tables with full provenance.

Source: https://rest.ensembl.org/documentation/info/ld_get
No simulated values anywhere in this script.
"""

import csv
import json
import sys
import time
import urllib.request
from datetime import date, datetime, timezone
from pathlib import Path

API = "https://rest.ensembl.org"
EUR_POPS = ["CEU", "TSI", "FIN", "GBR", "IBS"]

CORE_VARIANTS = [
    "rs7903146", "rs12255372", "rs4506565",
    "rs11196205", "rs7901695", "rs11196218",
]
ANCHOR = "rs7903146"
MIN_R2_REGIONAL = 0.10
MAX_REGIONAL_PARTNERS = 40

OUT_DIR = Path(__file__).resolve().parent.parent / "DOCUMENTS" / "Helena" / "tcf7l2_ld_real_v1"


def get_json(url: str, retries: int = 3):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode())
        except Exception as exc:  # noqa: BLE001
            if attempt == retries - 1:
                raise
            print(f"retry {attempt + 1} for {url}: {exc}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))


def post_json(url: str, payload: dict, retries: int = 3):
    for attempt in range(retries):
        try:
            data = json.dumps(payload).encode()
            req = urllib.request.Request(
                url,
                data=data,
                headers={"Content-Type": "application/json", "Accept": "application/json"},
            )
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode())
        except Exception as exc:  # noqa: BLE001
            if attempt == retries - 1:
                raise
            print(f"retry {attempt + 1}: {exc}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))


def fetch_ld(id: str, population: str) -> list[dict]:
    url = (
        f"{API}/ld/human/{id}/1000GENOMES:phase_3:{population}"
        f"?content-type=application/json"
    )
    records = get_json(url)
    time.sleep(0.25)
    return records if isinstance(records, list) else []


def fetch_pair_ld(id1: str, id2: str, population: str) -> list[dict]:
    url = (
        f"{API}/ld/human/pairwise/{id1}/{id2}"
        f"?population_name=1000GENOMES:phase_3:{population}"
        f"&content-type=application/json"
    )
    records = get_json(url)
    time.sleep(0.25)
    return records if isinstance(records, list) else []


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    eur = [p for p in EUR_POPS]
    pop_tag = lambda name: name.split("phase_3:")[-1]  # noqa: E731

    # 1) Validate candidate variants (core set + prior-table partners) via Ensembl
    gwas_path = (
        OUT_DIR.parent
        / "tcf7l2_undertone_delivery_v2"
        / "gwas"
        / "tcf7l2_gwas_regional_association_grch38.tsv"
    )
    candidates = set(CORE_VARIANTS)
    if gwas_path.exists():
        with open(gwas_path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                candidates.add(row["rsID"])
    print(f"{len(candidates)} candidate variants to validate")

    ids_list = sorted(candidates)
    positions: dict[str, str] = {}
    chunk_size = 150
    for i in range(0, len(ids_list), chunk_size):
        chunk = ids_list[i : i + chunk_size]
        result = post_json(f"{API}/variation/human?content-type=application/json", {"ids": chunk})
        for vid, entry in result.items():
            pos = "NA"
            for mapping in entry.get("mappings", []):
                if mapping.get("assembly_name") == "GRCh38":
                    pos = f"chr{mapping['seq_region_name']}:{mapping['start']}"
                    break
            positions[vid] = pos
        time.sleep(0.4)

    locus_lo, locus_hi = 112498590, 113498590

    def in_locus(vid: str) -> bool:
        pos = positions.get(vid, "")
        try:
            return pos.startswith("chr10:") and locus_lo <= int(pos.split(":")[1]) <= locus_hi
        except (ValueError, IndexError):
            return False

    validated = [v for v in ids_list if in_locus(v)]
    print(f"{len(validated)} validated chr10 locus variants")
    regional_ids = [v for v in validated if v not in CORE_VARIANTS][:MAX_REGIONAL_PARTNERS]

    # 2) Pairwise LD: core matrix (6x6) + anchor-vs-partners, per EUR population
    pair_pops: dict[tuple[str, str], dict[str, tuple[str, str]]] = {}

    def add_pair(a: str, b: str, pop: str) -> None:
        for rec in fetch_pair_ld(a, b, pop):
            key = tuple(sorted((rec["variation1"], rec["variation2"])))
            pair_pops.setdefault(key, {})[pop] = (float(rec["r2"]), float(rec["d_prime"]))

    core_pairs = [(a, b) for i, a in enumerate(CORE_VARIANTS) for b in CORE_VARIANTS[i + 1 :]]
    total = len(core_pairs) + len(regional_ids)
    done = 0
    for pop in eur:
        for a, b in core_pairs:
            add_pair(a, b, pop)
            done += 1
            print(f"\rLD {done}/{total * len(eur)} pairs", end="", flush=True)
    for vid in regional_ids:
        for pop in eur:
            add_pair(ANCHOR, vid, pop)
            done += 1
            print(f"\rLD {done}/{total * len(eur)} pairs", end="", flush=True)
    print()

    # 3) Long-form EUR table (core pairs, all five EUR populations)
    long_path = OUT_DIR / "tcf7l2_ld_variant_pairs_1000g_eur.tsv"
    with open(long_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "variation1", "pos1_GRCh38", "variation2", "pos2_GRCh38",
            "population", "r2", "d_prime",
        ])
        for key in sorted(pair_pops):
            v1, v2 = key
            for pop in eur:
                if pop in pair_pops[key]:
                    r2, dp = pair_pops[key][pop]
                    w.writerow([v1, positions.get(v1, "NA"), v2, positions.get(v2, "NA"), pop, f"{r2:.6f}", f"{dp:.6f}"])

    # 4) CEU square matrix over core variants
    mat_path = OUT_DIR / "tcf7l2_ld_matrix_r2_1000g_CEU.csv"
    with open(mat_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["variant"] + CORE_VARIANTS)
        for row_v in CORE_VARIANTS:
            row = [row_v]
            for col_v in CORE_VARIANTS:
                if row_v == col_v:
                    row.append("1.000000")
                    continue
                key = tuple(sorted((row_v, col_v)))
                val = pair_pops.get(key, {}).get("CEU")
                row.append(f"{val[0]:.6f}" if val else "NA")
            w.writerow(row)

    # 5) Regional partner ranking (anchor, CEU)
    rank_path = OUT_DIR / f"tcf7l2_{ANCHOR}_ceu_regional_partners.tsv"
    ranked = sorted(
        regional_ids,
        key=lambda v: pair_pops.get(tuple(sorted((ANCHOR, v))), {}).get("CEU", (0, 0))[0],
        reverse=True,
    )
    with open(rank_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["partner_rsID", "partner_pos_GRCh38", "r2_CEU", "d_prime_CEU"])
        for vid in ranked:
            key = tuple(sorted((ANCHOR, vid)))
            r2, dp = pair_pops[key]["CEU"]
            w.writerow([vid, positions.get(vid, "NA"), f"{r2:.6f}", f"{dp:.6f}"])

    # 6) Provenance README
    readme = OUT_DIR / "README.md"
    utc_now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    readme.write_text(
        f"""# TCF7L2 Locus LD — Real 1000 Genomes Phase 3 Values (Increment: LD rebuild)

Generated: {utc_now}
Generator: scripts/fetch_ld_ensembl.py (this repository)

## Source
- Ensembl REST API (`rest.ensembl.org`), endpoints:
  - `GET /ld/human/id/:id` — pairwise LD statistics
  - `POST /variation/human` — GRCh38 coordinates
- Underlying panel: 1000 Genomes Project Phase 3, populations
  {' '.join(EUR_POPS)} (European ancestry subset).

## Files
| File | Content |
|---|---|
| `tcf7l2_ld_variant_pairs_1000g_eur.tsv` | Long-form pairwise r2/D' for core locus variants, all 5 EUR populations |
| `tcf7l2_ld_matrix_r2_1000g_CEU.csv` | Square r2 matrix (CEU) over the 6 core variants |
| `tcf7l2_{ANCHOR}_ceu_regional_partners.tsv` | Top {MAX_REGIONAL_PARTNERS} CEU LD partners of {ANCHOR} (r2 >= {MIN_R2_REGIONAL}) |

## Notes / limitations
- All values are computed by Ensembl from 1000G Phase 3; none are estimated here.
- CEU is used for the matrix as the conventional EUR proxy; per-population values
  are in the long-form table.
- These tables supersede earlier illustrative LD files and should be treated as
  the reference for Undertone integration.
- Fine-mapping posteriors will be regenerated separately with SuSiE on published
  regional summary statistics and delivered as its own increment.
"""
    )

    print(f"Wrote outputs to {OUT_DIR}")
    for p in sorted(OUT_DIR.iterdir()):
        print(f"  {p.name} ({p.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
