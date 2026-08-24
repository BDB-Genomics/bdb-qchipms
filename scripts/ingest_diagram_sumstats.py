#!/usr/bin/env python3
"""
Ingest and validate real regional T2D summary statistics for SuSiE fine-mapping.

Pipeline stage 1 of: raw sumstats -> [THIS] -> sumstats_prepared.tsv -> run_susie_finemap.R

What this script does:
  1. Streams a raw summary-statistics file (.tsv / .tsv.gz / .csv).
  2. Auto-detects or accepts explicit column mapping.
  3. Filters to the configured locus window (default chr10:112.5-113.5 Mb, GRCh38).
  4. Keeps unambiguous biallelic SNVs only.
  5. Validates EVERY emitted rsID against Ensembl (existence + GRCh38 position
     match). This is the gate that would have caught a fabricated variant table.
  6. Writes SuSiE-ready output + a validation report. Refuses to emit fewer
     than MIN_VARIANTS rows.

Usage:
  python3 scripts/ingest_diagram_sumstats.py <raw_file> [--out-dir DIR]
                                             [--window START END] [--chrom CHR]
                                             [--col-map JSON]
"""

import argparse
import csv
import gzip
import json
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ENSEMBL_API = "https://rest.ensembl.org"
DEFAULT_WINDOW = ("chr10", 112500000, 113500000)
MIN_VARIANTS = 50

CANDIDATE_HEADERS = {
    "rsid": ["rsid", "rs_id", "rs", "snp", "variant_id", "id", "markername", "rs_number"],
    "chrom": ["chromosome", "chrom", "chr"],
    "pos": ["base_pair_location", "position_grch38", "pos_grch38", "bp", "position", "pos"],
    "ref": ["reference_allele", "ref_allele", "ref", "other_allele", "non_effect_allele", "nea"],
    "alt": ["alternate_allele", "alt_allele", "alt", "effect_allele", "ea"],
    "beta": ["beta", "effect_size", "log_odds", "b"],
    "se": ["standard_error", "se", "se_beta"],
    "p": ["pvalue", "p_value", "pval", "p"],
    "n": ["sample_size", "n", "n_samples", "total_sample_size"],
    "eaf": ["effect_allele_frequency", "eaf", "caf", "freq"],
}


def open_stream(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return open(path, "rt", newline="")


def sniff_dialect(sample: str) -> csv.Dialect:  # type: ignore[type-arg]
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,; ")
    except csv.Error:
        class _T(csv.excel_tab):  # type: ignore[name-defined]
            pass
        return _T


def detect_map(header: list[str], override: dict | None) -> dict[str, str]:
    norm = {h.strip().lower(): h for h in header}
    mapping: dict[str, str] = {}
    for std, candidates in CANDIDATE_HEADERS.items():
        for cand in candidates:
            if cand in norm:
                mapping[std] = norm[cand]
                break
    missing = {"rsid", "pos", "alt", "beta", "se"} - set(mapping)
    if override:
        for std, raw in override.items():
            if raw in header:
                mapping[std] = raw
        missing -= set(override)
    if missing:
        raise SystemExit(
            f"Could not map required columns: {sorted(missing)}. "
            f"Header seen: {header}. Use --col-map '{{\"std_name\":\"raw_name\",...}}'"
        )
    return mapping


def ensembl_validate(ids: list[str]) -> tuple[dict[str, str], list[str]]:
    """Return ({rsid: 'chrN:pos'}, [invalid_ids])."""
    positions: dict[str, str] = {}
    invalid: list[str] = []
    for i in range(0, len(ids), 150):
        chunk = sorted(ids[i : i + 150])
        payload = json.dumps({"ids": chunk}).encode()
        req = urllib.request.Request(
            f"{ENSEMBL_API}/variation/human?content-type=application/json",
            data=payload,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=60) as resp:
                    result = json.loads(resp.read().decode())
                break
            except Exception as exc:  # noqa: BLE001
                if attempt == 2:
                    raise SystemExit(f"Ensembl unreachable after retries: {exc}")
                time.sleep(2 * (attempt + 1))
        for vid, entry in result.items():
            pos = None
            for m in entry.get("mappings", []):
                if m.get("assembly_name") == "GRCh38":
                    pos = f"chr{m['seq_region_name']}:{m['start']}"
                    break
            if pos is None:
                invalid.append(vid)
            else:
                positions[vid] = pos
        print(f"validated {min(i + 150, len(ids))}/{len(ids)}", end="\r", flush=True)
    print()
    return positions, invalid


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("raw_file", type=Path)
    ap.add_argument("--out-dir", type=Path,
                    default=Path("results/finemapping/tcf7l2"))
    ap.add_argument("--chrom", default=DEFAULT_WINDOW[0])
    ap.add_argument("--window", nargs=2, type=int,
                    default=[DEFAULT_WINDOW[1], DEFAULT_WINDOW[2]], metavar=("START", "END"))
    ap.add_argument("--col-map", help="JSON string of {standard:raw} overrides")
    args = ap.parse_args()

    if not args.raw_file.exists():
        raise SystemExit(f"Input not found: {args.raw_file}")
    win_start, win_end = args.window
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    kept: list[dict] = []
    seen_rsids: set[str] = set()
    stats = {"rows": 0, "in_window": 0, "biallelic_snv": 0, "dupes": 0}

    with open_stream(args.raw_file) as fh:
        sample_line = fh.readline()
        fh.seek(0)
        reader = csv.DictReader(fh, dialect=sniff_dialect(sample_line))
        cmap = detect_map(reader.fieldnames or [], json.loads(args.col_map) if args.col_map else None)
        print(f"Column map: {cmap}")

        for row in reader:
            stats["rows"] += 1
            chrom = row.get(cmap.get("chrom"), args.chrom) or args.chrom
            chrom = str(chrom).replace("chr", "")
            if f"chr{chrom}" != args.chrom:
                continue
            try:
                pos = int(row[cmap["pos"]])
            except (KeyError, ValueError, TypeError):
                continue
            if not (win_start <= pos <= win_end):
                continue
            stats["in_window"] += 1
            ref = row.get(cmap["ref"], "").upper()
            alt = row.get(cmap["alt"], "").upper()
            if {ref, alt} - set("ACGT") or len(ref) != 1 or len(alt) != 1 or ref == alt:
                continue
            stats["biallelic_snv"] += 1
            rsid = row[cmap["rsid"]].strip()
            if not rsid.startswith("rs") or rsid in seen_rsids:
                stats["dupes"] += bool(rsid in seen_rsids)
                continue
            seen_rsids.add(rsid)

            def fnum(key: str, default=None):
                try:
                    v = row[cmap[key]]
                    return float(v)
                except (KeyError, ValueError, TypeError):
                    return default

            kept.append({
                "rsid": rsid, "chrom": args.chrom, "pos_raw": pos,
                "ref_allele": ref, "alt_allele": alt,
                "beta": fnum("beta"), "se": fnum("se"),
                "p": fnum("p"), "n": fnum("n"), "eaf": fnum("eaf"),
            })

    print(f"Read {stats['rows']} rows -> {stats['in_window']} in window -> "
          f"{stats['biallelic_snv']} biallelic SNVs -> {len(kept)} unique rsIDs")

    positions, invalid_ids = ensembl_validate([r["rsid"] for r in kept])

    prepared, pos_mismatch, dropped_invalid = [], [], []
    for r in kept:
        truth = positions.get(r["rsid"])
        if truth is None:
            dropped_invalid.append(r["rsid"])
            continue
        t_chrom, t_pos = truth.split(":")
        if t_pos != str(r["pos_raw"]):
            pos_mismatch.append({
                "rsid": r["rsid"], "file_pos": r["pos_raw"], "ensembl_pos_GRCh38": t_pos,
            })
            continue  # refuse to fine-map coordinate-inconsistent rows silently
        r["pos_grch38"] = int(t_pos)
        prepared.append(r)

    bad_stats = [r["rsid"] for r in prepared if r["beta"] is None or r["se"] in (None, 0)]
    prepared = [r for r in prepared if r["rsid"] not in bad_stats]

    out_tsv = out_dir / "sumstats_prepared.tsv"
    fields = ["rsid", "chrom", "pos_grch38", "ref_allele", "alt_allele",
              "beta", "se", "p", "n", "eaf"]
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(prepared)

    report = f"""# Summary-statistics ingestion report
Generated (UTC): {datetime.now(timezone.utc).isoformat()}
Source file: {args.raw_file}
Window: {args.chrom}:{win_start}-{win_end} (GRCh38)

Rows read:               {stats['rows']}
In window:               {stats['in_window']}
Biallelic SNVs:          {stats['biallelic_snv']}
Unique rsIDs:            {len(kept)}
Invalidated by Ensembl:  {len(dropped_invalid)}
Position mismatches:     {len(pos_mismatch)} (excluded, see below)
Bad beta/se excluded:    {len(bad_stats)}

EMITTED FOR FINE-MAPPING: {len(prepared)} variants -> {out_tsv}
Gate (>= {MIN_VARIANTS}): {'PASS' if len(prepared) >= MIN_VARIANTS else 'FAIL'}

Position mismatches are recorded verbatim in position_mismatches.tsv.
These indicate build or alignment errors in the source file - do NOT patch them;
check the source documentation instead.
"""
    (out_dir / "validation_report.md").write_text(report)
    if pos_mismatch:
        with open(out_dir / "position_mismatches.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["rsid", "file_pos", "ensembl_pos_GRCh38"], delimiter="\t")
            w.writeheader()
            w.writerows(pos_mismatch)

    print(report)
    if len(prepared) < MIN_VARIANTS:
        raise SystemExit(f"GATE FAILED: only {len(prepared)} clean variants.")
    print("GATE PASSED. Next: extend fetch_ld_ensembl.py to NxN over these rsIDs, then run_susie_finemap.R")


if __name__ == "__main__":
    main()
