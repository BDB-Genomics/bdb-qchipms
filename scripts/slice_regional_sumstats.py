#!/usr/bin/env python3
"""
BDB-Genomics: Single-pass regional slice extractor for the T2D atlas.

Streams the full DIAGRAM GWAS file once and emits every variant falling
inside any top-N locus window (from data/atlas/t2d_top60_loci_manifest.tsv)
to data/atlas/regional_slices/<locus_id>.tsv.gz.

Zero third-party dependencies (gzip, csv, os).
"""

import gzip
import csv
import os
import sys

MANIFEST = "data/atlas/t2d_top60_loci_manifest.tsv"
GWAS = "data/incoming/diagram_t2d_eur_gwas_catalog_raw.tsv.gz"
OUT_DIR = "data/atlas/regional_slices"

KEEP_COLS = [
    "hm_rsid", "hm_chrom", "hm_pos",
    "hm_effect_allele", "hm_other_allele",
    "beta", "standard_error", "p_value", "n",
]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    loci_by_chrom = {}
    with open(MANIFEST) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            chrom = row["chrom"].replace("chr", "")
            loci_by_chrom.setdefault(chrom, []).append(
                (int(row["window_start"]), int(row["window_end"]), row["locus_id"])
            )

    n_loci = sum(len(v) for v in loci_by_chrom.values())
    print(f"Slicing {n_loci} locus windows from {GWAS}")

    writers = {}

    def get_writer(locus_id):
        if locus_id not in writers:
            path = os.path.join(OUT_DIR, f"{locus_id}.tsv.gz")
            fh = gzip.open(path, "wt", newline="")
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["rsid", "chrom", "pos_grch38",
                        "effect_allele", "other_allele",
                        "beta", "se", "p_value", "n"])
            writers[locus_id] = (fh, w, 0)
        return writers[locus_id]

    with gzip.open(GWAS, "rt") as f:
        header = f.readline().strip().split("\t")
        idx = {name: i for i, name in enumerate(header)}
        cols = [idx[c] for c in KEEP_COLS]

        scanned = 0
        emitted = 0
        for line in f:
            scanned += 1
            if scanned % 1000000 == 0:
                print(f"  scanned {scanned:,} | emitted {emitted:,}")
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            chrom = parts[idx["hm_chrom"]].replace("chr", "")
            buckets = loci_by_chrom.get(chrom)
            if not buckets:
                continue
            try:
                pos = int(parts[idx["hm_pos"]])
            except ValueError:
                continue
            for start, end, locus_id in buckets:
                if start <= pos <= end:
                    fh, w, cnt = get_writer(locus_id)
                    vals = [parts[i] for i in cols]
                    # normalize chrom/pos types
                    vals[1] = str(chrom)
                    vals[2] = str(pos)
                    w.writerow(vals)
                    writers[locus_id] = (fh, w, cnt + 1)
                    emitted += 1

    total = 0
    for locus_id, (fh, _, cnt) in writers.items():
        fh.close()
        total += cnt
        print(f"  {locus_id}: {cnt:,} variants")

    print(f"DONE: {total:,} variants across {len(writers)} loci -> {OUT_DIR}")


if __name__ == "__main__":
    sys.exit(main())
