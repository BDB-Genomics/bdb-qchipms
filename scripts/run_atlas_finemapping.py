#!/usr/bin/env python3
"""
BDB-Genomics: T2D atlas fine-mapping driver.

For each locus in the manifest:
  1. Harmonize regional sumstats against the EUR panel (.bim + AF):
     rsID join, allele alignment to panel A1/A2, sign flip on swapped
     alleles, drop of palindromic SNPs in the EAF gray zone.
  2. Compute signed LD (r) matrix with plink2 (--r square, MAF>=0.01).
  3. Convert .sqr -> CSV for run_susie_finemap.R and invoke SuSiE.
  4. Collect per-locus summary into data/atlas/finemapping_summary.tsv.

Zero third-party Python dependencies.
"""

import csv
import glob
import gzip
import os
import subprocess
import sys

MANIFEST = "data/atlas/t2d_top60_loci_manifest.tsv"
SLICES = "data/atlas/regional_slices"
REF = "data/atlas/ref_panel"
BED_GLOB = f"{REF}/chr_bed/chr*_eur"
PLINK2 = "tools/plink2"
RSCRIPT = ["tools/bin/micromamba", "run", "-n", "qchipms_r", "Rscript",
           "scripts/run_susie_finemap.R"]
OUT_ROOT = "results/atlas_finemapping"

PALINDROMES = {frozenset(("A", "T")), frozenset(("C", "G"))}
EAF_GRAY = (0.42, 0.58)


def load_bim():
    variants = {}
    for bim_path in sorted(glob.glob(f"{BED_GLOB}.bim")):
        with open(bim_path) as f:
            for line in f:
                chrom, rsid, _, pos, a1, a2 = line.split()[:6]
                variants[rsid] = (chrom, int(pos), a1.upper(), a2.upper())
    return variants


def load_af():
    af = {}
    for af_path in sorted(glob.glob(f"{BED_GLOB}.afreq")):
        with open(af_path) as f:
            hdr = f.readline().split()
            i_id, i_alt, i_af = hdr.index("ID"), hdr.index("ALT"), hdr.index("ALT_FREQS")
            for line in f:
                p = line.split()
                af[p[i_id]] = float(p[i_af])
    return af


def bed_for(chrom):
    path = f"{REF}/chr_bed/chr{chrom}_eur.bed"
    if not os.path.exists(path):
        return None
    return f"{REF}/chr_bed/chr{chrom}_eur"


def harmonize(locus, slice_path, bim, af):
    """Returns (rows, stats) with rows aligned to panel A1."""
    out, dropped = [], {"no_match": 0, "allele_mismatch": 0,
                        "palindromic_gray": 0, "bad_fields": 0}
    with gzip.open(slice_path, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rsid = row["rsid"]
            if rsid not in bim:
                dropped["no_match"] += 1
                continue
            chrom, pos, a1, a2 = bim[rsid]
            try:
                beta = float(row["beta"])
                se = float(row["se"])
                n = float(row["n"])
                eaf = float(row.get("eaf_ss") or "nan")
            except ValueError:
                dropped["bad_fields"] += 1
                continue
            ea = row["effect_allele"].upper()
            oa = row["other_allele"].upper()

            if {ea, oa} != {a1, a2}:
                dropped["allele_mismatch"] += 1
                continue

            palindromic = frozenset((ea, oa)) in PALINDROMES
            if palindromic and EAF_GRAY[0] <= eaf <= EAF_GRAY[1]:
                dropped["palindromic_gray"] += 1
                continue

            if ea == a2:  # effect stated on panel A2; flip to A1
                beta = -beta
                ea, oa = a1, a2
                eaf = 1.0 - eaf if eaf == eaf else eaf

            out.append({"rsid": rsid, "chrom": chrom, "pos_grch38": pos,
                        "ref_allele": a2, "alt_allele": a1,
                        "beta": beta, "se": se, "n": n})
    return out, dropped


def sqr_to_csv(locus_dir):
    ids_path = os.path.join(locus_dir, "locus_ld.sqr.ids")
    sqr_path = os.path.join(locus_dir, "locus_ld.sqr")
    with open(ids_path) as f:
        ids = [line.strip() for line in f if line.strip()]
    rows = []
    with open(sqr_path) as f:
        for line in f:
            rows.append(line.split())
    n = len(ids)
    assert all(len(r) == n for r in rows), "sqr not square"
    csv_path = os.path.join(locus_dir, "ld_matrix.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variant"] + ids)
        for rid, vals in zip(ids, rows):
            w.writerow([rid] + vals)
    return csv_path, n


def main():
    os.makedirs(OUT_ROOT, exist_ok=True)
    bim = load_bim()
    print(f"Panel: {len(bim):,} SNPs loaded from .bim")
    try:
        af = load_af()
        print(f"AF file: {len(af):,} entries")
    except FileNotFoundError:
        af = {}
        print("WARNING: no AF file; palindrome check limited to EAF field")

    summary = []
    with open(MANIFEST) as mf:
        loci = list(csv.DictReader(mf, delimiter="\t"))

    only = set(sys.argv[1:])  # optional locus_ids to rerun
    for locus in loci:
        lid = locus["locus_id"]
        if only and lid not in only:
            continue
        slice_path = os.path.join(SLICES, f"{lid}.tsv.gz")
        if not os.path.exists(slice_path):
            print(f"{lid}: MISSING SLICE, skipped")
            continue

        locus_dir = os.path.join(OUT_ROOT, lid)
        os.makedirs(locus_dir, exist_ok=True)

        rows, dropped = harmonize(locus, slice_path, bim, af)
        if len(rows) < 50:
            print(f"{lid}: only {len(rows)} harmonized variants - skipping SuSiE")
            summary.append({"locus_id": lid, "lead_rsid": locus["lead_rsid"],
                            "n_harmonized": len(rows), "n_credible_sets": "",
                            "lead_cs_variant": "", "status": "too_few_variants"})
            continue

        ss_path = os.path.join(locus_dir, "sumstats.tsv")
        cols = ["rsid", "chrom", "pos_grch38", "ref_allele", "alt_allele",
                "beta", "se", "n"]
        with open(ss_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

        bed_prefix = bed_for(locus["chrom"].replace("chr", ""))
        if bed_prefix is None:
            print(f"{lid}: no reference BED for {locus['chrom']} - skipped")
            summary.append({"locus_id": lid, "lead_rsid": locus["lead_rsid"],
                            "n_harmonized": len(rows), "n_credible_sets": "",
                            "lead_cs_variant": "", "status": "no_reference"})
            continue

        subprocess.run([
            PLINK2, "--bfile", bed_prefix,
            "--from-bp", locus["window_start"],
            "--to-bp", locus["window_end"],
            "--maf", "0.01",
            "--r", "square",
            "--memory", "2000", "--threads", "8",
            "--out", os.path.join(locus_dir, "locus_ld"),
        ], check=True, capture_output=True)

        csv_path, n_ld = sqr_to_csv(locus_dir)
        gwas_n = int(round(sum(r["n"] for r in rows) / len(rows)))

        r = subprocess.run(RSCRIPT + [ss_path, csv_path, str(gwas_n), locus_dir],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"{lid}: SUSIE FAILED\n{r.stderr[-500:]}")
            summary.append({"locus_id": lid, "lead_rsid": locus["lead_rsid"],
                            "n_harmonized": len(rows), "n_credible_sets": "",
                            "lead_cs_variant": "", "status": "susie_error"})
            continue

        cs_path = os.path.join(locus_dir, "credible_sets_95.tsv")
        n_sets, lead_variants = 0, []
        if os.path.exists(cs_path):
            with open(cs_path) as f:
                cs_rows = list(csv.DictReader(f, delimiter="\t"))
                n_sets = len(cs_rows)
                lead_variants = [c["lead_rsID"] for c in cs_rows]

        pip_path = os.path.join(locus_dir, "susie_pips.tsv")
        top_pip, top_rsid = 0.0, ""
        if os.path.exists(pip_path):
            with open(pip_path) as f:
                first = next(iter(csv.DictReader(f, delimiter="\t")), None)
                if first:
                    top_pip, top_rsid = float(first["pip"]), first["rsid"]

        print(f"{lid}: harm={len(rows)} ld={n_ld} cs={n_sets} "
              f"topPIP={top_rsid}:{top_pip:.3f} drops={dropped}")
        summary.append({
            "locus_id": lid, "lead_rsid": locus["lead_rsid"],
            "n_harmonized": len(rows),
            "n_matched_ld": n_ld - len(rows) + len(rows),
            "n_credible_sets": n_sets,
            "lead_cs_variant": ";".join(lead_variants[:5]),
            "top_pip_variant": top_rsid, "top_pip": round(top_pip, 4),
            "status": "ok",
        })

    if summary:
        out_path = os.path.join("data/atlas/finemapping_summary.tsv")
        existing = {}
        if os.path.exists(out_path):
            with open(out_path) as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    existing[row["locus_id"]] = row
        for row in summary:
            existing[row["locus_id"]] = row
        merged = sorted(existing.values(), key=lambda r: r["locus_id"])
        with open(out_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(merged[0].keys()),
                               delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(merged)
        print(f"\nSummary -> {out_path} ({len(merged)} loci)")


if __name__ == "__main__":
    main()
