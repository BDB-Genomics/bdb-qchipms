#!/usr/bin/env python3
"""
Build a local NxN LD (r2) matrix for fine-mapping from 1000 Genomes genotypes.

Reference panel: 1000G 2504 high-coverage (2022 release, GRCh38) - build matches
the ingested GWAS Catalog summary statistics exactly. EUR-5 populations only
(CEU, TSI, FIN, GBR, IBS), matching the European-ancestry DIAMANTE meta-analysis.

Inputs : results/finemapping/tcf7l2/sumstats_prepared.tsv (from ingest stage)
Output : results/finemapping/tcf7l2/ld_matrix_1000g_hc_eur.csv
         results/finemapping/tcf7l2/ld_provenance.md
         results/finemapping/tcf7l2/ld_missing_variants.tsv
"""

import csv
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pysam

VCF_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
    "1kGP_high_coverage_Illumina.chr10.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
)
PANEL_URLS = [
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt",
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel",
]
EUR_POPS = {"CEU", "TSI", "FIN", "GBR", "IBS"}
EUR_SUPERPOP = "EUR"
MAX_MISSING_FRAC = 0.02

BASE = Path("results/finemapping/tcf7l2")
SUMSTATS = BASE / "sumstats_prepared.tsv"
OUT_CSV = BASE / "ld_matrix_1000g_hc_eur.csv"
OUT_PROV = BASE / "ld_provenance.md"
OUT_MISS = BASE / "ld_missing_variants.tsv"


def load_eur_samples() -> set[str]:
    last_err = None
    for url in PANEL_URLS:
        try:
            text = urllib.request.urlopen(url, timeout=60).read().decode()
            samples = set()
            mode = None
            for line in text.splitlines():
                t = line.split()
                if not t:
                    continue
                if t[0] == "FamilyID":  # ped format: FID IID PAT MAT SEX POP SUPERPOP
                    mode = "ped"
                    continue
                if t[0] == "sample":  # panel format: sample pop super_pop gender
                    mode = "panel"
                    continue
                if mode == "ped" and len(t) >= 7 and (t[6] == EUR_SUPERPOP or t[5] in EUR_POPS):
                    samples.add(t[1])
                elif mode == "panel" and len(t) >= 3 and (t[2] == EUR_SUPERPOP or t[1] in EUR_POPS):
                    samples.add(t[0])
            if samples:
                print(f"EUR sample map: {len(samples)} individuals from {url.rsplit('/', 1)[-1]}")
                return samples
        except Exception as exc:  # noqa: BLE001
            last_err = exc
    raise SystemExit(f"Could not load sample population panel: {last_err}")


def main() -> None:
    with open(SUMSTATS) as fh:
        variants = list(csv.DictReader(fh, delimiter="\t"))
    wanted = {}
    for v in variants:
        wanted[(int(v["pos_grch38"]), v["ref_allele"], v["alt_allele"])] = v["rsid"]
    print(f"{len(wanted)} target variants")

    eur = load_eur_samples()

    vcf = pysam.VariantFile(VCF_URL)
    all_samples = list(vcf.header.samples)
    eur_idx = [i for i, s in enumerate(all_samples) if s in eur]
    n_eur = len(eur_idx)
    print(f"Panel samples: {len(all_samples)} total, {n_eur} EUR-5")

    chrom = "chr10" if "chr10" in vcf.header.contigs else "10"
    positions = [k[0] for k in wanted]
    lo, hi = min(positions) - 1, max(positions)
    print(f"Fetching {chrom}:{lo}-{hi} (remote tabix slice)...")

    dosages: dict[tuple[int, str, str], np.ndarray] = {}
    n_records = 0
    for rec in vcf.fetch(chrom, lo, hi):
        n_records += 1
        if rec.alts is None or len(rec.alts) != 1:
            continue
        key = (rec.pos, rec.ref, rec.alts[0].upper())
        if key not in wanted:
            continue
        g = np.full(n_eur, np.nan)
        samples = rec.samples
        for j, i in enumerate(eur_idx):
            call = samples[all_samples[i]]
            a = call.allele_indices
            if a is None or a[0] is None or a[1] is None:
                continue
            g[j] = a[0] + a[1]
        dosages[key] = g
    vcf.close()
    print(f"Scanned {n_records} region records; matched {len(dosages)}/{len(wanted)} target variants")

    found_keys = set(dosages)
    # drop high-missing variants
    kept_keys = [k for k in found_keys if np.isnan(dosages[k]).mean() <= MAX_MISSING_FRAC]
    missing = [wanted[k] for k in wanted if k not in kept_keys]
    X = np.vstack([np.nan_to_num(dosages[k]) for k in kept_keys])
    rsids = [wanted[k] for k in kept_keys]
    print(f"LD computed over {len(rsids)} variants x {n_eur} EUR samples "
          f"({len(missing)} targets absent/high-missing)")

    r = np.corrcoef(X)
    r2 = np.round(np.clip(r, -1, 1) ** 2, 6)

    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["variant"] + rsids)
        for i, rid in enumerate(rsids):
            w.writerow([rid] + [f"{v:.6f}" for v in r2[i]])

    with open(OUT_MISS, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["rsid", "reason"])
        for rid in missing:
            w.writerow([rid, "not_found_or_high_missing_in_1000G_HC_panel"])

    OUT_PROV.write_text(f"""# LD reference provenance

Generated (UTC): {datetime.now(timezone.utc).isoformat()}
Generator: scripts/fetch_genotypes_1000g.py

## Panel
- 1000 Genomes Project 2504 high-coverage SNV/INDEL panel (2022 release)
- VCF: {VCF_URL}
- Build: GRCh38 (matches ingested summary statistics build exactly)
- Populations: {' '.join(sorted(EUR_POPS))} (EUR-5), n={n_eur} individuals
- Sample population map: {PANEL_URLS[0].rsplit('/', 1)[-1]}

## Method
- Region fetched remotely via tabix index (bgzf range requests); no full-file download.
- Variant matching by (POS, REF, ALT) against sumstats_prepared.tsv; VCF IDs ignored.
- Dosage = alternate-allele count (0/1/2); phased GT used, missing calls as NaN.
- Variants with >{int(MAX_MISSING_FRAC * 100)}% missing calls dropped (see ld_missing_variants.tsv).
- r2 = squared Pearson correlation of dosages across EUR samples (numpy.corrcoef).

## Output
- ld_matrix_1000g_hc_eur.csv: {len(rsids)} x {len(rsids)} r2 matrix
- ld_missing_variants.tsv: {len(missing)} target variants absent from panel
""")
    print(f"Wrote {OUT_CSV}, {OUT_PROV}")


if __name__ == "__main__":
    main()
