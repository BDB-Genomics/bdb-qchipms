#!/usr/bin/env bash
# Fixup: rebuild a chromosome's EUR reference using the annotated _rs pvar
# (fallback when the resources-page _rs_noannot link is corrupt).
# Usage: fixup_chr_pvar.sh <chrom> <pvar_share_id> <rlkey>
set -euo pipefail
cd "$(dirname "$0")/.."

N="$1"; PV_ID="$2"; PV_KEY="$3"
REF=data/atlas/ref_panel
BEDDIR=$REF/chr_bed
TOOLS=tools/plink2

OUT="$BEDDIR/chr${N}_eur"
if [ -f "$OUT.afreq" ]; then echo "chr$N: already built"; exit 0; fi

if [ ! -f "$REF/chr${N}_hg38.pgen" ]; then
  curl -fSsL --retry 15 --retry-all-errors -C - --speed-time 30 --speed-limit 102400 \
    -o "$REF/chr${N}_hg38.pgen.zst" "https://www.dropbox.com/s/PLACEHOLDER/chr${N}_hg38.pgen.zst?dl=1"
  zstd -d -f "$REF/chr${N}_hg38.pgen.zst" -o "$REF/chr${N}_hg38.pgen"
fi

curl -fSsL --retry 15 --retry-all-errors -C - --speed-time 30 --speed-limit 102400 \
  -o "$REF/chr${N}_hg38_rs.pvar.zst" \
  "https://www.dropbox.com/scl/fi/${PV_ID}/chr${N}_hg38_rs.pvar.zst?rlkey=${PV_KEY}&dl=1"
zstd -d -f "$REF/chr${N}_hg38_rs.pvar.zst" -o "$REF/chr${N}_hg38.pvar"
cp "$REF/hg38_corrected.psam" "$REF/chr${N}_hg38.psam"

"$TOOLS" --pfile "$REF/chr${N}_hg38" \
  --keep "$REF/eur_samples.txt" \
  --snps-only just-acgt --max-alleles 2 \
  --memory 2500 --threads 4 \
  --make-bed --out "$OUT"

"$TOOLS" --bfile "$OUT" --memory 2000 --threads 4 --freq --out "$OUT"

rm -f "$REF/chr${N}_hg38.pgen" "$REF/chr${N}_hg38.pvar" \
      "$REF/chr${N}_hg38.psam" "$REF/chr${N}_hg38.pgen.zst" \
      "$REF/chr${N}_hg38_rs.pvar.zst"
echo "chr$N FIXUP DONE"
