#!/usr/bin/env bash
# Per-chromosome EUR reference builder (RAM-safe: processes one chr at a time).
# Source: PLINK2 resources mirror, 1000G 2022 high-coverage callset (hg38, 3202 samples).
# Output: data/atlas/ref_panel/chr_bed/chrN_eur.{bed,bim,fam} + .afreq for chr 1-19 (loci-bearing).
set -uo pipefail
cd "$(dirname "$0")/.."

REF=data/atlas/ref_panel
BEDDIR=$REF/chr_bed
TOOLS=tools/plink2
mkdir -p "$BEDDIR"

declare -A PGEN PVAR
PGEN[1]=k5tul2pe5wytjlz;  PVAR[1]=0y3nyl1922gb3ni4vrgl1/suzj2r41pjguaar2ma73z6lnw
PGEN[2]=nsgqqahbfk4u9r8;  PVAR[2]=yac4wj78p9bar9c/qyp9iqkbdge7m5ug12mhbxbm8
PGEN[3]=q400ppmtps05147;  PVAR[3]=0ykeue1c72t5l46mftynv/0d9g94whig2z0f40qq6trzo14
PGEN[4]=63acikvz57nnc1z;  PVAR[4]=bkv4fyrmp6kesno8c52cy/623mkyra19qppr6aix1de1eug
PGEN[5]=wg01ye3cn06mm8e;  PVAR[5]=3trhniyvugb1gk8law4ap/9k6te666d02mbefnekvktlq35
PGEN[6]=7l4ozeaze53mlny;  PVAR[6]=10lulcz9cjw955rrp5lc4/cc3x3njk82scatsr2725gotbx
PGEN[7]=qiqiq91gjnlejj5;  PVAR[7]=pyqm2xhp3abd5j32904gtstyo/pyqm2xhp3abd5j32904gtstyo
PGEN[8]=dvtewdfrsa63orj;  PVAR[8]=tc0nxum3umqy46z8pk5wf/u986xnhq2ybbug9hqh7b4lpcd
PGEN[9]=321ercemdimyh8j;  PVAR[9]=ko7tm3zqmoy46q9rty8tk/2zn6tvkvb1d5jgiheuh9lojqf
PGEN[10]=767hbzu2qycx46q; PVAR[10]=jan6qqw17qpd3k3hcp4qk/gq11tsvqsla37ycqsn90ud1j5
PGEN[11]=jiykipnk2qleubn; PVAR[11]=d1iljtkfa7n51wj8gokmh/5wp6qtkrlgpabj09mbjte0oiu
PGEN[12]=gqbwlavrm2rfm3m; PVAR[12]=hmxmf32espmz3aq96xbz2/4ny14x48ngkb2ii07orvvorkg
PGEN[13]=f5xfqfc8g5lh70s; PVAR[13]=5vkbogy0jth6a5i9h1mvs/d4e9ow7o05waxac8uzlyeob9z
PGEN[15]=b3zp9p0bb21hrht; PVAR[15]=8g03djru6fzb7newpykef/eefk8ufbatx7nqwltnhy1cj13
PGEN[16]=bwyzc8qo7jtk8uj; PVAR[16]=324av0o5wxpppmc24z0uz/nvbbp1nqs6g0rj1scoqox7nxq
PGEN[17]=wnjclot6ksza8ru; PVAR[17]=citw5fjlszpkcnus2ync2/8vlt05v0i32boaii1bc299pvb
PGEN[18]=8g92e4f61yqyvas; PVAR[18]=kh6a3r4952h3952by04em/mz5w5i4hmh67i80jqr0uziblt
PGEN[19]=xfyp8uymflqts78; PVAR[19]=35ojjv4xc4pize9whtftp/f231kmks01mzz7lqisdc1uk6c

CURL="curl -fSsL --retry 15 --retry-all-errors -C - --speed-time 30 --speed-limit 102400"

for N in 1 2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19; do
  OUT="$BEDDIR/chr${N}_eur"
  if [ -f "$OUT.afreq" ]; then echo "chr$N: already built, skipping"; continue; fi
  echo "=== chr$N ==="

  pv_id="${PVAR[$N]%%/*}"; pv_key="${PVAR[$N]##*/}"
  [ -f "$REF/chr${N}_hg38.pgen" ] || { $CURL -o "$REF/chr${N}_hg38.pgen.zst" "https://www.dropbox.com/s/${PGEN[$N]}/chr${N}_hg38.pgen.zst?dl=1" || { echo "chr$N pgen DOWNLOAD FAILED"; continue; }; }
  [ -f "$REF/chr${N}_hg38.pvar" ] || { $CURL -o "$REF/chr${N}_hg38.pvar.zst" "https://www.dropbox.com/scl/fi/${pv_id}/chr${N}_hg38_rs_noannot.pvar.zst?rlkey=${pv_key}&dl=1" || { echo "chr$N pvar DOWNLOAD FAILED"; continue; }; }

  [ -f "$REF/chr${N}_hg38.pgen" ] || zstd -d -f "$REF/chr${N}_hg38.pgen.zst" -o "$REF/chr${N}_hg38.pgen" || { echo "chr$N pgen DECOMP FAILED"; continue; }
  [ -f "$REF/chr${N}_hg38.pvar" ] || zstd -d -f "$REF/chr${N}_hg38.pvar.zst" -o "$REF/chr${N}_hg38.pvar" || { echo "chr$N pvar DECOMP FAILED"; continue; }
  [ -f "$REF/chr${N}_hg38.psam" ] || cp "$REF/hg38_corrected.psam" "$REF/chr${N}_hg38.psam"

  "$TOOLS" --pfile "$REF/chr${N}_hg38" \
    --keep "$REF/eur_samples.txt" \
    --snps-only just-acgt --max-alleles 2 \
    --memory 2500 --threads 4 \
    --make-bed --out "$OUT" || { echo "chr$N BED BUILD FAILED"; continue; }

  "$TOOLS" --bfile "$OUT" --memory 2000 --threads 4 --freq --out "$OUT" || { echo "chr$N AF FAILED"; continue; }

  rm -f "$REF/chr${N}_hg38.pgen" "$REF/chr${N}_hg38.pvar" "$REF/chr${N}_hg38.psam"
  echo "chr$N: DONE ($(stat -c%s "$OUT.bed") bytes bed)"
done

echo "PERCHR_BUILD_COMPLETE"
ls -la "$BEDDIR" | awk '{print $5, $9}'
