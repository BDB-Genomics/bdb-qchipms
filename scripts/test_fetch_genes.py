import urllib.request
import json
import math
from pathlib import Path

# Function to fetch all annotated genes in a genomic window from Ensembl REST API
def fetch_ensembl_genes_in_window(chrom, start, end):
    # Normalize chrom
    c = str(chrom).replace('chr', '')
    url = f"https://rest.ensembl.org/overlap/region/human/{c}:{start}-{end}?feature=gene"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req) as resp:
            data = json.loads(resp.read().decode())
            genes = []
            for item in data:
                name = item.get("external_name") or item.get("id")
                strand = "+" if item.get("strand") == 1 else "-"
                g_start = item.get("start")
                g_end = item.get("end")
                tss = g_start if strand == "+" else g_end
                genes.append({
                    'id': item.get("id"),
                    'name': name,
                    'biotype': item.get("biotype"),
                    'chrom': f"chr{c}",
                    'start': g_start,
                    'end': g_end,
                    'strand': strand,
                    'tss': tss
                })
            return genes
    except Exception as e:
        print(f"Error fetching {chrom}:{start}-{end}: {e}")
        return []

if __name__ == '__main__':
    # Test 8q24.21 window
    genes = fetch_ensembl_genes_in_window("8", 127200000, 127850000)
    print(f"Retrieved {len(genes)} Ensembl genes in 8q24.21 window:")
    for g in sorted(genes, key=lambda x: x['tss']):
        print(f"{g['name']:<15} | {g['biotype']:<20} | {g['chrom']}:{g['start']}-{g['end']} ({g['strand']}) | TSS: {g['tss']}")
