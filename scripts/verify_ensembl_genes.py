import urllib.request
import json
import sys

genes = [
    ("MYC", "8"), ("PVT1", "8"), ("CASC8", "8"), ("CCAT2", "8"), ("POU5F1B", "8"), 
    ("CASC19", "8"), ("CASC11", "8"), ("LRATD2", "8"), ("CCDC26", "8"),
    ("FTO", "16"), ("IRX3", "16"), ("IRX5", "16"), ("RBL2", "16"),
    ("CDKN2A", "9"), ("CDKN2B", "9"), ("CDKN2B-AS1", "9"), ("MTAP", "9"), ("IFNA1", "9"),
    ("TCF7L2", "10"), ("VTI1A", "10"),
    ("UMOD", "16"), ("PDILT", "16")
]

print(f"{'Gene':<12} | {'Coordinates (GRCh38)':<32} | {'Strand':<6} | {'Gene_TSS':<10} | {'Canonical_Transcript':<18} | {'Canonical_TSS'}")
print("-" * 110)

for g, chrom in genes:
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{g}?expand=1"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req) as resp:
            data = json.loads(resp.read().decode())
            start = data.get("start")
            end = data.get("end")
            strand = "+" if data.get("strand") == 1 else "-"
            seq_region = data.get("seq_region_name")
            coord_str = f"chr{seq_region}:{start}-{end}"
            gene_tss = start if strand == "+" else end
            
            canonical_id = data.get("canonical_transcript", "")
            can_tss = None
            for tr in data.get("Transcript", []):
                if tr.get("id") == canonical_id or tr.get("is_canonical") == 1:
                    tr_strand = "+" if tr.get("strand") == 1 else "-"
                    can_tss = tr.get("start") if tr_strand == "+" else tr.get("end")
                    canonical_id = tr.get("id")
                    break
            if can_tss is None and data.get("Transcript"):
                tr = data.get("Transcript")[0]
                tr_strand = "+" if tr.get("strand") == 1 else "-"
                can_tss = tr.get("start") if tr_strand == "+" else tr.get("end")
                canonical_id = tr.get("id")
                
            print(f"{g:<12} | {coord_str:<32} | {strand:<6} | {gene_tss:<10} | {canonical_id:<18} | {can_tss}")
    except Exception as e:
        print(f"{g:<12} | Error: {e}")
