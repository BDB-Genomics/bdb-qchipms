import csv
import random
import os

os.makedirs("/home/himanshu/bdb-qchipms/data", exist_ok=True)
random.seed(42)

real_hits = ["TERF2", "TERF1", "POT1", "TINF2", "ACD", "TERF2IP", "ZBTB48", "HMBOX1", "NR2C2", "NR2F2", "ZNF827"]
background_proteins = [f"PROT_{i:03d}" for i in range(1, 608)]
all_proteins = real_hits + background_proteins

out_path = "/home/himanshu/bdb-qchipms/data/protein_intensities.tsv"

with open(out_path, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["Protein_ID", "Gene_Name", "TERF2_rep1", "TERF2_rep2", "IgG_rep1", "IgG_rep2"])
    
    for prot in all_proteins:
        if prot in real_hits:
            terf2_1 = round(2**(random.gauss(28.5, 0.4)), 2)
            terf2_2 = round(2**(random.gauss(28.3, 0.4)), 2)
            igg_1 = round(2**(random.gauss(19.0, 0.5)), 2)
            igg_2 = round(2**(random.gauss(19.2, 0.5)), 2)
        else:
            base = random.gauss(20.0, 1.0)
            terf2_1 = round(2**(base + random.gauss(0, 0.3)), 2)
            terf2_2 = round(2**(base + random.gauss(0, 0.3)), 2)
            igg_1 = round(2**(base + random.gauss(0, 0.3)), 2)
            igg_2 = round(2**(base + random.gauss(0, 0.3)), 2)
            
        writer.writerow([f"P_{prot}", prot, terf2_1, terf2_2, igg_1, igg_2])

print(f"Generated synthetic dataset with {len(all_proteins)} proteins at {out_path}")
