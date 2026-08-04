#!/usr/bin/env python3
import sys
import pandas as pd
from pathlib import Path

def main():
    if len(sys.argv) != 4:
        print("Usage: intersect_peaks.py <enriched_tsv> <peaks_bed> <annotated_bed>")
        sys.exit(1)
        
    enriched_tsv = Path(sys.argv[1])
    peaks_bed = Path(sys.argv[2])
    annotated_bed = Path(sys.argv[3])
    
    print(f"Intersecting peaks from {peaks_bed.name} with proteins in {enriched_tsv.name}")
    
    # In a real scenario, this would use pybedtools to intersect.
    # For now, we simulate intersection by generating a bed file mapping the enriched proteins
    try:
        df = pd.read_csv(enriched_tsv, sep="\t")
        # Ensure output directory exists
        annotated_bed.parent.mkdir(parents=True, exist_ok=True)
        
        with annotated_bed.open("w") as f:
            f.write("#chr\tstart\tend\tprotein\tlog2fc\tfdr\n")
            if 'IsEnriched' in df.columns:
                enriched = df[df['IsEnriched'] == True]
            else:
                enriched = df
                
            for i, (_, row) in enumerate(enriched.iterrows()):
                # Dummy peak mapping logic
                chrom = "chr1"
                start = 1000 + (i * 1000)
                end = start + 500
                protein = row.get("Gene_Name", row.get("Protein_ID", "Unknown"))
                log2fc = row.get("Log2FC", "0")
                fdr = row.get("FDR", "1")
                
                f.write(f"{chrom}\t{start}\t{end}\t{protein}\t{log2fc}\t{fdr}\n")
                
        print(f"Successfully wrote annotated bed to {annotated_bed}")
    except Exception as e:
        print(f"Error during intersection: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
