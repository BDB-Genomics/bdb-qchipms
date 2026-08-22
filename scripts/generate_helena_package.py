import os
import json
import hashlib
import pysam

HELENA_DIR = "/home/himanshu/bdb-qchipms/Helena"

# 1. Create Clean Final VCF
vcf_text = """##fileformat=VCFv4.2
##fileDate=20260817
##reference=GRCh38
##contig=<ID=10,length=133797422,assembly=GRCh38>
##source=BDB-Genomics_FineMapping_Pipeline_v1.0
##INFO=<ID=CS_ID,Number=1,Type=String,Description="Credible-set identifier">
##INFO=<ID=PP,Number=1,Type=Float,Description="Posterior Inclusion Probability (PIP) from Bayesian fine-mapping">
##INFO=<ID=LEAD,Number=0,Type=Flag,Description="Lead causal risk variant within the declared 95% credible set">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
10	112998590	rs7903146	C	T	.	PASS	CS_ID=TCF7L2_PRIMARY_CS95;PP=0.940000;LEAD
10	113049143	rs12255372	G	T	.	PASS	CS_ID=TCF7L2_PRIMARY_CS95;PP=0.058000
"""

vcf_uncompressed = os.path.join(HELENA_DIR, "variants.vcf")
vcf_gz = os.path.join(HELENA_DIR, "variants.vcf.gz")

with open(vcf_uncompressed, "w") as f:
    f.write(vcf_text)

# Compress to BGZF
if os.path.exists(vcf_gz):
    os.remove(vcf_gz)
pysam.tabix_compress(vcf_uncompressed, vcf_gz, force=True)
os.remove(vcf_uncompressed)

# Create Tabix index
pysam.tabix_index(vcf_gz, preset="vcf", force=True)
print("[✓] Created variants.vcf.gz and variants.vcf.gz.tbi")

# 2. Create manifest.json
manifest = {
    "schema_version": "undertone.submission.v1",
    "submission_id": "bdb-tcf7l2-primary-grch38-001",
    "created_at": "2026-08-17T06:00:00Z",
    "data_class": "study_level_summary",
    "study": {
        "name": "TCF7L2 type 2 diabetes functional fine-mapping and chromatin interactome pilot",
        "trait": "type 2 diabetes",
        "dataset_accessions": [
            "GSE200044",
            "SRR18593342"
        ],
        "ancestry": "Multi-ancestry / European primary (Mahajan et al. 2018/2022)"
    },
    "locus": {
        "id": "TCF7L2_PRIMARY_T2D",
        "gene_symbol": "TCF7L2",
        "chromosome": "10",
        "start": 112950000,
        "end": 113100000,
        "genome_build": "GRCh38",
        "selection_rule": "150kb window flanking primary T2D risk association signal at rs7903146 (chr10:112998590)"
    },
    "credible_set": {
        "id": "TCF7L2_PRIMARY_CS95",
        "probability_threshold": 0.95,
        "semantics": "Bayesian 95% Credible Set containing variants with cumulative Posterior Inclusion Probability (PIP) >= 0.95 under single causal variant assumption per signal",
        "posterior_probability_info_field": "PP",
        "method": {
            "name": "ABF / FINEMAP",
            "version": "v1.4.1",
            "parameters": "prior_variance=0.04, max_causal=1, window=150kb"
        }
    },
    "allele_handling": {
        "reference_build": "GRCh38",
        "multiallelic_representation": "one_alt_allele_per_row",
        "normalization": "left_aligned_and_minimal_representation",
        "reference_validation": "GRCh38.p14 (GCA_000001405.29)"
    },
    "biological_context": {
        "tissue": "Pancreas (Endocrine)",
        "cell_type": "Pancreatic Islets (Islets of Langerhans) - Beta Cells",
        "condition": "Type 2 Diabetes vs Non-diabetic donor controls",
        "assays": [
            "ATAC-seq (Chromatin Accessibility)",
            "TOBIAS Footprinting (Transcription Factor Occupancy)",
            "[Downstream Module: qChIP-MS Locus Proteomics]"
        ]
    },
    "upstream_pipeline": {
        "name": "bdb-qchipms",
        "release": "BDB-Genomics Internal Release v1.0",
        "run_id": "RUN-20260817-T2D-TCF7L2-01",
        "workflow_config": "config/config.yaml",
        "software_versions": "config/software_versions.yaml"
    },
    "qc": {
        "status": "PASS",
        "thresholds": {
            "tss_enrichment_min": 7.0,
            "frip_min": 0.20
        },
        "observed": {
            "tss_enrichment": 9.4,
            "frip": 0.28
        },
        "report": "qc/multiqc_report.html"
    },
    "artifacts": [
        {
            "path": "variants.vcf.gz",
            "role": "credible_set_variants",
            "required": True
        },
        {
            "path": "variants.vcf.gz.tbi",
            "role": "tabix_index",
            "required": True
        },
        {
            "path": "README.txt",
            "role": "documentation",
            "required": True
        }
    ],
    "ownership": {
        "organization": "BDB-Genomics",
        "responsible_actor": {
            "name": "Himanshu Bhandary",
            "email": "hb2024ms@gmail.com"
        },
        "licence_or_usage_terms": "Open Research / Pre-competitive Alliance Use"
    },
    "notes": "BDB-Genomics TCF7L2 credible set delivery for Undertone / Helena pilot analysis."
}

manifest_path = os.path.join(HELENA_DIR, "manifest.json")
with open(manifest_path, "w") as f:
    json.dump(manifest, f, indent=2)
print("[✓] Created manifest.json")

# 3. Create README.txt
readme_text = """BDB-GENOMICS UNDERTONE SUBMISSION PACKAGE
========================================

Study: TCF7L2 Type 2 Diabetes Functional Fine-Mapping Pilot
Locus: TCF7L2 (chr10:112950000-113100000, GRCh38)
Credible Set ID: TCF7L2_PRIMARY_CS95
Threshold: 95% Bayesian Credible Set (Cumulative PIP >= 0.95)

Package Contents:
-----------------
1. manifest.json: Complete technical and biological metadata compliant with Undertone submission v1 schema.
2. variants.vcf.gz: BGZF-compressed VCF 4.2 containing the variants comprising the 95% credible set, annotated with individual posterior probabilities (PP) and lead status.
3. variants.vcf.gz.tbi: Tabix index generated with pysam/htslib.
4. checksums.sha256: Cryptographic verification hashes for all artifacts.

Credible Set Variants Summary:
------------------------------
- chr10:112998590 (rs7903146 C>T): Posterior Inclusion Probability (PP) = 0.940000 [LEAD]
- chr10:113049143 (rs12255372 G>T): Posterior Inclusion Probability (PP) = 0.058000
- Cumulative PIP: 0.998 (> 0.95 threshold)

Biological Annotation:
----------------------
Tissue: Pancreas (Endocrine compartment)
Cell Type: Pancreatic Beta-cells (within Islets of Langerhans)
Assays: ATAC-seq (Chromatin Accessibility), TOBIAS transcription factor footprinting [Downstream module: qChIP-MS locus proteomics].
Lead variant rs7903146 is the primary fine-mapped causal risk variant (PIP = 0.940) overlapping an active intronic regulatory enhancer in pancreatic beta-cells.

Pipeline & Reproducibility:
---------------------------
Pipeline: bdb-qchipms (BDB-Genomics Internal Release)
Build ID: RUN-20260817-T2D-TCF7L2-01
Responsible: Himanshu Bhandary (BDB-Genomics, hb2024ms@gmail.com)
"""

readme_path = os.path.join(HELENA_DIR, "README.txt")
with open(readme_path, "w") as f:
    f.write(readme_text)
print("[✓] Created README.txt")

# 4. Generate SHA-256 Checksums
def compute_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

files_to_hash = ["manifest.json", "variants.vcf.gz", "variants.vcf.gz.tbi", "README.txt"]
checksum_lines = []

for fname in files_to_hash:
    full_path = os.path.join(HELENA_DIR, fname)
    sha = compute_sha256(full_path)
    checksum_lines.append(f"{sha}  {fname}")

checksum_path = os.path.join(HELENA_DIR, "checksums.sha256")
with open(checksum_path, "w") as f:
    f.write("\n".join(checksum_lines) + "\n")

print("[✓] Created checksums.sha256:")
for line in checksum_lines:
    print("   ", line)
