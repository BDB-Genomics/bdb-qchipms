#!/usr/bin/env python3
"""
Generate a publication-grade Academic CV for European, Canadian, UK & US
Graduate Admissions (Bioinformatics, Computational Biology, Human Genetics).
Zero orphan headers, zero corporate jargon, zero em dashes, exact 2-page academic proportion.
"""

import os
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, HRFlowable, Table, TableStyle, PageBreak
from reportlab.pdfgen import canvas

class NumberedCanvas(canvas.Canvas):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._saved_page_states = []

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        num_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self.draw_footer(num_pages)
            super().showPage()
        super().save()

    def draw_footer(self, page_count):
        self.saveState()
        self.setFont("Helvetica", 8)
        self.setFillColor(colors.HexColor("#64748B"))
        footer_text = f"Himanshu Bhandary  |  Curriculum Vitae  |  Page {self._pageNumber} of {page_count}"
        self.drawCentredString(letter[0] / 2.0, 24, footer_text)
        self.restoreState()

def build_academic_cv(output_path: str):
    doc = SimpleDocTemplate(
        output_path,
        pagesize=letter,
        leftMargin=34,
        rightMargin=34,
        topMargin=24,
        bottomMargin=34
    )

    styles = getSampleStyleSheet()

    name_style = ParagraphStyle(
        'CVName',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=18,
        leading=21,
        textColor=colors.HexColor('#0F172A'),
        alignment=1
    )

    contact_style = ParagraphStyle(
        'CVContact',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=8.5,
        leading=11.5,
        textColor=colors.HexColor('#334155'),
        alignment=1
    )

    section_heading = ParagraphStyle(
        'CVSectionHeading',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=10.0,
        leading=13.0,
        textColor=colors.HexColor('#0F172A'),
        spaceBefore=5,
        spaceAfter=2,
        keepWithNext=True
    )

    role_title = ParagraphStyle(
        'CVRoleTitle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=9.0,
        leading=12.0,
        textColor=colors.HexColor('#1E293B'),
        keepWithNext=True
    )

    role_meta = ParagraphStyle(
        'CVRoleMeta',
        parent=styles['Normal'],
        fontName='Helvetica-Oblique',
        fontSize=8.5,
        leading=11.5,
        textColor=colors.HexColor('#475569'),
        keepWithNext=True
    )

    bullet_style = ParagraphStyle(
        'CVBullet',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=8.5,
        leading=11.2,
        textColor=colors.HexColor('#1F2937'),
        leftIndent=10,
        firstLineIndent=-10,
        spaceAfter=1.8
    )

    body_style = ParagraphStyle(
        'CVBody',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=8.5,
        leading=11.2,
        textColor=colors.HexColor('#1F2937'),
        spaceAfter=2
    )

    story = []

    # =========================================================
    # PAGE 1: HEADER, EDUCATION, INTERESTS, RESEARCH EXPERIENCE
    # =========================================================
    story.append(Paragraph("HIMANSHU BHANDARY", name_style))
    story.append(Spacer(1, 2))
    
    contact_text = (
        "Email: 2032ushimanshu@gmail.com &nbsp;|&nbsp; "
        "GitHub: github.com/BDB-Genomics &nbsp;|&nbsp; "
        "LinkedIn: linkedin.com/in/himanshu2032 &nbsp;|&nbsp; "
        "Location: Bengaluru, India / Kathmandu, Nepal"
    )
    story.append(Paragraph(contact_text, contact_style))
    story.append(Spacer(1, 2))
    story.append(HRFlowable(width="100%", thickness=1.0, color=colors.HexColor("#0F172A"), spaceBefore=1, spaceAfter=3))

    # 1. EDUCATION
    story.append(Paragraph("EDUCATION", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=3))

    story.append(Paragraph("<b>Bachelor of Technology (B.Tech.) in Biotechnology</b> &nbsp;|&nbsp; GPA: 3.00 / 4.00 (First Division)", role_title))
    story.append(Paragraph("Kathmandu University, School of Science, Dhulikhel, Nepal &nbsp;|&nbsp; Sep 2019 - Jul 2024", role_meta))
    story.append(Paragraph("<b>Language of Instruction:</b> English (<b>Medium of Instruction / MOI Certified</b>)", bullet_style))
    story.append(Paragraph("<b>Undergraduate Thesis:</b> <i>Generation of an STR Profile Targeting Eight Combined DNA Index System (CODIS) Loci</i>", bullet_style))
    story.append(Paragraph("<b>Supervisor:</b> Dr. Hitesh Bhattrai, Assistant Professor of Biotechnology", bullet_style))
    story.append(Spacer(1, 2))

    # 2. RESEARCH INTERESTS
    story.append(Paragraph("RESEARCH INTERESTS", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=3))
    
    interests_text = (
        "Computational Genomics, Post-GWAS Statistical Fine-Mapping, Chromatin & Epigenomic Regulation (ATAC-seq, scRNA-seq), "
        "Transcription Factor Dynamics, 3D Enhancer-Promoter Looping (ABC Modeling), and Reproducible Pipeline Engineering (Snakemake, Nextflow)."
    )
    story.append(Paragraph(interests_text, body_style))
    story.append(Spacer(1, 2))

    # 3. RESEARCH EXPERIENCE
    story.append(Paragraph("RESEARCH EXPERIENCE", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=3))

    # Item 1: BDB-Genomics / Helena Collaboration
    story.append(Paragraph("<b>Lead Computational Biologist & Platform Architect</b>", role_title))
    story.append(Paragraph("BDB-Genomics (Collaborator: Helena Bioinformatics, Bulgaria) &nbsp;|&nbsp; May 2026 - Present", role_meta))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Multi-Omic Workflow Engineering:</b> Developed automated, modular Snakemake pipelines for ATAC-seq, RNA-seq, and locus-specific quantitative proteomics (`bdb-qchipms`), executing end-to-end data ingestion, alignment, and statistical QC.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Transcription Factor Footprint Modeling:</b> Applied TOBIAS and position weight matrix (PWM) scoring engines against JASPAR and the 2026 Nature Human TF Codebook to quantify allelic binding affinity disruption deltas.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Statistical Fine-Mapping & 3D Enhancer Linking:</b> Fine-mapped metabolic disease 95% credible sets (DIAMANTE, 1000 Genomes LD) and integrated ENCODE-rE2G Activity-by-Contact (ABC) models to link distal non-coding enhancers to disease promoters.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Global Health Platform Prototyping:</b> Prototyped a modular, resource-efficient multiomics analysis platform for low-resource settings, formally endorsed via a Letter of Support by the JIDCuk Charitable Trust (UK, Jul 2026).",
        bullet_style
    ))
    story.append(Spacer(1, 2))

    # Item 2: Shodhaka Life Sciences & IBAB
    story.append(Paragraph("<b>Bioinformatics Researcher</b>", role_title))
    story.append(Paragraph("Shodhaka Life Sciences Pvt. Ltd. (In collaboration with IBAB, Bengaluru) &nbsp;|&nbsp; Mar 2025 - Apr 2026", role_meta))
    story.append(Paragraph("<i>Supervisor: Dr. Kshitish Acharya; Mentor: Jessica Evangeline KC</i>", role_meta))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Pipeline Automation:</b> Built a configurable Snakemake ATAC-seq pipeline for automated processing from raw FASTQ files through alignment and quality control metrics (FRiP score, TSS enrichment).",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Post-Alignment Quality Control:</b> Integrated ENCODE blacklist region filtering, mitochondrial read removal, and Tn5 coordinate shifting to generate publication-grade consensus peak files.",
        bullet_style
    ))
    story.append(Spacer(1, 2))

    # Item 3: St. Jude Children's Research Hospital
    story.append(Paragraph("<b>Collaborative Researcher (Single-Cell Transcriptomics)</b>", role_title))
    story.append(Paragraph("St. Jude Children's Research Hospital &nbsp;|&nbsp; Sep 2024 - Feb 2025 &nbsp;|&nbsp; <i>Mentor: Dr. Yogesh Dhungana</i>", role_meta))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Single-Cell Transcriptomic Analysis:</b> Analyzed single-cell RNA sequencing data to examine cellular heterogeneity in human pancreatic islet and acinar endothelial cells across healthy and diabetic states.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Custom Analysis Pipelines:</b> Wrote modular R scripts utilizing Seurat and SingleCellExperiment for quality control filtering, sctransform normalization, high-variance gene selection, PCA, and UMAP visualization.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Differential Expression & Enrichment:</b> Performed differential expression analysis across endothelial subsets and conducted Gene Ontology and pathway enrichment to characterize vascular gene expression changes.",
        bullet_style
    ))
    story.append(Spacer(1, 2))

    # Item 4: Kathmandu University (Undergraduate Thesis)
    story.append(Paragraph("<b>Undergraduate Research Assistant (Thesis)</b>", role_title))
    story.append(Paragraph("Molecular Biology Laboratory, Kathmandu University, Nepal &nbsp;|&nbsp; Jan 2024 - Jul 2024 &nbsp;|&nbsp; <i>Supervisor: Dr. Hitesh Bhattrai</i>", role_meta))
    story.append(Paragraph(
        "&bull;&nbsp; <b>PCR Protocol Optimization:</b> Extracted genomic DNA from whole blood and optimized singleplex PCR conditions for 8 core CODIS STR loci (TH01, TPOX, D16S539, D7S820, D13S317, D21S11, D18S51, FGA).",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Experimental Validation:</b> Verified amplicon sizing (100-140 bp) via 1% agarose gel electrophoresis with 100% amplification success, documenting results in an undergraduate thesis.",
        bullet_style
    ))

    # Explicit page break to guarantee Page 2 layout
    story.append(PageBreak())

    # =========================================================
    # PAGE 2: SOFTWARE, SKILLS, LEADERSHIP, TRAINING, REFERENCES
    # =========================================================

    # 4. OPEN-SOURCE SOFTWARE & PACKAGES
    story.append(Paragraph("OPEN-SOURCE SOFTWARE & WORKFLOWS", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=3))
    
    story.append(Paragraph(
        "&bull;&nbsp; <b>bdb-qchipms:</b> End-to-end Snakemake workflow for quantitative locus-specific chromatin composition and label-free mass spectrometry proteomics interactomes. (github.com/BDB-Genomics/bdb-qchipms)",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>bdb-atacseq:</b> Production-grade Snakemake pipeline for reproducible ATAC-seq data processing, Tn5 shift correction, MACS2 peak calling, and TSS/FRiP quality metrics. (github.com/BDB-Genomics/bdb-atacseq)",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>bdb-rnaseq:</b> Modular Snakemake workflow for bulk RNA-seq alignment (STAR), quantification (featureCounts), DESeq2 differential analysis, and MultiQC reporting. (github.com/BDB-Genomics/bdb-rnaseq)",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>AlphaGenomeR:</b> R/Bioconductor interface bridging multimodal deep learning genomic predictions across chromatin accessibility, expression, and 3D contact models. (github.com/BDB-Genomics/AlphaGenomeR)",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>EvoR:</b> R package interfacing Evo 2 foundational models via NVIDIA BioNeMo for zero-shot variant effect prediction in non-coding regulatory sequences. (github.com/BDB-Genomics/EvoR)",
        bullet_style
    ))
    story.append(Spacer(1, 3))

    # 5. TECHNICAL SKILLS
    story.append(Paragraph("TECHNICAL SKILLS", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=3))

    skills_data = [
        [
            Paragraph("<b>Computational Genomics:</b>", bullet_style),
            Paragraph("ATAC-seq, scRNA-seq (Seurat), bulk RNA-seq (STAR, DESeq2), ChIP-seq, peak calling (MACS2), TF motif analysis (TOBIAS, JASPAR), Bayesian fine-mapping (SuSiE), Activity-by-Contact (ABC) 3D modeling.", body_style)
        ],
        [
            Paragraph("<b>Languages & HPC:</b>", bullet_style),
            Paragraph("Python (NumPy, SciPy, Pandas, PyRanges), R (Bioconductor, GenomicRanges, ggplot2), Bash/Shell scripting, Snakemake, Nextflow, Docker, Conda, Git/GitHub, Linux HPC (SLURM).", body_style)
        ],
        [
            Paragraph("<b>Molecular & Wet-Lab:</b>", bullet_style),
            Paragraph("Genomic DNA/RNA isolation, plasmid extraction, PCR optimization (gradient/touchdown), restriction digestion, agarose gel electrophoresis, SDS-PAGE, Bradford protein assay, bacterial transformation.", body_style)
        ]
    ]

    skills_table = Table(skills_data, colWidths=[120, 424])
    skills_table.setStyle(TableStyle([
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1.2),
        ('TOPPADDING', (0,0), (-1,-1), 1.2),
        ('LEFTPADDING', (0,0), (-1,-1), 0),
        ('RIGHTPADDING', (0,0), (-1,-1), 0),
    ]))
    story.append(skills_table)
    story.append(Spacer(1, 3))

    # 6. LEADERSHIP & RECOGNITION
    story.append(Paragraph("LEADERSHIP & RECOGNITION", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=3))

    story.append(Paragraph(
        "&bull;&nbsp; <b>Project Endorsement & Letter of Support</b>, JIDCuk Charitable Trust, UK (Jul 2026): Formally recognized and supported for developing open-source, resource-efficient multiomics workflows for LMIC research settings.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Founding Member, TP53NexusIndia</b> (2025 - Present): Participating in a national collaborative network of researchers, clinicians, and data scientists studying TP53 biology.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Joint Treasurer, Kathmandu University Youth Red Cross Circle</b> (2020 - 2021): Managed budgeting and logistics for campus health awareness programs and blood donation drives.",
        bullet_style
    ))
    story.append(Spacer(1, 3))

    # 7. WORKSHOPS & TRAINING
    story.append(Paragraph("WORKSHOPS & TRAINING", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=3))

    story.append(Paragraph(
        "&bull;&nbsp; <b>Metabolic Modelling: Tools, Techniques, and Applications:</b> RIBB and University of Nepal Board (Sep 2024). Computational metabolic network modeling and flux balance analysis.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Introduction to Metabolomics Data Science:</b> Omics Logic (Jun 2024). Mass spectrometry data processing and metabolic classification algorithms.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Academic Research Methodologies:</b> Nepal Research and Collaboration Center (Jun - Sep 2024). Experimental design, scientific documentation, and proposal authoring.",
        bullet_style
    ))
    story.append(Spacer(1, 3))

    # 8. ACADEMIC REFERENCES
    story.append(Paragraph("ACADEMIC REFERENCES", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=3))

    ref_text = (
        "<b>1. Dr. Kshitish Acharya</b>, Director & Research Scientist, Shodhaka Life Sciences / IBAB Bengaluru &nbsp;|&nbsp; "
        "<b>2. Dr. Hitesh Bhattrai</b>, Assistant Professor, Department of Biotechnology, Kathmandu University &nbsp;|&nbsp; "
        "<b>3. Dr. Yogesh Dhungana</b>, Postdoctoral Researcher, St. Jude Children's Research Hospital (Full contact details available upon request)."
    )
    story.append(Paragraph(ref_text, body_style))

    # Build the document
    doc.build(story, canvasmaker=NumberedCanvas)
    print(f"[SUCCESS] Flawless academic CV generated at: {output_path}")

if __name__ == "__main__":
    out_pdf = "/home/himanshu/bdb-qchipms/Himanshu_Bhandary_Academic_CV.pdf"
    doc_pdf = "/home/himanshu/bdb-qchipms/DOCUMENTS/Himanshu_Bhandary_Academic_CV.pdf"
    os.makedirs(os.path.dirname(doc_pdf), exist_ok=True)
    build_academic_cv(out_pdf)
    build_academic_cv(doc_pdf)
