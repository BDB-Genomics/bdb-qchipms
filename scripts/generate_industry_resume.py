#!/usr/bin/env python3
"""
Generate an Industry / Corporate Resume for HR & Recruiter submissions
(Neuberg Diagnostics, clinical variant analyst roles, industry genomics).
Features: Zero GitHub links, zero repo URLs, full clinical/industry alignment,
clean 2-page layout, zero em dashes.
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
        footer_text = f"Himanshu Bhandary  |  Resume  |  Page {self._pageNumber} of {page_count}"
        self.drawCentredString(letter[0] / 2.0, 26, footer_text)
        self.restoreState()

def build_industry_resume(output_path: str):
    doc = SimpleDocTemplate(
        output_path,
        pagesize=letter,
        leftMargin=36,
        rightMargin=36,
        topMargin=26,
        bottomMargin=36
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
        fontSize=10.5,
        leading=13.5,
        textColor=colors.HexColor('#0F172A'),
        spaceBefore=6,
        spaceAfter=2,
        keepWithNext=True
    )

    role_title = ParagraphStyle(
        'CVRoleTitle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=9.5,
        leading=12.5,
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
        leading=11.5,
        textColor=colors.HexColor('#1F2937'),
        leftIndent=10,
        firstLineIndent=-10,
        spaceAfter=2
    )

    body_style = ParagraphStyle(
        'CVBody',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=8.5,
        leading=11.5,
        textColor=colors.HexColor('#1F2937'),
        spaceAfter=2
    )

    story = []

    # =========================================================
    # PAGE 1: HEADER, EDUCATION, SUMMARY, EXPERIENCE
    # =========================================================
    story.append(Paragraph("HIMANSHU BHANDARY", name_style))
    story.append(Spacer(1, 2))
    
    contact_text = (
        "Muzaffarpur, Bihar, India &nbsp;|&nbsp; "
        "Email: 2032ushimanshu@gmail.com &nbsp;|&nbsp; "
        "LinkedIn: linkedin.com/in/himanshu2032"
    )
    story.append(Paragraph(contact_text, contact_style))
    story.append(Spacer(1, 3))
    story.append(HRFlowable(width="100%", thickness=1.0, color=colors.HexColor("#0F172A"), spaceBefore=1, spaceAfter=4))

    # 1. EDUCATION
    story.append(Paragraph("EDUCATION", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=3))

    story.append(Paragraph("<b>Bachelor of Technology (B.Tech.) in Biotechnology</b> &nbsp;|&nbsp; GPA: 3.00 / 4.00 (75%, First Division)", role_title))
    story.append(Paragraph("Department of Biotechnology, Kathmandu University, Dhulikhel, Nepal &nbsp;|&nbsp; Sep 2019 - Jul 2024", role_meta))
    story.append(Paragraph("<b>Undergraduate Thesis:</b> <i>Generation of a Short Tandem Repeat (STR) Profile Targeting Eight Combined DNA Index System (CODIS) Loci</i>", bullet_style))
    story.append(Paragraph("<b>Supervisor:</b> Dr. Hitesh Bhattrai, Assistant Professor of Biotechnology", bullet_style))
    story.append(Spacer(1, 3))

    # 2. PROFESSIONAL SUMMARY
    story.append(Paragraph("PROFESSIONAL SUMMARY", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=3))
    
    summary_text = (
        "Biotechnology professional with over 1.5 years of hands-on experience in genomic variant curation, "
        "NGS data analysis, and reproducible bioinformatic workflows. Proficient in interpreting germline and non-coding variants "
        "using ClinVar, dbSNP, and Ensembl, with deep technical competency in Python, R, and automated pipeline execution."
    )
    story.append(Paragraph(summary_text, body_style))
    story.append(Spacer(1, 3))

    # 3. PROFESSIONAL & RESEARCH EXPERIENCE
    story.append(Paragraph("PROFESSIONAL & RESEARCH EXPERIENCE", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=3))

    # Item 1: BDB-Genomics / Helena Collaboration
    story.append(Paragraph("<b>Bioinformatics Specialist & Platform Developer</b>", role_title))
    story.append(Paragraph("BDB-Genomics (Collaborator: Helena Bioinformatics, Bulgaria) &nbsp;|&nbsp; May 2026 - Present", role_meta))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Non-Coding Variant Fine-Mapping:</b> Evaluated genetic variant effects at the <i>TCF7L2</i> locus, fine-mapping the 95% credible set (rs7903146, posterior probability = 0.940; rs12255372, posterior probability = 0.058) in human pancreatic beta-cell open chromatin regions.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Genomic Workflow Engineering:</b> Developed automated, modular Snakemake workflows for high-throughput sequencing data, incorporating adapter trimming, alignment, QC filtering, and variant peak calling.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Transcription Factor Dynamics:</b> Applied TOBIAS BINDetect with JASPAR motif databases to investigate chromatin accessibility and transcription factor occupancy across disease loci.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Platform Prototyping:</b> Prototyped a modular, resource-efficient multiomics analysis platform for resource-limited settings, formally endorsed via a Letter of Support by the JIDCuk Charitable Trust (UK, Jul 2026).",
        bullet_style
    ))
    story.append(Spacer(1, 3))

    # Item 2: Shodhaka Life Sciences & IBAB
    story.append(Paragraph("<b>Bioinformatics Researcher</b>", role_title))
    story.append(Paragraph("Shodhaka Life Sciences Pvt. Ltd. (In collaboration with IBAB, Bengaluru) &nbsp;|&nbsp; Mar 2025 - Apr 2026", role_meta))
    story.append(Paragraph("<i>Supervisor: Dr. Kshitish Acharya; Mentor: Jessica Evangeline KC</i>", role_meta))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Pipeline Automation:</b> Built a configurable Snakemake ATAC-seq pipeline for automated processing from raw FASTQ files through alignment and quality control metrics (FRiP score, TSS enrichment).",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Post-Alignment Processing:</b> Integrated blacklist region filtering, mitochondrial read removal, and coordinate shifting to generate analysis-ready genomic peak files.",
        bullet_style
    ))
    story.append(Spacer(1, 3))

    # Item 3: St. Jude Children's Research Hospital
    story.append(Paragraph("<b>Collaborative Researcher (Single-Cell Transcriptomics)</b>", role_title))
    story.append(Paragraph("St. Jude Children's Research Hospital &nbsp;|&nbsp; Sep 2024 - Feb 2025 &nbsp;|&nbsp; <i>Mentor: Dr. Yogesh Dhungana</i>", role_meta))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Single-Cell Transcriptomic Analysis:</b> Analyzed single-cell RNA sequencing data to examine cellular heterogeneity in human pancreatic islet and acinar endothelial cells across healthy and diabetic states.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Custom Analysis Pipelines:</b> Wrote modular R scripts utilizing Seurat and SingleCellExperiment for quality control filtering, sctransform normalization, high-variance gene selection, PCA, and UMAP/t-SNE visualization.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Differential Expression & Enrichment:</b> Performed differential expression analysis across endothelial subsets and conducted Gene Ontology and pathway enrichment to characterize vascular gene expression changes.",
        bullet_style
    ))
    story.append(Spacer(1, 3))

    # Item 4: Kathmandu University (Undergraduate Thesis)
    story.append(Paragraph("<b>Undergraduate Research Assistant (Thesis)</b>", role_title))
    story.append(Paragraph("Molecular Biology Laboratory, Kathmandu University, Nepal &nbsp;|&nbsp; Jan 2024 - Jul 2024 &nbsp;|&nbsp; <i>Supervisor: Dr. Hitesh Bhattrai</i>", role_meta))
    story.append(Paragraph(
        "&bull;&nbsp; <b>PCR Protocol Optimization:</b> Extracted genomic DNA from whole blood and optimized singleplex PCR conditions for 8 core CODIS STR loci (TH01, TPOX, D16S539, D7S820, D13S317, D21S11, D18S51, FGA).",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Experimental Validation:</b> Achieved 100% amplification across all 8 loci, verified amplicon sizing (100-140 bp) via 1% agarose gel electrophoresis, and documented results in an undergraduate thesis.",
        bullet_style
    ))
    story.append(Spacer(1, 3))

    # Item 5: MutS Evolutionary Analysis
    story.append(Paragraph("<b>Undergraduate Research Assistant</b>", role_title))
    story.append(Paragraph("Department of Biotechnology, Kathmandu University, Nepal &nbsp;|&nbsp; Mar 2022 - May 2022 &nbsp;|&nbsp; <i>Supervisor: Dr. Hitesh Bhattrai</i>", role_meta))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Sequence & Domain Analysis:</b> Retrieved 140 MutS homolog sequences from NCBI and UniProt to evaluate domain conservation in the DNA mismatch repair pathway across prokaryotes and eukaryotes.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Phylogenetic Reconstruction:</b> Conducted multiple sequence alignments (MUSCLE, SeaView), constructed bootstrap-supported maximum parsimony trees, and identified conserved terminal motifs.",
        bullet_style
    ))

    # Clean explicit page break
    story.append(PageBreak())

    # =========================================================
    # PAGE 2: SOFTWARE PACKAGES, SKILLS, RECOGNITION, TRAINING
    # =========================================================

    # 4. BIOINFORMATIC WORKFLOWS & SOFTWARE TOOLS
    story.append(Paragraph("BIOINFORMATIC WORKFLOWS & SOFTWARE TOOLS", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=4))
    
    story.append(Paragraph(
        "&bull;&nbsp; <b>bdb-atacseq:</b> Production-grade Snakemake pipeline for reproducible ATAC-seq data processing, Tn5 shift correction, peak calling (MACS2), and TSS/FRiP quality metrics.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>bdb-rnaseq:</b> Modular Snakemake workflow for end-to-end RNA-seq alignment (STAR), quantification (featureCounts), DESeq2 differential analysis, and MultiQC reporting.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>AlphaGenomeR:</b> R client interface bridging multimodal deep learning genomic predictions across chromatin accessibility, expression, and 3D structure.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>EvoR:</b> R package interfacing Evo 2 foundational models via NVIDIA BioNeMo for zero-shot variant effect prediction in non-coding regulatory sequences.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>bdb-qchipms:</b> Workflow for identifying locus-specific chromatin composition and label-free quantitative proteomics interactomes.",
        bullet_style
    ))
    story.append(Spacer(1, 5))

    # 5. TECHNICAL & LABORATORY SKILLS
    story.append(Paragraph("TECHNICAL & LABORATORY SKILLS", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=4))

    skills_data = [
        [
            Paragraph("<b>Genomic Variant Analysis:</b>", bullet_style),
            Paragraph("Germline variant interpretation, variant curation (ClinVar, dbSNP, Ensembl, UniProt), non-coding variant fine-mapping, TF motif analysis (TOBIAS, JASPAR), differential accessibility.", body_style)
        ],
        [
            Paragraph("<b>NGS & Computational Tools:</b>", bullet_style),
            Paragraph("Python, R (Bioconductor, ggplot2), Bash scripting, Snakemake, ATAC-seq, scRNA-seq (Seurat), bulk RNA-seq (STAR, DESeq2), ChIP-seq, peak calling (MACS2), Linux HPC.", body_style)
        ],
        [
            Paragraph("<b>Molecular & Wet-Lab:</b>", bullet_style),
            Paragraph("Genomic DNA/RNA isolation, plasmid extraction, PCR optimization (gradient/touchdown), restriction digestion, agarose gel electrophoresis, SDS-PAGE, protein quantification (Bradford, UV-Vis), aseptic cell culture.", body_style)
        ]
    ]

    skills_table = Table(skills_data, colWidths=[130, 414])
    skills_table.setStyle(TableStyle([
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1.5),
        ('TOPPADDING', (0,0), (-1,-1), 1.5),
        ('LEFTPADDING', (0,0), (-1,-1), 0),
        ('RIGHTPADDING', (0,0), (-1,-1), 0),
    ]))
    story.append(skills_table)
    story.append(Spacer(1, 5))

    # 6. PROFESSIONAL AFFILIATIONS & RECOGNITION
    story.append(Paragraph("PROFESSIONAL AFFILIATIONS & RECOGNITION", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=4))

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
    story.append(Spacer(1, 5))

    # 7. WORKSHOPS & CERTIFICATIONS
    story.append(Paragraph("WORKSHOPS & CERTIFICATIONS", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#CBD5E1"), spaceBefore=1, spaceAfter=4))

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

    # Build the document
    doc.build(story, canvasmaker=NumberedCanvas)
    print(f"Industry Resume successfully generated at: {output_path}")

if __name__ == "__main__":
    out_pdf = "/home/himanshu/bdb-qchipms/Himanshu_Bhandary_Resume_Industry.pdf"
    build_industry_resume(out_pdf)
