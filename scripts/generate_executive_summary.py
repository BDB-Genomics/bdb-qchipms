#!/usr/bin/env python3
"""
Generate a strict Black & White, 1-Page Publication-Grade Executive Summary
for BDB-Genomics (Investor & Grant Submission).
Follows exact same typography, geometry, and layout metrics as the Academic CV.
Zero em dashes, exact 1-page fit.
"""

import os
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, HRFlowable, Table, TableStyle
from reportlab.pdfgen import canvas

class SinglePageCanvas(canvas.Canvas):
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
        self.setFont("Helvetica", 7.5)
        self.setFillColor(colors.HexColor("#475569"))
        footer_text = "BDB-Genomics  |  Executive Summary  |  Confidential"
        self.drawCentredString(letter[0] / 2.0, 16, footer_text)
        self.restoreState()

def build_executive_summary(output_path: str):
    doc = SimpleDocTemplate(
        output_path,
        pagesize=letter,
        leftMargin=36,
        rightMargin=36,
        topMargin=22,
        bottomMargin=26
    )

    styles = getSampleStyleSheet()

    company_title = ParagraphStyle(
        'ESCompanyTitle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=15,
        leading=17,
        textColor=colors.HexColor('#000000'),
        alignment=1
    )

    subtitle_style = ParagraphStyle(
        'ESSubtitle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=8.5,
        leading=11,
        textColor=colors.HexColor('#1E293B'),
        alignment=1
    )

    contact_style = ParagraphStyle(
        'ESContact',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=7.5,
        leading=10,
        textColor=colors.HexColor('#334155'),
        alignment=1
    )

    section_heading = ParagraphStyle(
        'ESSectionHeading',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=9,
        leading=11,
        textColor=colors.HexColor('#000000'),
        spaceBefore=3,
        spaceAfter=1,
        keepWithNext=True
    )

    body_style = ParagraphStyle(
        'ESBody',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=7.8,
        leading=10.2,
        textColor=colors.HexColor('#000000'),
        spaceAfter=1.5
    )

    bullet_style = ParagraphStyle(
        'ESBullet',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=7.8,
        leading=10.2,
        textColor=colors.HexColor('#000000'),
        leftIndent=8,
        firstLineIndent=-8,
        spaceAfter=1.5
    )

    table_cell = ParagraphStyle(
        'ESTableCell',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=7.2,
        leading=9.2,
        textColor=colors.HexColor('#000000')
    )

    table_header = ParagraphStyle(
        'ESTableHeader',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=7.5,
        leading=9.5,
        textColor=colors.HexColor('#000000')
    )

    story = []

    # HEADER
    story.append(Paragraph("BDB-GENOMICS", company_title))
    story.append(Spacer(1, 1))
    story.append(Paragraph("Computational Epigenomics & Chromatin Interactomics Infrastructure", subtitle_style))
    story.append(Spacer(1, 1))
    
    contact_text = (
        "Founder & Lead Developer: Himanshu Bhandary &nbsp;|&nbsp; "
        "Email: 2032ushimanshu@gmail.com &nbsp;|&nbsp; "
        "LinkedIn: linkedin.com/in/himanshu2032 &nbsp;|&nbsp; "
        "GitHub: github.com/BDB-Genomics"
    )
    story.append(Paragraph(contact_text, contact_style))
    story.append(Spacer(1, 2))
    story.append(HRFlowable(width="100%", thickness=1.0, color=colors.HexColor("#000000"), spaceBefore=1, spaceAfter=2.5))

    # 1. PROBLEM & INDUSTRY CONTEXT
    story.append(Paragraph("1. PROBLEM & INDUSTRY CONTEXT", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#000000"), spaceBefore=1, spaceAfter=2))
    story.append(Paragraph(
        "Understanding disease etiology requires analyzing not only where regulatory proteins bind along 1D chromatin, but resolving the physical multi-protein complexes and interactomes assembled at specific genomic loci. While sequencing assays (ChIP-seq, ATAC-seq) identify linear coordinates, quantitative ChIP coupled with Mass Spectrometry (qChIP-MS) is necessary to determine local protein composition. Currently, integrating label-free quantitative (LFQ) proteomics with peak calls relies on fragmented in-house scripts, causing severe reproducibility bottlenecks and high overhead for translational teams.",
        body_style
    ))
    story.append(Spacer(1, 1.5))

    # 2. CORE TECHNOLOGY & INFRASTRUCTURE
    story.append(Paragraph("2. CORE TECHNOLOGY & SUITE", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#000000"), spaceBefore=1, spaceAfter=2))
    story.append(Paragraph(
        "&bull;&nbsp; <b>bdb-qchipms:</b> Automated Snakemake workflow merging label-free quantitative proteomics intensities with MACS2 ChIP-seq peaks, generating locus-annotated tracks, volcano statistical enrichment plots, and automated QC reports.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>bdb-atacseq & bdb-rnaseq:</b> Production-grade containerized pipelines for high-throughput sequencing with automated adapter trimming, alignment, Tn5 shift correction, and TSS enrichment metrics.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>AlphaGenomeR & EvoR:</b> High-performance R interface packages connecting multimodal deep learning and Evo 2 foundation models via NVIDIA BioNeMo for zero-shot variant effect prediction in non-coding regulatory elements.",
        bullet_style
    ))
    story.append(Spacer(1, 1.5))

    # 3. TARGET USERS & VALUE PROPOSITION
    story.append(Paragraph("3. TARGET USERS & VALUE PROPOSITION", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#000000"), spaceBefore=1, spaceAfter=2))

    table_data = [
        [
            Paragraph("<b>Target Audience</b>", table_header),
            Paragraph("<b>Key Industry Pain Point</b>", table_header),
            Paragraph("<b>BDB-Genomics Solution</b>", table_header)
        ],
        [
            Paragraph("<b>Academic Labs</b>", table_cell),
            Paragraph("Lack dedicated bioinformatics engineers to process complex mass spectrometry and sequencing data.", table_cell),
            Paragraph("Turnkey, reproducible Snakemake workflows reducing turnaround from weeks to hours.", table_cell)
        ],
        [
            Paragraph("<b>Biopharma & Discovery</b>", table_cell),
            Paragraph("Need rapid identification of novel chromatin remodelers, co-factors, and target validation.", table_cell),
            Paragraph("Standardized locus-level protein interactome mapping with statistical enrichment rigor.", table_cell)
        ],
        [
            Paragraph("<b>Core Facilities & CROs</b>", table_cell),
            Paragraph("High maintenance compute costs and inconsistent ad-hoc scripting pipelines.", table_cell),
            Paragraph("Containerized, version-controlled workflows with automated multi-sample reporting.", table_cell)
        ]
    ]

    summary_table = Table(table_data, colWidths=[95, 215, 230])
    summary_table.setStyle(TableStyle([
        ('BOX', (0,0), (-1,-1), 0.5, colors.HexColor('#000000')),
        ('INNERGRID', (0,0), (-1,-1), 0.3, colors.HexColor('#CBD5E1')),
        ('TOPPADDING', (0,0), (-1,-1), 2),
        ('BOTTOMPADDING', (0,0), (-1,-1), 2),
        ('LEFTPADDING', (0,0), (-1,-1), 4),
        ('RIGHTPADDING', (0,0), (-1,-1), 4),
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('LINEBELOW', (0,0), (-1,0), 0.8, colors.HexColor('#000000')),
    ]))
    story.append(summary_table)
    story.append(Spacer(1, 2))

    # 4. CURRENT TRACTION & STRATEGIC VALIDATION
    story.append(Paragraph("4. TRACTION & STRATEGIC VALIDATION", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#000000"), spaceBefore=1, spaceAfter=2))
    story.append(Paragraph(
        "&bull;&nbsp; <b>European Commercial Pilot:</b> Active research collaboration with Helena Bioinformatics (Bulgaria), integrating upstream variant fine-mapping and chromatin accessibility into their clinical Undertone platform.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>International Support:</b> Formal Letter of Support and project endorsement from JIDCuk Charitable Trust (UK, Jul 2026) for developing resource-efficient, accessible multiomics workflows.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Production-Ready Codebases:</b> Open-source repositories actively maintained with schema validation, automated testing, and comprehensive technical documentation.",
        bullet_style
    ))
    story.append(Spacer(1, 1.5))

    # 5. STRATEGIC OBJECTIVES & PRE-SEED ROADMAP
    story.append(Paragraph("5. STRATEGIC OBJECTIVES & PRE-SEED ROADMAP", section_heading))
    story.append(HRFlowable(width="100%", thickness=0.4, color=colors.HexColor("#000000"), spaceBefore=1, spaceAfter=2))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Corporate Incorporation:</b> Completing formal entity registration (UK / India) and incubator onboarding.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>High-Performance Compute:</b> Securing dedicated GPU/CPU cloud infrastructure for foundation model inference.",
        bullet_style
    ))
    story.append(Paragraph(
        "&bull;&nbsp; <b>Commercial Expansion:</b> Scaling biopharma pilot partnerships and submitting non-dilutive grant applications.",
        bullet_style
    ))

    # Build the document
    doc.build(story, canvasmaker=SinglePageCanvas)
    print(f"Publication-Grade B&W Executive Summary successfully generated at: {output_path}")

if __name__ == "__main__":
    out_pdf = "/home/himanshu/bdb-qchipms/BDB_Genomics_Executive_Summary.pdf"
    build_executive_summary(out_pdf)
