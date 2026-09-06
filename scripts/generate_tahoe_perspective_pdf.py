"""
Publication-grade, human-written, monochrome 1-page Technical Perspective PDF for Tahoe Therapeutics & Functional Genomics Teams.
Focuses on the 4 real screening failure modes (Ghost Peaks, Promoter Tug-of-War, Pre-Wired Quiescence, Bystander Trap).
Zero AI buzzwords, zero physics formulas, clean booktabs styling, perfectly balanced on 1 page.
"""

import sys
from pathlib import Path
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle

def generate_pdf():
    output_dir = Path("/home/himanshu/bdb-qchipms/DOCUMENTS/Presentation")
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = output_dir / "BDB_Genomics_Tahoe_Technical_Perspective.pdf"

    # Letter size: 612 x 792 pt. Margins: 36 pt
    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize=letter,
        leftMargin=36,
        rightMargin=36,
        topMargin=32,
        bottomMargin=32
    )

    styles = getSampleStyleSheet()
    
    # Pure Monochrome Palette
    c_black = colors.HexColor("#000000")
    c_charcoal = colors.HexColor("#1a1a1a")
    c_gray_dark = colors.HexColor("#333333")
    c_line = colors.HexColor("#000000")

    title_style = ParagraphStyle(
        'DocTitle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=13.5,
        leading=16.5,
        textColor=c_black
    )
    
    meta_style = ParagraphStyle(
        'DocMeta',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=8,
        leading=11,
        textColor=c_gray_dark
    )

    meta_style_right = ParagraphStyle(
        'DocMetaRight',
        parent=meta_style,
        alignment=2
    )

    section_heading = ParagraphStyle(
        'SectionHead',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=9.5,
        leading=12.5,
        textColor=c_black
    )

    body_text = ParagraphStyle(
        'BodyText',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=8,
        leading=11.5,
        textColor=c_charcoal
    )

    table_header = ParagraphStyle(
        'TableHead',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=7.5,
        leading=10,
        textColor=c_black
    )

    table_cell = ParagraphStyle(
        'TableCell',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=7.5,
        leading=10.5,
        textColor=c_charcoal
    )

    table_cell_bold = ParagraphStyle(
        'TableCellBold',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=7.5,
        leading=10.5,
        textColor=c_black
    )

    elements = []

    # 1. Header
    header_data = [
        [
            Paragraph("<b>BDB-GENOMICS</b> &nbsp;|&nbsp; Technical Perspective", title_style),
            Paragraph("FUNCTIONAL GENOMICS BRIEF", ParagraphStyle('RHead', parent=meta_style_right, fontName='Helvetica-Bold', textColor=c_black))
        ],
        [
            Paragraph("Beyond 1D Epigenomics: Why High-Confidence Enhancers Fail in CRISPRi Screens", ParagraphStyle('Sub', parent=meta_style, fontName='Helvetica-Bold', textColor=c_gray_dark)),
            Paragraph("September 2026 &nbsp;|&nbsp; Prepared for Target Discovery Teams", meta_style_right)
        ]
    ]
    header_table = Table(header_data, colWidths=[365, 175])
    header_table.setStyle(TableStyle([
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1),
        ('TOPPADDING', (0,0), (-1,-1), 0),
        ('LEFTPADDING', (0,0), (-1,-1), 0),
        ('RIGHTPADDING', (0,0), (-1,-1), 0),
    ]))
    elements.append(header_table)
    elements.append(Spacer(1, 4))
    elements.append(HRFlowable(width="100%", thickness=1.0, color=c_line, spaceBefore=2, spaceAfter=8))

    # 2. Executive Problem Statement
    intro_text = """
    <b>The Core Screening Challenge:</b> High-throughput CRISPRi and Perturb-seq platforms routinely filter target libraries using standard 1D epigenomic profiles: high ATAC accessibility, prominent H3K27ac enrichment, and strong transcription factor motifs. Yet in primary screens, <b>over 70% of these high-confidence candidate enhancers produce zero transcriptional or phenotypic consequence when repressed</b>. Relying on 1D peak calls conflates permissive chromatin with productive regulatory coupling, consuming substantial synthesis and validation bandwidth on biological dead-ends.
    """
    elements.append(Paragraph(intro_text, body_text))
    elements.append(Spacer(1, 8))

    # 3. The 4 Fundamental Mechanisms of Screening Failure
    elements.append(Paragraph("<b>Four Mechanisms Driving Enhancer Disconnection in Perturbation Screens</b>", section_heading))
    elements.append(Spacer(1, 3))

    points_data = [
        [
            Paragraph("<b>1. The 'Ghost Peak' (Permissive vs. Coupled):</b><br/>Accessibility is permissive, not functional. Without stable 3D loop anchoring to a target promoter, an open enhancer is functionally poised or redundant. Synthesizing guides across uncoupled accessible peaks accounts for the majority of inert hits.", body_text),
            Paragraph("<b>2. The Promoter Tug-of-War (Intra-Domain Competition):</b><br/>In multi-gene loop domains, enhancers do not distribute signal equally—promoters actively compete for contact. Repressing one target promoter often causes an unexpected compensatory surge in a flanking gene, frequently misidentified as off-target toxicity.", body_text)
        ],
        [
            Paragraph("<b>3. Pre-Wired Quiescence (Basal vs. Triggered):</b><br/>Many disease-critical chromatin loops are physically pre-formed in basal cell architecture, but the promoter remains silent until activated by acute stress (cytokines, hypoxia, metabolic shock). Screening solely in basal states discards the most potent disease targets.", body_text),
            Paragraph("<b>4. The Bystander Phenotype Trap:</b><br/>In single-cell Perturb-seq readouts, an adjacent gene often displays minor transcriptional fluctuation while the true causal phenotype is driven by a distal transcriptional hub. 1D nearest-gene assignment causes teams to develop assays around the wrong protein.", body_text)
        ]
    ]
    points_table = Table(points_data, colWidths=[264, 264])
    points_table.setStyle(TableStyle([
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('TOPPADDING', (0,0), (-1,-1), 4),
        ('BOTTOMPADDING', (0,0), (-1,-1), 5),
        ('LEFTPADDING', (0,0), (-1,-1), 4),
        ('RIGHTPADDING', (0,0), (-1,-1), 4),
        ('LINEBELOW', (0,0), (-1,0), 0.5, colors.HexColor("#dddddd")),
    ]))
    elements.append(points_table)
    elements.append(Spacer(1, 10))

    # 4. Comparative Framework Table
    elements.append(Paragraph("<b>Operational Impact: 1D Epigenomics vs. 3D Domain Engagement Filtering</b>", section_heading))
    elements.append(Spacer(1, 3))

    comp_data = [
        [Paragraph("<b>Screening Parameter</b>", table_header),
         Paragraph("<b>Standard 1D Pipeline (ATAC + ChIP)</b>", table_header),
         Paragraph("<b>3D Regulatory Domain Pre-Filter</b>", table_header)],
        [Paragraph("<b>Library Candidate Selection</b>", table_cell_bold),
         Paragraph("Tiles all broad ATAC / H3K27ac peaks indiscriminately", table_cell),
         Paragraph("<b>Filters exclusively for physically coupled, stable loop anchors</b>", table_cell_bold)],
        [Paragraph("<b>Multi-Promoter Ambiguity</b>", table_cell_bold),
         Paragraph("Defaults to nearest TSS; blind to domain competition", table_cell),
         Paragraph("<b>Models competitive contact dynamics and compensatory surge</b>", table_cell_bold)],
        [Paragraph("<b>Basal State Screening</b>", table_cell_bold),
         Paragraph("Discards silent promoters that lack basal transcription", table_cell),
         Paragraph("<b>Identifies pre-wired structural loops poised for inducible firing</b>", table_cell_bold)],
        [Paragraph("<b>Screening Efficiency</b>", table_cell_bold),
         Paragraph("High attrition; ~70% of synthesized guides yield no phenotype", table_cell),
         Paragraph("<b>Eliminates inert 'ghost' peaks prior to guide library synthesis</b>", table_cell_bold)]
    ]
    t_comp = Table(comp_data, colWidths=[120, 210, 210])
    t_comp.setStyle(TableStyle([
        ('LINEABOVE', (0,0), (-1,0), 1.0, c_black),
        ('LINEBELOW', (0,0), (-1,0), 0.5, c_black),
        ('LINEBELOW', (0,-1), (-1,-1), 1.0, c_black),
        ('TOPPADDING', (0,0), (-1,-1), 3.5),
        ('BOTTOMPADDING', (0,0), (-1,-1), 3.5),
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
    ]))
    elements.append(t_comp)
    elements.append(Spacer(1, 10))

    # 5. Technical Scope & Contact Box
    scope_text = """
    <b>Pilot Evaluation for Target Discovery Teams:</b> If you have a panel of 10 high-confidence ATAC/H3K27ac candidate peaks that failed to produce expected transcriptional or phenotypic effects in CRISPRi screens, BDB-Genomics can evaluate their 3D loop coupling, domain insulation, and promoter competition profiles under a standard mutual CDA.
    <br/><br/>
    <b>Inquiries & Technical Pilot:</b> Himanshu Bhandary &nbsp;|&nbsp; BDB-Genomics &nbsp;|&nbsp; <b>2032ushimanshu@gmail.com</b>
    """
    scope_table = Table([[Paragraph(scope_text, ParagraphStyle('ScopeStyle', parent=body_text, fontSize=7.5, leading=11))]], colWidths=[540])
    scope_table.setStyle(TableStyle([
        ('LINEABOVE', (0,0), (-1,-1), 0.75, c_line),
        ('LINEBELOW', (0,0), (-1,-1), 0.75, c_line),
        ('TOPPADDING', (0,0), (-1,-1), 5),
        ('BOTTOMPADDING', (0,0), (-1,-1), 5),
        ('LEFTPADDING', (0,0), (-1,-1), 8),
        ('RIGHTPADDING', (0,0), (-1,-1), 8),
    ]))
    elements.append(scope_table)

    doc.build(elements)
    print(f"✅ Generated 1-page Tahoe Technical Perspective PDF: {pdf_path}")

if __name__ == "__main__":
    generate_pdf()
