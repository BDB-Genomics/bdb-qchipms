"""
Publication-grade, monochrome (black & white), high-contrast 1-page PDF document for LinkedIn posting targeting Tahoe Therapeutics & functional genomics teams.
Zero AI buzzwords, zero physics formula leaks, large readable typography for LinkedIn document preview.
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
    pdf_path = output_dir / "BDB_Genomics_Tahoe_Target_Brief.pdf"

    # Letter size: 612 x 792 pt. Margins: 36 pt
    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize=letter,
        leftMargin=36,
        rightMargin=36,
        topMargin=34,
        bottomMargin=34
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
        fontSize=15,
        leading=18,
        textColor=c_black
    )
    
    meta_style = ParagraphStyle(
        'DocMeta',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=8.5,
        leading=11.5,
        textColor=c_gray_dark
    )

    meta_style_right = ParagraphStyle(
        'DocMetaRight',
        parent=meta_style,
        alignment=2
    )

    main_heading = ParagraphStyle(
        'MainHead',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=11,
        leading=14.5,
        textColor=c_black
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
        fontSize=8.5,
        leading=12,
        textColor=c_charcoal
    )

    table_header = ParagraphStyle(
        'TableHead',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=8,
        leading=10.5,
        textColor=c_black
    )

    table_cell = ParagraphStyle(
        'TableCell',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=8,
        leading=11,
        textColor=c_charcoal
    )

    table_cell_bold = ParagraphStyle(
        'TableCellBold',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=8,
        leading=11,
        textColor=c_black
    )

    elements = []

    # 1. Header (Institutional & High Contrast)
    header_data = [
        [
            Paragraph("<b>BDB-GENOMICS</b> &nbsp;|&nbsp; Technical Perspective", title_style),
            Paragraph("FUNCTIONAL GENOMICS BRIEF", ParagraphStyle('RHead', parent=meta_style_right, fontName='Helvetica-Bold', textColor=c_black))
        ],
        [
            Paragraph("Why Nearest-Gene Heuristics Fail in High-Throughput Perturbation Screens", ParagraphStyle('Sub', parent=meta_style, fontName='Helvetica-Bold', textColor=c_gray_dark)),
            Paragraph("Evaluating Distal Regulatory Targets Beyond 1D Proximity", meta_style_right)
        ]
    ]
    header_table = Table(header_data, colWidths=[370, 170])
    header_table.setStyle(TableStyle([
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1),
        ('TOPPADDING', (0,0), (-1,-1), 0),
        ('LEFTPADDING', (0,0), (-1,-1), 0),
        ('RIGHTPADDING', (0,0), (-1,-1), 0),
    ]))
    elements.append(header_table)
    elements.append(Spacer(1, 4))
    elements.append(HRFlowable(width="100%", thickness=1.2, color=c_line, spaceBefore=2, spaceAfter=8))

    # 2. The Core Problem Statement
    problem_text = """
    <b>The Core Failure Mode:</b> In primary CRISPR perturbation screens, over 60% of functional non-coding hits do not regulate the nearest transcription start site (TSS). Standard pipelines that default to 1D linear distance frequently assign regulatory elements to inactive or bystander promoters, leading teams to expend substantial target validation budgets on incorrect genes while missing the true distal driver.
    """
    elements.append(Paragraph(problem_text, body_text))
    elements.append(Spacer(1, 10))

    # 3. Canonical Case Study Box: FTO vs IRX3
    elements.append(Paragraph("<b>Case Study: Resolving Distal Target Selection at Canonical Loci</b>", main_heading))
    elements.append(Spacer(1, 4))

    case_text = """
    A prime example of this failure mode is the <b>FTO obesity-associated locus (rs1421085)</b>:
    <br/>
    • <b>The Linear Assignment (65 kb):</b> The non-coding variant resides in intron 1 of <i>FTO</i>. Linear proximity models assign <i>FTO</i> as the causal target. In metabolic target cells (preadipocytes), the <i>FTO</i> promoter is closed and transcriptionally uncoupled.
    <br/>
    • <b>The 3D Ground Truth (531 kb & 1.1 Mb):</b> Chromosome conformation capture (4C-seq / Capture-C) established that the functional enhancer loops across 531 kb to activate <i>IRX3</i>, with co-regulation extending to <i>IRX5</i> at 1.12 Mb within the same loop domain.
    <br/>
    • <b>Boundary Insulation:</b> Flanking genes located outside the loop anchor (such as <i>RBL2</i> at 317 kb) exhibit zero interaction due to strict domain boundary insulation.
    """
    case_table = Table([[Paragraph(case_text, body_text)]], colWidths=[540])
    case_table.setStyle(TableStyle([
        ('LINEABOVE', (0,0), (-1,-1), 0.75, c_line),
        ('LINEBELOW', (0,0), (-1,-1), 0.75, c_line),
        ('TOPPADDING', (0,0), (-1,-1), 6),
        ('BOTTOMPADDING', (0,0), (-1,-1), 6),
        ('LEFTPADDING', (0,0), (-1,-1), 8),
        ('RIGHTPADDING', (0,0), (-1,-1), 8),
    ]))
    elements.append(case_table)
    elements.append(Spacer(1, 12))

    # 4. Comparative Evaluation Table
    elements.append(Paragraph("<b>Comparison of Target Assignment Methodologies</b>", section_heading))
    elements.append(Spacer(1, 4))

    comp_data = [
        [Paragraph("<b>Evaluation Dimension</b>", table_header),
         Paragraph("<b>Linear Proximity (Nearest TSS)</b>", table_header),
         Paragraph("<b>Statistical eQTL Colocalization</b>", table_header),
         Paragraph("<b>3D Regulatory Domain Mapping</b>", table_header)],
        [Paragraph("<b>Distal Contacts (&gt;100 kb)</b>", table_cell_bold),
         Paragraph("Fails; systematically defaults to adjacent promoter regardless of conformation", table_cell),
         Paragraph("Constrained by tissue availability, sample size, and LD confounding", table_cell),
         Paragraph("<b>Accurately identifies long-range looping contacts to distal promoters</b>", table_cell_bold)],
        [Paragraph("<b>Boundary Insulation</b>", table_cell_bold),
         Paragraph("Blind to domain architecture and insulation boundaries", table_cell),
         Paragraph("Indirect; cannot establish physical boundary stops", table_cell),
         Paragraph("<b>Explicitly enforces domain insulation to eliminate bystander genes</b>", table_cell_bold)],
        [Paragraph("<b>Validation Risk</b>", table_cell_bold),
         Paragraph("High risk of pursuing false-positive proximal genes", table_cell),
         Paragraph("Moderate risk of false positives in shared LD blocks", table_cell),
         Paragraph("<b>Prioritizes verified physical interactors for downstream assay</b>", table_cell_bold)]
    ]
    t_comp = Table(comp_data, colWidths=[120, 140, 140, 140])
    t_comp.setStyle(TableStyle([
        ('LINEABOVE', (0,0), (-1,0), 1.0, c_black),
        ('LINEBELOW', (0,0), (-1,0), 0.5, c_black),
        ('LINEBELOW', (0,-1), (-1,-1), 1.0, c_black),
        ('TOPPADDING', (0,0), (-1,-1), 4.5),
        ('BOTTOMPADDING', (0,0), (-1,-1), 4.5),
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
    ]))
    elements.append(t_comp)
    elements.append(Spacer(1, 14))

    # 5. Strategic Takeaway & Contact Box
    takeaway_text = """
    <b>Strategic Question for Discovery Teams:</b> When screening thousands of non-coding perturbations, linear distance leaves over half of true causal targets uncharacterized. Incorporating 3D domain architecture and boundary insulation prior to experimental follow-up prevents misallocating validation budgets on bystander genes.
    <br/><br/>
    <b>Inquiries & Technical Collaboration:</b> Himanshu Bhandary &nbsp;|&nbsp; BDB-Genomics &nbsp;|&nbsp; <b>2032ushimanshu@gmail.com</b>
    """
    takeaway_table = Table([[Paragraph(takeaway_text, ParagraphStyle('TakeawayStyle', parent=body_text, fontSize=8, leading=11.5))]], colWidths=[540])
    takeaway_table.setStyle(TableStyle([
        ('LINEABOVE', (0,0), (-1,-1), 0.75, c_line),
        ('LINEBELOW', (0,0), (-1,-1), 0.75, c_line),
        ('TOPPADDING', (0,0), (-1,-1), 6),
        ('BOTTOMPADDING', (0,0), (-1,-1), 6),
        ('LEFTPADDING', (0,0), (-1,-1), 8),
        ('RIGHTPADDING', (0,0), (-1,-1), 8),
    ]))
    elements.append(takeaway_table)

    doc.build(elements)
    print(f"✅ Generated high-impact 1-page LinkedIn PDF for Tahoe: {pdf_path}")

if __name__ == "__main__":
    generate_pdf()
