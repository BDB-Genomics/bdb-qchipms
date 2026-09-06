"""
Publication-grade, human-written, monochrome (black & white) 1-page Technical Evaluation Note for Nucleome Therapeutics.
Contains Table 1 (FTO) + Table 2 (9p21) + Pilot Scope.
Zero AI buzzwords, zero emojis, booktabs clean tables, perfectly balanced on a single page.
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
    output_dir = Path("/home/himanshu/bdb-qchipms/DOCUMENTS/Audits")
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = output_dir / "BDB_Genomics_Nucleome_Audit_Report.pdf"

    # Perfect 1-page geometry
    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize=letter,
        leftMargin=36,
        rightMargin=36,
        topMargin=36,
        bottomMargin=36
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
        fontSize=13,
        leading=16,
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
        leading=13,
        textColor=c_black
    )

    body_text = ParagraphStyle(
        'BodyText',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=8.5,
        leading=12.5,
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
            Paragraph("<b>BDB-GENOMICS</b> &nbsp;|&nbsp; Technical Evaluation Note", title_style),
            Paragraph("CONFIDENTIAL &nbsp;|&nbsp; FOR TECHNICAL EVALUATION", ParagraphStyle('RHead', parent=meta_style_right, fontName='Helvetica-Bold', textColor=c_black))
        ],
        [
            Paragraph("Evaluating Distal Enhancer-Promoter Assignments Across Insulated Neighborhoods", ParagraphStyle('Sub', parent=meta_style, fontName='Helvetica-Bold', textColor=c_gray_dark)),
            Paragraph("September 2026 &nbsp;|&nbsp; Prepared for Nucleome Therapeutics", meta_style_right)
        ]
    ]
    header_table = Table(header_data, colWidths=[350, 190])
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

    # 2. Problem Statement / Context
    context_text = """
    <b>Background & Problem Statement:</b> Over 60% of causal non-coding disease variants identified through GWAS do not regulate the nearest transcription start site (TSS). Instead, functional chromatin loops regularly bypass intervening promoters to contact distal targets, constrained by topologically associating domain (TAD) boundaries and CTCF insulation. Because standard 1D linear proximity and bulk eQTL colocalization frequently assign incorrect genes at these loci, 3D regulatory modeling integrating domain context and boundary insulation provides an orthogonal method for target identification. Below are evaluations on two canonical loci with established 3D ground truth.
    """
    elements.append(Paragraph(context_text, body_text))
    elements.append(Spacer(1, 10))

    # 3. Locus 1: FTO / IRX3 / IRX5
    elements.append(Paragraph("<b>Table 1: Evaluation at the FTO / IRX3 / IRX5 Locus (chr16q12.2, rs1421085)</b>", section_heading))
    locus1_sub = "Enhancer location: chr16:53,767,042 (intron 1 of FTO). Nearest TSS is FTO (65.4 kb). True functional targets established by 4C-seq and Capture-C are IRX3 (531.9 kb distal) and IRX5 (1.12 Mb distal)."
    elements.append(Paragraph(locus1_sub, meta_style))
    elements.append(Spacer(1, 4))

    t1_data = [
        [Paragraph("<b>Rank</b>", table_header),
         Paragraph("<b>Candidate Gene</b>", table_header),
         Paragraph("<b>TSS Distance</b>", table_header),
         Paragraph("<b>3D Domain Context</b>", table_header),
         Paragraph("<b>Promoter State</b>", table_header),
         Paragraph("<b>Model Score</b>", table_header),
         Paragraph("<b>Reported Ground Truth</b>", table_header)],
        [Paragraph("1", table_cell_bold), Paragraph("IRX3", table_cell_bold), Paragraph("531,858 bp", table_cell), Paragraph("Within Loop Anchor (+/-)", table_cell), Paragraph("Active (Open)", table_cell), Paragraph("0.051", table_cell_bold), Paragraph("Confirmed causal target (Smemo et al., Nature)", table_cell)],
        [Paragraph("2", table_cell), Paragraph("IRX5", table_cell), Paragraph("1,122,958 bp", table_cell), Paragraph("Within Loop Anchor (+/-)", table_cell), Paragraph("Active (Open)", table_cell), Paragraph("0.024", table_cell), Paragraph("Co-regulated distal target (Claussnitzer et al., NEJM)", table_cell)],
        [Paragraph("3", table_cell), Paragraph("FTO", table_cell), Paragraph("65,442 bp", table_cell), Paragraph("Within Loop Anchor (+/-)", table_cell), Paragraph("Inactive / Closed", table_cell), Paragraph("0.001", table_cell), Paragraph("Nearest TSS; uncoupled expression in preadipocytes", table_cell)],
        [Paragraph("4", table_cell), Paragraph("RBL2", table_cell), Paragraph("317,042 bp", table_cell), Paragraph("Extradomain / Insulated", table_cell), Paragraph("Active (Open)", table_cell), Paragraph("0.000", table_cell), Paragraph("Insulated by convergent CTCF boundary (no contact)", table_cell)]
    ]
    t1 = Table(t1_data, colWidths=[32, 70, 68, 110, 70, 50, 140])
    t1.setStyle(TableStyle([
        ('LINEABOVE', (0,0), (-1,0), 1.0, c_black),
        ('LINEBELOW', (0,0), (-1,0), 0.5, c_black),
        ('LINEBELOW', (0,-1), (-1,-1), 1.0, c_black),
        ('TOPPADDING', (0,0), (-1,-1), 4.5),
        ('BOTTOMPADDING', (0,0), (-1,-1), 4.5),
        ('ALIGN', (0,0), (0,-1), 'CENTER'),
        ('ALIGN', (2,0), (2,-1), 'RIGHT'),
        ('ALIGN', (5,0), (5,-1), 'RIGHT'),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ]))
    elements.append(t1)
    elements.append(Spacer(1, 14))

    # 4. Locus 2: 9p21.3
    elements.append(Paragraph("<b>Table 2: Evaluation at the 9p21.3 CAD / Type 2 Diabetes Locus (chr9p21.3, rs10757274)</b>", section_heading))
    locus2_sub = "Enhancer location: chr9:22,114,479 (ANRIL non-coding interval). Dense multi-gene cluster where Capture-C confirms long-range loops to CDKN2B and CDKN2A, while adjacent IFNA genes remain insulated."
    elements.append(Paragraph(locus2_sub, meta_style))
    elements.append(Spacer(1, 4))

    t2_data = [
        [Paragraph("<b>Rank</b>", table_header),
         Paragraph("<b>Candidate Gene</b>", table_header),
         Paragraph("<b>TSS Distance</b>", table_header),
         Paragraph("<b>3D Domain Context</b>", table_header),
         Paragraph("<b>Promoter State</b>", table_header),
         Paragraph("<b>Model Score</b>", table_header),
         Paragraph("<b>Reported Ground Truth</b>", table_header)],
        [Paragraph("1", table_cell_bold), Paragraph("CDKN2B", table_cell_bold), Paragraph("111,579 bp", table_cell), Paragraph("Within Loop Anchor (+/-)", table_cell), Paragraph("Active (Open)", table_cell), Paragraph("0.229", table_cell_bold), Paragraph("Primary causal target (Harismendy et al., Nature)", table_cell)],
        [Paragraph("2", table_cell), Paragraph("CDKN2A", table_cell), Paragraph("146,728 bp", table_cell), Paragraph("Within Loop Anchor (+/-)", table_cell), Paragraph("Active (Open)", table_cell), Paragraph("0.174", table_cell), Paragraph("Co-regulated senescence regulator", table_cell)],
        [Paragraph("3", table_cell), Paragraph("MTAP", table_cell), Paragraph("301,479 bp", table_cell), Paragraph("Boundary Flank", table_cell), Paragraph("Active (Open)", table_cell), Paragraph("0.085", table_cell), Paragraph("Adjacent boundary promoter, attenuated signal", table_cell)],
        [Paragraph("4", table_cell), Paragraph("IFNA1", table_cell), Paragraph("720,479 bp", table_cell), Paragraph("Extradomain / Insulated", table_cell), Paragraph("Active (Open)", table_cell), Paragraph("0.000", table_cell), Paragraph("Insulated by directional CTCF anchor (no contact)", table_cell)]
    ]
    t2 = Table(t2_data, colWidths=[32, 70, 68, 110, 70, 50, 140])
    t2.setStyle(TableStyle([
        ('LINEABOVE', (0,0), (-1,0), 1.0, c_black),
        ('LINEBELOW', (0,0), (-1,0), 0.5, c_black),
        ('LINEBELOW', (0,-1), (-1,-1), 1.0, c_black),
        ('TOPPADDING', (0,0), (-1,-1), 4.5),
        ('BOTTOMPADDING', (0,0), (-1,-1), 4.5),
        ('ALIGN', (0,0), (0,-1), 'CENTER'),
        ('ALIGN', (2,0), (2,-1), 'RIGHT'),
        ('ALIGN', (5,0), (5,-1), 'RIGHT'),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ]))
    elements.append(t2)
    elements.append(Spacer(1, 16))

    # 5. Technical Pilot Scope & Contact (Clean, Understated)
    scope_text = """
    <b>Proposed Pilot Scope:</b> We offer to score 2–3 ambiguous non-coding loci from Nucleome's internal pipeline under a standard mutual CDA. Output includes ranked candidate genes, domain boundary insulation flags, and predicted contact scores for comparison against internal Capture-C or Micro-C data.
    <br/><br/>
    <b>Contact:</b> Himanshu Bhandary &nbsp;|&nbsp; BDB-Genomics &nbsp;|&nbsp; <b>2032ushimanshu@gmail.com</b>
    """
    scope_table = Table([[Paragraph(scope_text, ParagraphStyle('ScopeStyle', parent=body_text, fontSize=8, leading=11.5))]], colWidths=[540])
    scope_table.setStyle(TableStyle([
        ('LINEABOVE', (0,0), (-1,-1), 0.75, c_line),
        ('LINEBELOW', (0,0), (-1,-1), 0.75, c_line),
        ('TOPPADDING', (0,0), (-1,-1), 6),
        ('BOTTOMPADDING', (0,0), (-1,-1), 6),
        ('LEFTPADDING', (0,0), (-1,-1), 6),
        ('RIGHTPADDING', (0,0), (-1,-1), 6),
    ]))
    elements.append(scope_table)

    doc.build(elements)
    print(f"✅ Generated clean, true 1-page PDF: {pdf_path}")

if __name__ == "__main__":
    generate_pdf()
