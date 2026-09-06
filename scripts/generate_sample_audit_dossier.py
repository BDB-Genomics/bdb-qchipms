"""
Script to generate an authentic, publication-grade 1-page Sample Screening Audit Dossier.
Outputs to DOCUMENTS/Audits/BDB_Genomics_Sample_Audit_Dossier.pdf.
Strictly zero em-dashes, fully auditable GRCh38/Ensembl coordinates and formulas.
"""

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
    pdf_path = output_dir / "BDB_Genomics_Sample_Audit_Dossier.pdf"

    # Exact 1-page fit on letter size (612 x 792 pt)
    # Available width = 612 - 20 - 20 = 572 pt
    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize=letter,
        leftMargin=20,
        rightMargin=20,
        topMargin=18,
        bottomMargin=18
    )

    styles = getSampleStyleSheet()

    # Color Palette
    primary_color = colors.HexColor("#0f172a") # Slate 900
    accent_blue = colors.HexColor("#1e40af")   # Blue 800
    teal_dark = colors.HexColor("#0f766e")     # Teal 700
    text_dark = colors.HexColor("#1e293b")     # Slate 800
    text_muted = colors.HexColor("#64748b")    # Slate 500
    bg_light = colors.HexColor("#f8fafc")      # Slate 50
    border_color = colors.HexColor("#cbd5e1")  # Slate 300
    green_bg = colors.HexColor("#f0fdf4")      # Green 50
    amber_bg = colors.HexColor("#fffbeb")      # Amber 50

    title_style = ParagraphStyle(
        'DocTitle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=12,
        leading=15,
        textColor=primary_color
    )

    subtitle_style = ParagraphStyle(
        'DocSubtitle',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=7.5,
        leading=10,
        textColor=accent_blue
    )

    meta_style = ParagraphStyle(
        'DocMeta',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=6.5,
        leading=8.5,
        textColor=text_muted
    )

    section_heading = ParagraphStyle(
        'SectionHead',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=7.8,
        leading=9.8,
        textColor=accent_blue
    )

    body_text = ParagraphStyle(
        'BodyText',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=6.5,
        leading=8.5,
        textColor=text_dark
    )

    th_style = ParagraphStyle(
        'THStyle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=5.8,
        leading=7.2,
        textColor=colors.white
    )

    th_center = ParagraphStyle(
        'THCenter',
        parent=th_style,
        alignment=1
    )

    td_style = ParagraphStyle(
        'TDStyle',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=5.6,
        leading=7.2,
        textColor=text_dark
    )

    td_center = ParagraphStyle(
        'TDCenter',
        parent=td_style,
        alignment=1
    )

    elements = []

    # 1. Header Table
    header_data = [
        [
            Paragraph("<b>BDB-GENOMICS</b> | FAST IN SILICO ABC REGULATORY ENGINE", title_style),
            Paragraph("<b>SAMPLE AUDIT DOSSIER</b><br/><font color='#64748b'>Commercial Deliverable Specification</font>", ParagraphStyle('RHead', parent=meta_style, alignment=2))
        ],
        [
            Paragraph("Screening Target Audit: Resolving Gene Desert Skips (>100 kb), CTCF Boundaries, and Multi-Promoter Competition", subtitle_style),
            Paragraph("Turnaround: 14 Business Days | Inquiries: <b>2032ushimanshu@gmail.com</b>", ParagraphStyle('RMeta', parent=meta_style, alignment=2))
        ]
    ]
    header_table = Table(header_data, colWidths=[402, 170])
    header_table.setStyle(TableStyle([
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1),
        ('TOPPADDING', (0,0), (-1,-1), 0),
        ('LEFTPADDING', (0,0), (-1,-1), 0),
        ('RIGHTPADDING', (0,0), (-1,-1), 0),
    ]))
    elements.append(header_table)
    elements.append(Spacer(1, 2))
    elements.append(HRFlowable(width="100%", thickness=1.2, color=accent_blue, spaceBefore=1, spaceAfter=2))

    # 2. Executive Scope & Objective Box
    callout_text = """
    <b>AUDIT OBJECTIVE:</b> CRISPR screening and GWAS hits frequently map to non-coding deserts where 1D distance heuristics (Nearest-TSS) misassign the functional gene. BDB-Genomics executes a 14-day computational audit using a fast In Silico Activity-by-Contact (ABC) framework (Fulco et al. 2019 Nat Genet) enhanced with parameterized CTCF boundary gating (Rao et al. 2014) and cell-type promoter competence (CAGE/FANTOM5) to rank causal targets without requiring expensive per-cell Hi-C assays.
    <br/><b>BENCHMARK GROUND TRUTH:</b> On held-out CRISPRi distance-failure loci (N=67 where 1D nearest-TSS has 0% accuracy), BDB-Genomics achieves <b>83.58% Top-3 accuracy (56/67)</b>, outperforming the Broad Institute ABC model (79.10%, 53/67). All coordinates anchored to Ensembl GRCh38 Release 113 / dbSNP Build 156.
    """
    callout_data = [[Paragraph(callout_text, ParagraphStyle('Callout', parent=body_text, fontSize=6.1, leading=8.0))]]
    callout_table = Table(callout_data, colWidths=[572])
    callout_table.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), colors.HexColor("#eff6ff")),
        ('BOX', (0,0), (-1,-1), 0.8, colors.HexColor("#bfdbfe")),
        ('TOPPADDING', (0,0), (-1,-1), 2.5),
        ('BOTTOMPADDING', (0,0), (-1,-1), 2.5),
        ('LEFTPADDING', (0,0), (-1,-1), 5),
        ('RIGHTPADDING', (0,0), (-1,-1), 5),
    ]))
    elements.append(callout_table)
    elements.append(Spacer(1, 3))

    # 3. Table 1: 8q24.21 Screen De-risking Matrix
    elements.append(Paragraph("<b>TABLE 1: 8q24.21 Oncology Desert Screen Audit (Fulco et al. 2016 Science e1-e4 + rs6983267 Enhancer)</b>", section_heading))
    locus1_sub = "<i>Context: Core 8q24.21 cluster in K562. 1D nearest-TSS assigns inactive lncRNAs (CASC8, CCAT2). BDB tracks active polymer contact to MYC (244-334 kb distal).</i>"
    elements.append(Paragraph(locus1_sub, meta_style))
    elements.append(Spacer(1, 1.5))

    # Col Widths: 572 total -> [46, 92, 70, 44, 68, 44, 40, 40, 128]
    t1_headers = [
        Paragraph("Hit ID", th_center),
        Paragraph("Genomic Coord (GRCh38)", th_style),
        Paragraph("1D Nearest TSS", th_style),
        Paragraph("1D Dist", th_center),
        Paragraph("BDB #1 Target", th_style),
        Paragraph("Loop Dist", th_center),
        Paragraph("P_contact", th_center),
        Paragraph("Share", th_center),
        Paragraph("Audit De-Risking Action & Rationale", th_style)
    ]

    t1_rows = [
        t1_headers,
        [
            Paragraph("<b>CR-8q-e1</b>", td_center),
            Paragraph("chr8:127465656-655", td_style),
            Paragraph("CASC8 (lncRNA)", td_style),
            Paragraph("15.5 kb", td_center),
            Paragraph("<b>MYC</b> (Oncogene)", td_style),
            Paragraph("268.8 kb", td_center),
            Paragraph("0.0183", td_center),
            Paragraph("67.5%", td_center),
            Paragraph("<font color='#166534'><b>[PROMOTE MYC / DEMOTE CASC8]</b></font> CASC8 promoter closed (CAGE 0.05 TPM). CRISPRi proves MYC proliferation loss.", td_style)
        ],
        [
            Paragraph("<b>CR-8q-e2</b>", td_center),
            Paragraph("chr8:127471556-055", td_style),
            Paragraph("CASC8 (lncRNA)", td_style),
            Paragraph("9.8 kb", td_center),
            Paragraph("<b>MYC</b> (Oncogene)", td_style),
            Paragraph("263.1 kb", td_center),
            Paragraph("0.0186", td_center),
            Paragraph("67.6%", td_center),
            Paragraph("<font color='#166534'><b>[PROMOTE MYC / DEMOTE CASC8]</b></font> Bypasses non-expressed CASC8; loops directly to active MYC core.", td_style)
        ],
        [
            Paragraph("<b>CR-8q-e3</b>", td_center),
            Paragraph("chr8:127482856-855", td_style),
            Paragraph("CASC8 (lncRNA)", td_style),
            Paragraph("1.7 kb", td_center),
            Paragraph("<b>MYC</b> (Oncogene)", td_style),
            Paragraph("251.6 kb", td_center),
            Paragraph("0.0195", td_center),
            Paragraph("67.8%", td_center),
            Paragraph("<font color='#166534'><b>[PROMOTE MYC / DEMOTE CASC8]</b></font> 1.7 kb proximal TSS is epigenetically inert; 3D contact drives MYC.", td_style)
        ],
        [
            Paragraph("<b>CR-8q-e4</b>", td_center),
            Paragraph("chr8:127490056-055", td_style),
            Paragraph("CASC8 (lncRNA)", td_style),
            Paragraph("8.9 kb", td_center),
            Paragraph("<b>MYC</b> (Oncogene)", td_style),
            Paragraph("244.4 kb", td_center),
            Paragraph("0.0200", td_center),
            Paragraph("67.9%", td_center),
            Paragraph("<font color='#166534'><b>[PROMOTE MYC / DEMOTE CASC8]</b></font> CRISPRi titration validates 78% MYC knockdown; CASC8 unaffected.", td_style)
        ],
        [
            Paragraph("<b>CR-8q-rs698</b>", td_center),
            Paragraph("chr8:127401060 (SNP)", td_style),
            Paragraph("CCAT2 (lncRNA)", td_style),
            Paragraph("661 bp", td_center),
            Paragraph("<b>MYC</b> (Oncogene)", td_style),
            Paragraph("334.4 kb", td_center),
            Paragraph("0.0147", td_center),
            Paragraph("66.7%", td_center),
            Paragraph("<font color='#166534'><b>[PROMOTE MYC / DEMOTE CCAT2]</b></font> Sits in CCAT2; CCAT2 promoter silenced in K562. 3D contact drives MYC (Pomerantz 2009 3C).", td_style)
        ]
    ]

    t1 = Table(t1_rows, colWidths=[46, 92, 70, 44, 68, 44, 40, 40, 128])
    t1.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), primary_color),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1.2),
        ('TOPPADDING', (0,0), (-1,-1), 1.2),
        ('LEFTPADDING', (0,0), (-1,-1), 2.5),
        ('RIGHTPADDING', (0,0), (-1,-1), 2.5),
        ('BACKGROUND', (0,1), (-1,5), green_bg),
        ('GRID', (0,0), (-1,-1), 0.5, border_color),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ]))
    elements.append(t1)
    elements.append(Spacer(1, 3))

    # 4. Table 2: Multi-Promoter Tug-of-War Matrix
    elements.append(Paragraph("<b>TABLE 2: Multi-Promoter Kinetic Competition & Boundary Insulation (GWAS & 3C/4C Benchmarks)</b>", section_heading))
    locus2_sub = "<i>Context: Enhancers frequently bypass proximal promoters due to closed chromatin, or fail to contact neighboring genes due to CTCF insulation.</i>"
    elements.append(Paragraph(locus2_sub, meta_style))
    elements.append(Spacer(1, 1.5))

    t2_headers = [
        Paragraph("Locus & Hit", th_style),
        Paragraph("Genomic Coord", th_style),
        Paragraph("1D Nearest TSS", th_style),
        Paragraph("1D Dist", th_center),
        Paragraph("BDB Rank #1 Target", th_style),
        Paragraph("Target Dist", th_center),
        Paragraph("Comp. Share", th_center),
        Paragraph("Tug-of-War Mechanism & Experimental Validation", th_style)
    ]

    t2_rows = [
        t2_headers,
        [
            Paragraph("<b>FTO / IRX3</b><br/><font color='#64748b'>Obesity / Browning</font>", td_style),
            Paragraph("chr16:53767042<br/>(rs1421085)", td_style),
            Paragraph("FTO (Intron 1)", td_style),
            Paragraph("65.4 kb", td_center),
            Paragraph("<b>IRX3</b> (Homeobox)", td_style),
            Paragraph("520.4 kb", td_center),
            Paragraph("70.8%<br/>(IRX5: 29.2%)", td_center),
            Paragraph("<font color='#166534'><b>[PROMOTE IRX3 / DEMOTE FTO]</b></font> FTO promoter closed in preadipocytes. 4C/CRISPR confirms long-range loop to IRX3/IRX5 (Claussnitzer 2015 NEJM).", td_style)
        ],
        [
            Paragraph("<b>9p21.3 / CAD</b><br/><font color='#64748b'>Coronary Disease</font>", td_style),
            Paragraph("chr9:22124478<br/>(rs10757278)", td_style),
            Paragraph("CDKN2B (TSS)<br/>(in ANRIL body)", td_style),
            Paragraph("115.2 kb<br/>(130 kb ANRIL)", td_center),
            Paragraph("<b>CDKN2B / CDKN2A</b>", td_style),
            Paragraph("115.2 kb /<br/>129.2 kb", td_center),
            Paragraph("CDKN2B: 36.1%<br/>CDKN2A: 28.1%", td_center),
            Paragraph("<font color='#166534'><b>[MULTI-PROMOTER TARGET]</b></font> Harismendy 2011 Nature 3C/4C in HUVECs + IFN-gamma proves physical interaction with CDKN2A/B and MTAP; IFNA cluster insulated.", td_style)
        ],
        [
            Paragraph("<b>TCF7L2 / T2D</b><br/><font color='#64748b'>Type 2 Diabetes</font>", td_style),
            Paragraph("chr10:112998590<br/>(rs7903146)", td_style),
            Paragraph("TCF7L2 (Intron 4)", td_style),
            Paragraph("48.6 kb", td_center),
            Paragraph("<b>TCF7L2</b> (Islet TF)", td_style),
            Paragraph("48.6 kb", td_center),
            Paragraph("100.0%<br/>(VTI1A: 0%)", td_center),
            Paragraph("<font color='#166534'><b>[CONCORDANT / INSULATED]</b></font> 3C in human islets confirms intronic loop to TCF7L2 (Gaulton 2010 Nat Genet / Greenwald 2019 Nat Commun). VTI1A (551 kb) is insulated across CTCF boundary.", td_style)
        ]
    ]

    t2 = Table(t2_rows, colWidths=[68, 70, 72, 44, 76, 46, 52, 144])
    t2.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#1e293b")),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1.2),
        ('TOPPADDING', (0,0), (-1,-1), 1.2),
        ('LEFTPADDING', (0,0), (-1,-1), 2.5),
        ('RIGHTPADDING', (0,0), (-1,-1), 2.5),
        ('GRID', (0,0), (-1,-1), 0.5, border_color),
        ('BACKGROUND', (0,1), (-1,-1), bg_light),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ]))
    elements.append(t2)
    elements.append(Spacer(1, 3))

    # 5. Full Per-Row Candidate Triage Breakdown (Side-by-Side Loci)
    elements.append(Paragraph("<b>AUDIT BREAKDOWN: Full Candidate Triage Matrices (8q24.21 vs 9p21.3 CAD Competition Fields)</b>", section_heading))
    locus3_sub = "<i>Candidate field deconvolution: 8q24.21 filters inactive lncRNAs to promote MYC (67.5%). 9p21.3 allocates competitive share across active cluster (CDKN2B 36.1%, CDKN2A 28.1%, ANRIL 25.0%, MTAP 10.8%).</i>"
    elements.append(Paragraph(locus3_sub, meta_style))
    elements.append(Spacer(1, 1.5))

    # Table 3A: 8q24.21 Screen Candidate Triage (Left Table, width = 283 pt)
    t3a_headers = [
        Paragraph("Rank", th_center),
        Paragraph("Gene (8q)", th_style),
        Paragraph("Dist", th_center),
        Paragraph("CAGE", th_center),
        Paragraph("Score", th_center),
        Paragraph("Share", th_center),
        Paragraph("Triage Status", th_style)
    ]
    t3a_rows = [
        t3a_headers,
        [
            Paragraph("<b>#1</b>", td_center),
            Paragraph("<b>MYC</b>", td_style),
            Paragraph("268.8k", td_center),
            Paragraph("184.2", td_center),
            Paragraph("0.0417", td_center),
            Paragraph("<b>67.5%</b>", td_center),
            Paragraph("<font color='#166534'><b>TARGET</b></font> (Open/Loop)", td_style)
        ],
        [
            Paragraph("#2", td_center),
            Paragraph("PVT1", td_style),
            Paragraph("327.9k", td_center),
            Paragraph("12.8", td_center),
            Paragraph("0.0201", td_center),
            Paragraph("32.5%", td_center),
            Paragraph("<font color='#b45309'>DECOY</font> (Active prom.)", td_style)
        ],
        [
            Paragraph("#3", td_center),
            Paragraph("CASC8", td_style),
            Paragraph("15.5k", td_center),
            Paragraph("0.05", td_center),
            Paragraph("0.0012", td_center),
            Paragraph("0.0%", td_center),
            Paragraph("<font color='#64748b'>FILTERED</font> (Closed TSS)", td_style)
        ],
        [
            Paragraph("#4", td_center),
            Paragraph("CCAT2", td_style),
            Paragraph("66.3k", td_center),
            Paragraph("0.00", td_center),
            Paragraph("0.0004", td_center),
            Paragraph("0.0%", td_center),
            Paragraph("<font color='#64748b'>FILTERED</font> (Closed TSS)", td_style)
        ],
        [
            Paragraph("#5", td_center),
            Paragraph("POU5F1B", td_style),
            Paragraph("144.5k", td_center),
            Paragraph("0.20", td_center),
            Paragraph("0.0002", td_center),
            Paragraph("0.0%", td_center),
            Paragraph("<font color='#64748b'>FILTERED</font> (Closed prom.)", td_style)
        ],
        [
            Paragraph("#6", td_center),
            Paragraph("LRATD2", td_style),
            Paragraph("908.2k", td_center),
            Paragraph("0.00", td_center),
            Paragraph("0.0000", td_center),
            Paragraph("0.0%", td_center),
            Paragraph("<font color='#64748b'>FILTERED</font> (Insulated)", td_style)
        ]
    ]
    t3a = Table(t3a_rows, colWidths=[20, 48, 36, 30, 36, 32, 81])
    t3a.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#334155")),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1.0),
        ('TOPPADDING', (0,0), (-1,-1), 1.0),
        ('LEFTPADDING', (0,0), (-1,-1), 2.0),
        ('RIGHTPADDING', (0,0), (-1,-1), 2.0),
        ('GRID', (0,0), (-1,-1), 0.5, border_color),
        ('BACKGROUND', (0,1), (-1,1), green_bg),
        ('BACKGROUND', (0,2), (-1,2), amber_bg),
        ('BACKGROUND', (0,3), (-1,-1), bg_light),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ]))

    # Table 3B: 9p21.3 CAD Multi-Promoter Candidate Triage (Right Table, width = 283 pt)
    t3b_headers = [
        Paragraph("Rank", th_center),
        Paragraph("Gene (9p)", th_style),
        Paragraph("Dist", th_center),
        Paragraph("CAGE", th_center),
        Paragraph("Score", th_center),
        Paragraph("Share", th_center),
        Paragraph("Triage Status", th_style)
    ]
    t3b_rows = [
        t3b_headers,
        [
            Paragraph("<b>#1</b>", td_center),
            Paragraph("<b>CDKN2B</b>", td_style),
            Paragraph("115.2k", td_center),
            Paragraph("38.2", td_center),
            Paragraph("0.0697", td_center),
            Paragraph("<b>36.1%</b>", td_center),
            Paragraph("<font color='#166534'><b>TARGET</b></font> (Primary prom.)", td_style)
        ],
        [
            Paragraph("#2", td_center),
            Paragraph("CDKN2A", td_style),
            Paragraph("129.2k", td_center),
            Paragraph("19.4", td_center),
            Paragraph("0.0542", td_center),
            Paragraph("28.1%", td_center),
            Paragraph("<font color='#166534'><b>CO-TARGET</b></font> (3C active)", td_style)
        ],
        [
            Paragraph("#3", td_center),
            Paragraph("ANRIL", td_style),
            Paragraph("130.3k", td_center),
            Paragraph("11.2", td_center),
            Paragraph("0.0482", td_center),
            Paragraph("25.0%", td_center),
            Paragraph("<font color='#b45309'>CO-TARGET</font> (lncRNA)", td_style)
        ],
        [
            Paragraph("#4", td_center),
            Paragraph("MTAP", td_style),
            Paragraph("321.9k", td_center),
            Paragraph("14.1", td_center),
            Paragraph("0.0209", td_center),
            Paragraph("10.8%", td_center),
            Paragraph("<font color='#b45309'>CO-TARGET</font> (3C loop)", td_style)
        ],
        [
            Paragraph("#5", td_center),
            Paragraph("IFNA1", td_style),
            Paragraph("684.0k", td_center),
            Paragraph("0.00", td_center),
            Paragraph("0.0000", td_center),
            Paragraph("0.0%", td_center),
            Paragraph("<font color='#64748b'>FILTERED</font> (Insulated)", td_style)
        ],
        [
            Paragraph("-", td_center),
            Paragraph("<i>Domain</i>", td_style),
            Paragraph("21.7-22M", td_center),
            Paragraph("HUVEC", td_center),
            Paragraph("Sum", td_center),
            Paragraph("100%", td_center),
            Paragraph("Harismendy 2011 3C", td_style)
        ]
    ]
    t3b = Table(t3b_rows, colWidths=[20, 48, 36, 30, 36, 32, 81])
    t3b.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#334155")),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1.0),
        ('TOPPADDING', (0,0), (-1,-1), 1.0),
        ('LEFTPADDING', (0,0), (-1,-1), 2.0),
        ('RIGHTPADDING', (0,0), (-1,-1), 2.0),
        ('GRID', (0,0), (-1,-1), 0.5, border_color),
        ('BACKGROUND', (0,1), (-1,1), green_bg),
        ('BACKGROUND', (0,2), (-1,4), amber_bg),
        ('BACKGROUND', (0,5), (-1,5), bg_light),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ]))

    # Side-by-Side Container Table (572 pt wide)
    t3_container = Table([[t3a, t3b]], colWidths=[283, 283])
    t3_container.setStyle(TableStyle([
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('BOTTOMPADDING', (0,0), (-1,-1), 0),
        ('TOPPADDING', (0,0), (-1,-1), 0),
        ('LEFTPADDING', (0,0), (-1,-1), 0),
        ('RIGHTPADDING', (0,0), (-1,-1), 0),
    ]))
    elements.append(t3_container)
    elements.append(Spacer(1, 3))

    # 6. Pilot Specifications & Scoring Formula Footnote Box
    footer_text = """
    <b>SCORING FORMULA & METHODS:</b> In Silico ABC formulation: Score = P_contact(d) * F_TAD * A_promoter. Polymer contact prior: P_contact = [5000 / (dist_bp + 5000)]^gamma (gamma=1.0 inside TAD, 1.75 across boundary). Boundary factor: F_TAD = 1.0 (in convergent CTCF loop), 0.08 (insulated across boundary). Promoter competence: A_promoter = log10(CAGE_TPM + 9) if active (CAGE >= 1 TPM / ATAC peak), else 0.005 (silenced). Competitive share = Score / sum(Active_Scores_in_Domain).
    <br/><b>14-DAY PILOT AUDIT SPECIFICATION:</b> 20 to 50 candidate screening hits audited per engagement. Deliverables include: (1) Target De-Risking TSV with full candidate rankings, (2) Promoter Competition & Decoy Matrix, (3) Vector-mapped CTCF boundary profiles. Turnaround: 14 Business Days. Contact: <b>2032ushimanshu@gmail.com</b>
    """
    footer_data = [[Paragraph(footer_text, ParagraphStyle('Footer', parent=body_text, fontSize=5.5, leading=7.2))]]
    footer_table = Table(footer_data, colWidths=[572])
    footer_table.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), colors.HexColor("#f1f5f9")),
        ('BOX', (0,0), (-1,-1), 0.8, colors.HexColor("#cbd5e1")),
        ('TOPPADDING', (0,0), (-1,-1), 2.0),
        ('BOTTOMPADDING', (0,0), (-1,-1), 2.0),
        ('LEFTPADDING', (0,0), (-1,-1), 5),
        ('RIGHTPADDING', (0,0), (-1,-1), 5),
    ]))
    elements.append(footer_table)

    doc.build(elements)
    print(f"✅ Publication-grade 1-page Sample Audit Dossier generated: {pdf_path}")

if __name__ == "__main__":
    generate_pdf()
