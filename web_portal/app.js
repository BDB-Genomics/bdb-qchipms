// Locus Benchmark Database
const locusData = {
    tcf7l2: {
        id: "T2D_LOCUS_001",
        rsid: "rs7903146",
        chrom: "chr10",
        pos: "112,998,590",
        pip: "PIP = 0.998 (SuSiE Single Causal Variant)",
        tissue: "Pancreatic Islet (&beta;-cell)",
        mechanism: "FOXA2 Pioneer Motif Creation (T-risk allele)",
        helicalPhasing: "10.5 bp exact stereospecific co-factor distance",
        targetGene: "ACSL5 / Distal TCF7L2 Promoters (150 kb 3D loop)",
        legacyFailure: "Assumed proximal TCF7L2 promoter was primary driver without enhancer looping evidence.",
        causalProof: "T-allele creates high-affinity FOXA2 binding site, opening chromatin to recruit SWI/SNF complex and establishing an ABC 3D loop to downregulate insulin secretion."
    },
    mtnr1b: {
        id: "T2D_LOCUS_002",
        rsid: "rs10830963",
        chrom: "chr11",
        pos: "92,708,710",
        pip: "PIP = 0.984 (SuSiE 95% CS1)",
        tissue: "Human Islet & Hypothalamus",
        mechanism: "NEUROD1 / FOXA2 Disruption (G-risk allele)",
        helicalPhasing: "21.0 bp dual-turn stereospecific phase",
        targetGene: "MTNR1B Circadian Promoters",
        legacyFailure: "Discarded as non-functional intronic noise due to weak whole-blood eQTL expression.",
        causalProof: "G-allele drives nocturnal overexpression of MTNR1B specifically in &beta;-cells, blunting morning glucose-stimulated insulin secretion."
    },
    slc30a8: {
        id: "T2D_LOCUS_003",
        rsid: "rs13266634",
        chrom: "chr8",
        pos: "117,172,514",
        pip: "PIP = 0.991 (SuSiE 95% CS1)",
        tissue: "Pancreatic Endocrine Granules",
        mechanism: "PDX1 / MAFA Co-occupancy Enhancer Remodeling",
        helicalPhasing: "10.5 bp pitch alignment with adjacent CTCF boundary",
        targetGene: "SLC30A8 / Zinc Transporter 8 Operon",
        legacyFailure: "Failed to map regulatory insulator boundary using standard linear genomic distance.",
        causalProof: "Disrupts zinc crystallization in insulin secretory granules through stereospecific loss of PDX1 master transcriptional drive."
    },
    hnf4a: {
        id: "T2D_LOCUS_004",
        rsid: "rs1800574",
        chrom: "chr20",
        pos: "44,414,890",
        pip: "PIP = 0.962 (SuSiE 95% CS1)",
        tissue: "Liver Hepatocyte & Pancreatic Islet",
        mechanism: "HNF1A / GATA4 Co-factor Complex Failure",
        helicalPhasing: "31.5 bp triple-turn stereospecific phase",
        targetGene: "HNF4A P2 Promoter (&beta;-cell isoform)",
        legacyFailure: "Assumed liver P1 promoter was therapeutic driver; clinical trials failed on liver targets.",
        causalProof: "Variant specifically impairs islet P2 isoform transcription, proving &beta;-cell dysfunction is the causal driver."
    }
};

function loadLocus(key) {
    const data = locusData[key];
    if (!data) return;

    // Update active tab styling
    document.querySelectorAll('.locus-btn').forEach(btn => btn.classList.remove('active'));
    if (event && event.currentTarget) {
        event.currentTarget.classList.add('active');
    }

    const container = document.getElementById('locusDetailView');
    container.innerHTML = `
        <div class="locus-summary-card">
            <div class="summary-header">
                <span class="summary-coords">${data.chrom}:${data.pos} &bull; ${data.rsid}</span>
                <span class="pip-badge">${data.pip}</span>
            </div>
            <div class="summary-details">
                <div class="detail-row"><span>Tissue:</span> <strong>${data.tissue}</strong></div>
                <div class="detail-row"><span>Causal Switch:</span> <strong>${data.mechanism}</strong></div>
                <div class="detail-row"><span>Stereospecific Phase:</span> <strong>${data.helicalPhasing}</strong></div>
                <div class="detail-row"><span>Validated Target:</span> <strong>${data.targetGene}</strong></div>
            </div>
        </div>

        <div class="comparison-box-grid">
            <div class="comp-box box-old">
                <span class="box-label label-red">Legacy Nearest-Gene Prediction</span>
                <p>${data.legacyFailure}</p>
            </div>
            <div class="comp-box box-new">
                <span class="box-label label-green">BDB-Genomics Causal Proof</span>
                <p>${data.causalProof}</p>
            </div>
        </div>
    `;
}

// Modal Handlers
function openModal() {
    document.getElementById('pilotModal').classList.add('active');
}

function closeModal() {
    document.getElementById('pilotModal').classList.remove('active');
}

function handleModalOverlay(e) {
    if (e.target.id === 'pilotModal') {
        closeModal();
    }
}

// Form Submission
function handleFormSubmit(e) {
    e.preventDefault();
    const payload = {
        name: document.getElementById('leadName').value,
        email: document.getElementById('leadEmail').value,
        organization: document.getElementById('leadOrg').value,
        collabType: document.getElementById('leadType').value,
        notes: document.getElementById('leadNotes').value,
        timestamp: new Date().toISOString()
    };

    // Save lead locally
    const leads = JSON.parse(localStorage.getItem('bdb_pilot_leads') || '[]');
    leads.push(payload);
    localStorage.setItem('bdb_pilot_leads', JSON.stringify(leads));

    document.getElementById('evaluationForm').style.display = 'none';
    document.getElementById('modalSuccessState').style.display = 'block';

    setTimeout(() => {
        closeModal();
        document.getElementById('evaluationForm').reset();
        document.getElementById('evaluationForm').style.display = 'block';
        document.getElementById('modalSuccessState').style.display = 'none';
    }, 2800);
}

// Initialize on page load
document.addEventListener('DOMContentLoaded', () => {
    loadLocus('tcf7l2');
});
