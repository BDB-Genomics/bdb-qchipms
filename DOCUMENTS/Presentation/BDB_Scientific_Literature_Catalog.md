# BDB-Genomics — Curated Scientific Literature & Benchmark Repository

A living, annotated catalog of peer-reviewed literature, preprints, benchmark datasets, and open atlases shaping multi-omics, chromatin biophysics, single-cell epigenomics, and regulatory syntax.

---

## 1. Single-Cell Epigenomics, Multi-Omics & Regulatory Syntax

| Title | Authors / Lab | Journal / Year | Key Contribution / Resource | Relevance to BDB-Genomics | Link |
|---|---|:---:|---|---|:---:|
| **Multiomics and deep learning dissect regulatory syntax in human development** | Greenleaf Lab (Stanford) | *Nature* (2024) | High-resolution single-cell multiomic atlas (*HDMA*) decoding how TF motif spacing, order, and orientation dictate developmental gene expression. | Direct validation of BDB's Scale 2 ($10.5\text{ bp}$ helical rotational phasing $\theta$) governing TF syntax. | [Nature Link](https://doi.org/10.1038/s41586-024-07110-x) / [GitHub: GreenleafLab/HDMA](https://github.com/GreenleafLab/HDMA) |
| **An expanded codebook of human transcription factor DNA-binding specificity** | Jolma, Laverty, Hughes et al. (Toronto) | *Nature* (Aug 2026) | Comprehensive experimental determination of DNA-binding specificity for 332 previously uncharacterized human TFs via HT-SELEX, ChIP-seq, and PBMs. | Massively expands BDB's Scale 1 & Scale 2 TF motif dictionary for non-coding regulatory variant scoring. | [Nature Link](https://doi.org/10.1038/s41586-026-10798-9) |
| **Human nuclear ribomes reveal DNA-embedded ribonucleotides as epigenetic modulators** | Kundnani, Storici et al. (Georgia Tech) | *Cell* (Aug 2026) | Maps ~1 million DNA-embedded ribonucleotides (rNMPs) across the human genome; proves rNMPs modulate DNA supercoiling at active TSSs. | Validates BDB's polymer mechanics & torsional strain modeling linking DNA architecture to active transcription. | [Cell Link](https://doi.org/10.1016/j.cell.2026.07.031) / [demo.glass.bio](https://demo.glass.bio/datasets) |
| **Mitochondrial metabolism and epigenetic crosstalk drive SASP** | Martini, H., et al. | *Nature* (2026) | Demonstrates how mitochondrial metabolic flux (acetyl-CoA / SAM) directly alters histone acetylation (H3K27ac) and chromatin opening at senescence-associated secretory phenotype (SASP) loci. | Validates BDB's dynamic epigenomic activity modeling ($A = \sqrt{\text{ATAC} \times \text{H3K27ac}}$) across metabolic and aging stress states. | [Nature 656, 980–992 (2026)](https://doi.org/10.1038/nature) |
| **Ultrafast and reference-free sequence discovery in single-cell data (Malva)** | Awasthi et al. | *Nature* (2026) | Reference-free single-cell sequence discovery platform bypassing rigid alignment references for ultrafast variation calling. | Validates BDB's fast C-accelerated configuration and hardware-tuned algorithms over heavy legacy pipelines. | [Nature Link](https://doi.org/10.1038/nature) |
| **A large-scale single-cell CRISPRi screen defines the grammar of human enhancers** | Gasperini et al. (Shendure Lab) | *Cell* (2019) | 470 high-confidence enhancer–gene pairs mapped via pooled single-cell CRISPRi perturbations and scRNA-seq readouts. | Core independent gold-standard benchmark cohort ($n = 470$) for BDB-Genomics. | [Cell Link](https://doi.org/10.1016/j.cell.2018.11.029) |
| **Activity-by-Contact model of enhancer–promoter regulation from thousands of CRISPR perturbations** | Fulco et al. (Engreitz / Lander Labs) | *Nature Genetics* (2019) | Mathematical formulation of the Activity-by-Contact ($ABC$) model; CRISPRi-FlowFISH benchmark across 46 distal loci (Table S3a). | Primary baseline model ($ABC$) and calibration dataset ($n = 46$) for BDB biophysical comparisons. | [Nat Genet Link](https://doi.org/10.1038/s41588-019-0538-0) |

---

## 2. Protein Language Models, Biophysics & Fitness Landscapes

| Title | Authors / Lab | Journal / Year | Key Contribution / Resource | Relevance to BDB-Genomics | Link |
|---|---|:---:|---|---|:---:|
| **Protein language models are overly constrained by covariation** | Wilke Lab (UT Austin) | *Preprint / Blog* (2024/2026) | Demonstrates that sequence-only pLMs get trapped in conservative evolutionary fitness landscapes, struggling to propose novel functional mutations. | Strong conceptual justification for why BDB uses fundamental biophysical forces rather than pure sequence covariation. | [Blog / Paper](https://blog.genesmindsmachines.com) |

---

## 3. Epigenomics, Methylation & Regional Association (MRWAS)

| Title | Authors / Lab | Journal / Year | Key Contribution / Resource | Relevance to BDB-Genomics | Link |
|---|---|:---:|---|---|:---:|
| **Methylated Region-Wide Association Studies (MRWAS)** | Sergey Lobanov et al. (NeoGenomics) | *Patent / Methodology* (WO 2026/161559) | Fast rolling regional aggregation of CpG methylation signals to detect true Differentially Methylated Regions (DMRs) without prior gene annotation. | Mirrors BDB's domain-level biophysical aggregation over 1D point tests. | [Patent WO 2026/161559](https://patentscope.wipo.int) |

---

## 4. Targeted Therapeutics, Surfaceome & Targeted Delivery

| Title | Authors / Initiative | Organization / Year | Key Contribution / Resource | Relevance to BDB-Genomics | Link |
|---|---|:---:|---|---|:---:|
| **The Deliverome Project: An Open Atlas of the Human Surfaceome** | Deliverome Bio (Astera / Renaissance / Thermo Fisher) | *Focused Research Organization (FRO)* (2024–2026) | Open atlas mapping receptor abundance (mass spec), internalization (CRISPR screens), and intracellular routing to enable targeted delivery beyond the liver. | Potential prospective wet-lab validation partner for BDB's target prioritization framework. | [deliverome.org](https://deliverome.org/) |
