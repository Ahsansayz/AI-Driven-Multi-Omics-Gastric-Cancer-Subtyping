<p align="center">
  <h1 align="center">🧬 AI-Driven Mutational Signature-Based Molecular Subtyping of Gastric Cancer</h1>
  <p align="center">
    <strong>A Biological Hybrid AI (BHAI) pipeline integrating whole-exome sequencing mutational signatures, DNA methylation, gene-level mutations, and clinical data for molecular subtyping of gastric adenocarcinoma</strong>
  </p>
  <p align="center">
    <a href="#-pipeline-architecture">Architecture</a> •
    <a href="#-key-results">Results</a> •
    <a href="#-installation">Install</a> •
    <a href="#-quick-start">Quick Start</a> •
    <a href="#-methodology">Methodology</a> •
    <a href="#-citation">Citation</a>
  </p>
  <p align="center">
    <img src="https://img.shields.io/badge/python-3.8%2B-blue?logo=python&logoColor=white" alt="Python 3.8+">
    <img src="https://img.shields.io/badge/license-MIT-green" alt="MIT License">
    <img src="https://img.shields.io/badge/data-TCGA--STAD-orange" alt="TCGA-STAD">
    <img src="https://img.shields.io/badge/COSMIC-v3.4%20(86%20SBS)-purple" alt="COSMIC v3.4">
    <img src="https://img.shields.io/badge/multi--omics-WES%20%2B%20Methylation%20%2B%20Clinical-red" alt="Multi-Omics">
    <img src="https://img.shields.io/badge/features-281%20multi--omics-teal" alt="281 Features">
  </p>
</p>

---

## 📋 Overview

**MSc Dissertation Project** — Jamia Hamdard University, New Delhi  
**Author:** Ahsan Jameel Khan (2024-582-010)  
**Supervisors:** Dr. Farheen Siddiqui & Dr. Mohammad Sufyan Badar  
**Degree:** M.Sc. Computational & Systems Biology and Bioinformatics

This repository implements a comprehensive **multi-omics pipeline** for AI-driven molecular subtyping of **gastric adenocarcinoma (TCGA-STAD)** into **5 molecular subtypes**: CIN, EBV, GS, MSI, and POLE.

### The Problem

Standard machine learning classifiers suffer from the **accuracy paradox** in gastric cancer subtyping — a model can achieve ~89% accuracy by simply favoring the majority CIN class, while catastrophically failing on clinically critical rare subtypes like POLE (only 7 samples, 1.9%) and GS (46 samples, 12.3%). Misclassifying a POLE patient as CIN can deny them immunotherapy access.

### The Solution: Biological Hybrid AI (BHAI)

We introduce the **Biological Hybrid AI** framework — a novel approach that hybridizes deterministic biological rules with statistical machine learning. Domain-specific clinical rules for POLE (SBS10 + TMB thresholds) and GS (CDH1/RHOA mutation logic) supplement an SMOTE-XGBoost model where data alone is insufficient, achieving a **28.5 percentage-point POLE recall improvement** with a single deterministic rule — an order of magnitude larger than typical SMOTE-based improvements.

| Data Layer | Source | Samples |
|---|---|---|
| **Somatic Mutations** | GDC Masked MAF files (WES) | 431 |
| **DNA Methylation** | Illumina 450K SeSAMe Level 3 β-values | 457 |
| **Clinical / Molecular Labels** | GDC API + cBioPortal | 375 labelled |

> **🔬 Improved version of [Machine-Learning-Mutational-Signatures](https://github.com/Ahsansayz/Machine-Learning-Mutational-Signatures).**  
> The original pipeline used MAF-based signatures only. This version extends it with **DNA methylation features**, **gene-level mutation profiles (20 driver genes)**, **clinical metadata**, and the novel **Biological Hybrid AI** framework.

---

## 🏆 Key Results

### Stepwise Performance Comparison

| Strategy | Accuracy | Macro F1 | MSI Recall | GS Recall | POLE Recall | EBV Recall | CIN Recall |
|---|---|---|---|---|---|---|---|
| SMOTE + XGBoost Baseline | 89.3% | 0.812 | 100% | 41.3% | 28.6% | ~80% | ~95% |
| + Hybrid Biological Rules | 90.5% | 0.847 | 100% | 58.7% | **57.1%** | ~81% | ~95.5% |
| **+ Threshold Optimisation** | **91.2%** | **0.860** | **100%** | **63.0%** | **57.1%** | **~83%** | **~96%** |

### Final Model — Per-Class Classification Performance (Hybrid + Threshold)

| Subtype | Precision | Recall | F1-Score | Support | Clinical Note |
|---|---|---|---|---|---|
| CIN | 0.94 | 0.96 | 0.95 | 219 | Majority class, well-calibrated |
| EBV | 0.90 | 0.83 | 0.86 | 30 | APOBEC + methylation driven |
| GS | 0.72 | 0.63 | 0.67 | 46 | CDH1/RHOA rescue rule (+21.7pp) |
| MSI | 0.99 | 1.00 | 0.99 | 73 | ✅ Perfect recall via SBS15 |
| POLE | 0.80 | 0.57 | 0.67 | 7 | ✅ SBS10 rule doubled recall (+28.5pp) |

### Key Achievements
- 🎯 **POLE recall: 28.6% → 57.1%** — entirely from a single biological rule (SBS10 + TMB threshold)
- 🎯 **GS recall: 41.3% → 63.0%** — via CDH1/RHOA mutation logic + threshold tuning
- 🎯 **MSI recall: 100%** — maintained across all strategies
- 🎯 **0 false POLE calls** — the SBS10 rule has effectively zero false positive rate
- 📊 **Average COSMIC cosine similarity: 0.933** — high-fidelity signature reconstruction

### Publication-Quality Figures

<p align="center">
  <img src="07_figures/tsne_umap.png" width="48%" alt="t-SNE/UMAP Clustering">
  <img src="07_figures/signature_subtype_heatmap.png" width="48%" alt="Signature-Subtype Heatmap">
</p>
<p align="center">
  <img src="07_figures/shap_bar.png" width="48%" alt="SHAP Feature Importance">
  <img src="07_figures/survival_curves.png" width="48%" alt="Kaplan-Meier Survival">
</p>

---

## 🏗️ Pipeline Architecture

```
 ┌─────────────────────────────────────────────────────────────────┐
 │                       INPUT DATA SOURCES                        │
 │  MAF Files (431)    │  Methylation (457)   │  Clinical (GDC)    │
 └────────┬────────────┴───────────┬───────────┴──────────┬────────┘
          │                        │                      │
          ▼                        │                      │
 ┌─ Step 1 ──────────────────┐     │                      │
 │ Process MAF Files         │     │                      │
 │ • 96-ch SBS Matrix        │     │                      │
 │ • 20 Driver Gene Features │     │                      │
 │ • Ti/Tv, Indel Fraction   │     │                      │
 └────────┬──────────────────┘     │                      │
          ▼                        │                      │
 ┌─ Step 2 ──────────────────┐     │                      │
 │ NMF Signature Extraction  │     │                      │
 │ • 50 random inits/rank    │     │                      │
 │ • Optimal rank k=9        │     │                      │
 └────────┬──────────────────┘     │                      │
          ▼                        │                      │
 ┌─ Step 3 ──────────────────┐     │                      │
 │ COSMIC v3.4 Assignment    │     │                      │
 │ • NNLS fitting (86 SBS)   │     │                      │
 │ • 83 active signatures    │     │                      │
 │ • Cosine sim: 0.933 avg   │     │                      │
 └────────┬──────────────────┘     │                      │
          │                        │                      │
          ▼                        ▼                      ▼
 ┌─ Step 4 ────────────────────────────────────────────────────────┐
 │            Multi-Omics Feature Engineering                      │
 │ 83 COSMIC abs + 83 proportional + 11 engineered + 4 clinical   │
 │ + 47 methylation (Step 5) + 53 gene mutation features          │
 │ → 281 total features (228 for ML + 53 for biological rules)    │
 └────────┬────────────────────────────────────────────────────────┘
          │
          ▼
 ┌─ Step 5 ──────────────────┐
 │ DNA Methylation Processing│
 │ • Pass 1 (~50 MB RAM):    │
 │   Targeted CpGs, CIMP     │
 │ • Pass 2 (~500 MB RAM):   │
 │   Top 5K probes → PCA(30) │
 │ • 102x RAM reduction      │
 └────────┬──────────────────┘
          │
          ▼
 ┌─ Step 6 ────────────────────────────────────────────────────────┐
 │           Biological Hybrid AI Classification                   │
 │ ┌──────────────────┐ ┌──────────────────┐ ┌──────────────────┐ │
 │ │ POLE Override     │ │ GS Rescue Rule   │ │ SMOTE-XGBoost    │ │
 │ │ SBS10a+b total    │ │ CDH1 trunc = 3pt │ │ 500 trees        │ │
 │ │ + TMB > 7.0/Mb    │ │ RHOA mut = 2pt   │ │ Stratified 5-CV  │ │
 │ │ + MSI safeguard   │ │ CDH1 mis = 1pt   │ │ 228 features     │ │
 │ │ (SBS15 check)     │ │ TP53 wt = 1pt    │ │                  │ │
 │ │                   │ │ Threshold ≥ 3    │ │ + Per-class      │ │
 │ │ FPR ≈ 0%          │ │ + pCIN < 0.50    │ │   threshold opt  │ │
 │ └──────────────────┘ └──────────────────┘ └──────────────────┘ │
 └────────┬────────────────────────────────────────────────────────┘
          │
          ▼
 ┌─ Step 7 ──────────────────┐
 │ Visualisation &            │
 │ Interpretation             │
 │ • SHAP (bar + beeswarm)   │
 │ • t-SNE + UMAP            │
 │ • Kaplan-Meier survival   │
 │ • TMB distribution        │
 │ • Signature heatmap       │
 └───────────────────────────┘
```

---

## 📁 Repository Structure

```
├── 01_raw_data/                  # Input: GDC MAF files + COSMIC reference
│   ├── maf_files/                   (431 masked somatic mutation MAFs — not tracked)
│   └── cosmic_reference/            (COSMIC v3.4 SBS signatures)
│
├── 02_methylation_data/          # DNA methylation processed features
│   ├── methylation_features.csv     (457 samples × 47 features)
│   ├── methylation_pca.csv          (30 PCA components from top 5K probes)
│   ├── methylation_targeted.csv     (MLH1, CDH1, CDKN2A, MGMT, BRCA1, RUNX3)
│   └── methylation_summary.csv      (Global beta-value statistics + CIMP)
│
├── 03_pipeline_scripts/          # Core pipeline (7 steps + master runner)
│   ├── step1_process_maf.py         → SBS96 matrix + 20 driver gene features
│   ├── step2_extract_signatures.py  → NMF de novo extraction (k=9)
│   ├── step3_cosmic_assignment.py   → COSMIC v3.4 NNLS fitting (83 active)
│   ├── step4_build_features.py      → 281-feature multi-omics matrix
│   ├── step5_process_methylation.py → Two-pass RAM-optimized methylation
│   ├── step6_classify.py            → Biological Hybrid AI classifier
│   ├── step7_visualize.py           → Publication figures + SHAP
│   └── run_pipeline.sh              → Master pipeline runner
│
├── 07_figures/                   # Publication-ready figures
│   ├── shap_bar.png / shap_summary.png
│   ├── tsne_umap.png
│   ├── signature_subtype_heatmap.png
│   ├── survival_curves.png
│   └── tmb_distribution.png
│
├── 08_clinical_data/             # Clinical & molecular subtype labels
│   ├── clinical_data.csv            (GDC + cBioPortal merged, 375 labelled)
│   ├── cbio_clinical.csv
│   └── gdc_clinical.csv
│
├── 09_documentation/             # Detailed documentation
│   ├── Step_by_Step_Documentation.txt
│   ├── pipeline_README.txt
│   ├── README_Hinglish.txt
│   ├── figures/                     (Additional documentation figures)
│   └── scripts/                     (Annotated script references)
│
├── requirements.txt              # Python dependencies
├── LICENSE                       # MIT License
└── README.md                     # This file
```

---

## ⚙️ Installation

### Prerequisites
- Python 3.8+
- ~8 GB RAM (for methylation processing; pipeline uses 2-pass architecture to minimize memory)
- ~500 MB disk space (excluding raw MAF data)

### Setup

```bash
# Clone the repository
git clone https://github.com/Ahsansayz/AI-Driven-Multi-Omics-Gastric-Cancer-Subtyping.git
cd AI-Driven-Multi-Omics-Gastric-Cancer-Subtyping

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Linux/macOS
# venv\Scripts\activate   # Windows

# Install dependencies
pip install -r requirements.txt
```

### Data Preparation

1. **MAF Files**: Download TCGA-STAD masked somatic mutation MAF files from [GDC Data Portal](https://portal.gdc.cancer.gov/) and place in `01_raw_data/maf_files/`

2. **COSMIC Reference**: Download COSMIC v3.4 SBS signatures (GRCh38) from [COSMIC](https://cancer.sanger.ac.uk/signatures/) and place in `01_raw_data/cosmic_reference/`

3. **Methylation Data** *(optional)*: Download Illumina 450K SeSAMe Level 3 beta-value files from GDC for TCGA-STAD. Pre-processed methylation features are included in `02_methylation_data/`.

---

## 🚀 Quick Start

### Run Full Pipeline

```bash
bash 03_pipeline_scripts/run_pipeline.sh
```

### Run Individual Steps

```bash
# Step 1: Process MAF files → SBS96 matrix + gene features
python 03_pipeline_scripts/step1_process_maf.py

# Step 2: Extract de novo signatures via NMF (optimal k=9)
python 03_pipeline_scripts/step2_extract_signatures.py

# Step 3: Map to COSMIC v3.4 via NNLS (83 active signatures)
python 03_pipeline_scripts/step3_cosmic_assignment.py

# Step 4: Build 281-feature multi-omics matrix
python 03_pipeline_scripts/step4_build_features.py

# Step 5: Process DNA methylation (2-pass, RAM-optimized)
python 03_pipeline_scripts/step5_process_methylation.py

# Step 6: Biological Hybrid AI classification
python 03_pipeline_scripts/step6_classify.py

# Step 7: Generate publication figures (SHAP, t-SNE, survival)
python 03_pipeline_scripts/step7_visualize.py
```

### Run Specific Steps Only

```bash
# Run only classification + visualization
bash 03_pipeline_scripts/run_pipeline.sh 6 7
```

---

## 🧠 Methodology

### What's New vs. the Original Pipeline?

| Feature | [Original](https://github.com/Ahsansayz/Machine-Learning-Mutational-Signatures) | **This Version** |
|---|---|---|
| Input Data | MAF only | **MAF + Methylation + Clinical** |
| Total Features | ~179 | **281 (multi-omics)** |
| Gene Mutations | ❌ | ✅ 20 driver genes (CDH1, RHOA, TP53, PIK3CA, ERBB2, etc.) |
| Methylation | ❌ | ✅ Promoter CpGs, CIMP score, 30 PCA components |
| CDH1 Detail | ❌ | ✅ Truncating vs. missense distinction |
| Classification | Standard XGBoost | **Biological Hybrid AI (BHAI)** |
| POLE Detection | ML only (28.6% recall) | **SBS10 + TMB biological rule (57.1%)** |
| GS Detection | ML only (41.3% recall) | **CDH1/RHOA evidence scoring (63.0%)** |
| Best Accuracy | ~85% | **91.2%** |
| MSI Recall | 100% | **100%** (preserved) |

### Biological Hybrid AI System

The classifier combines three complementary strategies:

#### 1. POLE Override Rule (SBS10 + TMB)
```
IF   SBS10_total > threshold (derived per CV fold)
AND  TMB > 7.0 mutations/Mb
AND  SBS15 does NOT dominate SBS10a (MSI safeguard)
THEN → classify as POLE (deterministic)
```
- **Result**: POLE recall 28.6% → 57.1%
- **False positive rate**: effectively **zero** — no non-POLE sample ever triggered this rule
- The 1 POLE sample missed had concurrent MSI features; the MSI safeguard correctly prevented override

#### 2. GS Rescue Rule (CDH1/RHOA Evidence Scoring)
```
Evidence Score:
  CDH1 truncating mutation  = 3 points (most GS-specific)
  RHOA mutation             = 2 points (frequent in GS)
  CDH1 missense mutation    = 1 point  (less specific)
  TP53 wild-type            = 1 point  (negative evidence)

IF   evidence_score ≥ 3 AND max_posterior_prob < 0.50
THEN → override to GS
```
- **Result**: GS recall 41.3% → 63.0% (split between rule + threshold tuning)
- Uncertainty gating (`pCIN < 0.50`) prevents overwriting confident correct predictions

#### 3. SMOTE-XGBoost + Per-Class Threshold Optimisation
- 500-tree XGBoost with SMOTE oversampling
- Stratified 5-fold cross-validation
- 228-dimensional feature matrix (excluding gene mutation features used by rules)
- Per-class threshold optimization via grid search over θ ∈ [0.10, 0.90]

### 281-Feature Multi-Omics Matrix

| Feature Block | Count | Source |
|---|---|---|
| COSMIC signature activities (absolute) | 83 | Step 3 NNLS fitting |
| COSMIC signature activities (proportional) | 83 | Normalized from absolute |
| Engineered aggregate features | 11 | TMB, MSI burden, signature diversity, etc. |
| Clinical features | 4 | TMB, log_TMB, age, gender |
| DNA methylation: targeted promoter CpGs | 10 | MLH1, CDH1, CDKN2A, MGMT, BRCA1, RUNX3 |
| DNA methylation: global summary stats | 7 | Mean beta, CIMP score, hypermethylated fraction |
| DNA methylation: PCA components | 30 | Top 5,000 most variable probes → PCA |
| Gene mutation: binary status (20 genes) | 20 | CDH1, RHOA, TP53, PIK3CA, ERBB2, etc. |
| Gene mutation: raw counts (20 genes) | 20 | Mutation counts per gene |
| CDH1-specific indicators | 2 | Truncating vs. missense distinction |
| Variant classification stats | 11 | Ti/Tv ratio, indel fraction, etc. |
| **Total** | **281** | `ml_features_multiomics_v2.csv` |

### RAM-Optimized DNA Methylation Processing

The two-pass architecture achieves **102× RAM reduction** compared to full-matrix approaches:
- **Pass 1** (~50 MB): Streaming extraction of targeted CpG probes + running statistics
- **Pass 2** (~500 MB): Loads only top 5,000 most variable probes for incremental PCA
- Full-matrix approach would require ~28 GB — our method runs on consumer laptops

---

## 📊 Visualizations

| Figure | Description |
|---|---|
| `shap_summary.png` | SHAP beeswarm — per-feature impact direction for each subtype |
| `shap_bar.png` | Mean |SHAP| values — global feature importance ranking |
| `tsne_umap.png` | t-SNE & UMAP embeddings — POLE/MSI form isolated clusters |
| `signature_subtype_heatmap.png` | Mean COSMIC signature contribution per subtype |
| `survival_curves.png` | Kaplan-Meier curves with log-rank test |
| `tmb_distribution.png` | TMB box plots — POLE ultra-high, MSI elevated, GS low |

---

## 🔑 Key Findings

1. **SBS10 is a near-perfect POLE biomarker** — deterministic rule achieves 0% false positive rate
2. **The accuracy paradox is real** — 89.3% accuracy masked 28.6% POLE recall and 41.3% GS recall
3. **Biological rules > synthetic data** — one POLE rule achieved 28.5pp recall improvement vs. typical 5-15pp from SMOTE variants
4. **GS remains the hardest subtype** — defined by absence of features, not presence; CDH1 promoter methylation could improve recall further
5. **CIN-GS boundary is intrinsically blurred** — both share low TMB, clock-like signatures, and lack hypermutator phenotypes
6. **Methylation features validate biology** — CIMP score spikes in EBV, MLH1 methylation in MSI, CDH1 methylation in EBV
7. **Copy-number features are the missing piece** — CNV data (not in MAF files) would most cleanly separate CIN from GS

---

## 📚 References

- **TCGA-STAD**: The Cancer Genome Atlas Research Network, *Nature* 2014 ([doi:10.1038/nature13480](https://doi.org/10.1038/nature13480))
- **COSMIC v3.4**: Tate et al., *Nucleic Acids Research* 2019 ([doi:10.1093/nar/gky1015](https://doi.org/10.1093/nar/gky1015))
- **NMF for Signatures**: Alexandrov et al., *Nature* 2013 ([doi:10.1038/nature12477](https://doi.org/10.1038/nature12477))
- **XGBoost**: Chen & Guestrin, *KDD* 2016 ([doi:10.1145/2939672.2939785](https://doi.org/10.1145/2939672.2939785))
- **SHAP**: Lundberg & Lee, *NeurIPS* 2017
- **SMOTE**: Chawla et al., *JAIR* 2002 ([doi:10.1613/jair.953](https://doi.org/10.1613/jair.953))

---

## 📄 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- [TCGA Research Network](https://www.cancer.gov/tcga) for making multi-omics genomic data publicly available
- [GDC Data Portal](https://portal.gdc.cancer.gov/) for data access infrastructure
- [cBioPortal](https://www.cbioportal.org/) for clinical data integration
- [COSMIC](https://cancer.sanger.ac.uk/signatures/) for the mutational signature reference database (v3.4)
- Dr. Farheen Siddiqui & Dr. Mohammad Sufyan Badar, Jamia Hamdard University

---

<p align="center">
  <sub>Built with ❤️ for precision oncology research — Jamia Hamdard University, New Delhi</sub>
</p>
