<p align="center">
  <h1 align="center">🧬 Multi-Omics Mutational Signature Analysis of Gastric Cancer</h1>
  <p align="center">
    <strong>An end-to-end bioinformatics & machine learning pipeline for molecular subtyping of Gastric Cancer using Whole Exome Sequencing, DNA Methylation, and Clinical Data</strong>
  </p>
  <p align="center">
    <a href="#-pipeline-architecture">Architecture</a> •
    <a href="#-key-results">Results</a> •
    <a href="#-installation">Install</a> •
    <a href="#-quick-start">Quick Start</a> •
    <a href="#-citation">Citation</a>
  </p>
  <p align="center">
    <img src="https://img.shields.io/badge/python-3.8%2B-blue?logo=python&logoColor=white" alt="Python 3.8+">
    <img src="https://img.shields.io/badge/license-MIT-green" alt="MIT License">
    <img src="https://img.shields.io/badge/data-TCGA--STAD-orange" alt="TCGA-STAD">
    <img src="https://img.shields.io/badge/signatures-COSMIC%20v3.4-purple" alt="COSMIC v3.4">
    <img src="https://img.shields.io/badge/multi--omics-WES%20%2B%20Methylation%20%2B%20Clinical-red" alt="Multi-Omics">
  </p>
</p>

---

## 📋 Overview

This repository implements a comprehensive **multi-omics pipeline** for AI-driven molecular subtyping of **gastric cancer (TCGA-STAD)** by integrating:

| Data Layer | Source | Samples |
|---|---|---|
| **Somatic Mutations** | GDC Masked MAF files (WES) | 434 |
| **DNA Methylation** | Illumina 450K SeSAMe Level 3 β-values | 457 |
| **Clinical / Molecular** | GDC API + cBioPortal | 375 labeled |

The pipeline extracts de novo mutational signatures via NMF, maps them to **COSMIC v3.4** references, integrates methylation epigenomic features and gene-level mutation profiles, and classifies tumors into **5 molecular subtypes** (CIN, EBV, GS, MSI, POLE) using a novel **Hybrid Biological AI** system that combines rule-based biological logic with machine learning.

> **🔬 This is an improved version of [Machine-Learning-Mutational-Signatures](https://github.com/Ahsansayz/Machine-Learning-Mutational-Signatures).**  
> The original pipeline used only MAF-based mutation data. This version extends it with **DNA methylation data**, **gene-level mutation features**, and **clinical metadata** for a true multi-omics approach — achieving significantly improved classification performance.

---

## 🏆 Key Results

| Metric | Value |
|---|---|
| **Overall Accuracy** | **87.2%** |
| **Weighted F1-Score** | **0.856** |
| **Matthews Correlation Coefficient** | **0.820** |
| **MSI Recall** | **100%** (73/73) |
| **POLE Recall** | **100%** (8/8) |
| **EBV Recall** | **93%** |

### Classification Performance by Subtype

| Subtype | Precision | Recall | F1-Score | Support |
|---|---|---|---|---|
| CIN | 0.83 | 0.89 | 0.86 | 171 |
| EBV | 0.96 | 0.93 | 0.94 | 28 |
| GS | 0.73 | 0.67 | 0.70 | 45 |
| MSI | 0.93 | 1.00 | 0.97 | 73 |
| POLE | 1.00 | 1.00 | 1.00 | 8 |

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
 ┌────────────────────────────────────────────────────────────┐
 │                    INPUT DATA SOURCES                       │
 │  MAF Files (434)  │  Methylation (457)  │  Clinical (GDC)  │
 └───────┬───────────┴──────────┬──────────┴──────────┬───────┘
         │                      │                     │
         ▼                      │                     │
 ┌─ Step 1 ─────────────┐      │                     │
 │ Process MAF Files     │      │                     │
 │ • SBS96 Matrix        │      │                     │
 │ • Gene-Level Features │      │                     │
 └───────┬───────────────┘      │                     │
         ▼                      │                     │
 ┌─ Step 2 ─────────────┐      │                     │
 │ NMF Signature         │      │                     │
 │ Extraction (9 sigs)   │      │                     │
 └───────┬───────────────┘      │                     │
         ▼                      │                     │
 ┌─ Step 3 ─────────────┐      │                     │
 │ COSMIC v3.4 NNLS      │      │                     │
 │ Assignment (83 sigs)  │      │                     │
 └───────┬───────────────┘      │                     │
         │                      │                     │
         ▼                      ▼                     ▼
 ┌─ Step 4 ──────────────────────────────────────────────────┐
 │ Multi-Omics Feature Engineering                            │
 │ Signatures + Methylation + Clinical + Gene Mutations       │
 │ → 279 Features                                             │
 └───────┬───────────────────────────────────────────────────┘
         │
         ▼
 ┌─ Step 5 ─────────────┐
 │ DNA Methylation       │
 │ Processing (PCA,      │
 │ CIMP, Targeted CpGs)  │
 └───────┬───────────────┘
         │
         ▼
 ┌─ Step 6 ──────────────────────────────────────────────────┐
 │ Hybrid Biological AI Classification                        │
 │ ┌─────────────┐  ┌──────────────┐  ┌────────────────────┐│
 │ │ POLE Rules   │  │ GS Biology   │  │ XGBoost + SMOTE   ││
 │ │ SBS10+TMB    │  │ CDH1/RHOA/   │  │ Threshold Tuning  ││
 │ │ Detection    │  │ TP53 Logic   │  │ 500 trees, CV=5   ││
 │ └─────────────┘  └──────────────┘  └────────────────────┘│
 └───────┬───────────────────────────────────────────────────┘
         │
         ▼
 ┌─ Step 7 ─────────────┐
 │ Visualization &       │
 │ Interpretation        │
 │ SHAP, t-SNE, UMAP,   │
 │ Survival, TMB         │
 └───────────────────────┘
```

---

## 📁 Repository Structure

```
├── 01_raw_data/                  # Input: GDC MAF files + COSMIC reference
│   ├── maf_files/                   (434 masked somatic mutation MAFs)
│   └── cosmic_reference/            (COSMIC v3.4 SBS signatures)
│
├── 02_methylation_data/          # DNA methylation processed features
│   ├── methylation_features.csv     (457 samples × 44 features)
│   ├── methylation_pca.csv          (PCA-reduced methylation)
│   ├── methylation_targeted.csv     (Gene-promoter CpG features)
│   └── methylation_summary.csv      (Global methylation stats)
│
├── 03_pipeline_scripts/          # Core pipeline (7 steps + runner)
│   ├── step1_process_maf.py         → SBS96 matrix + gene features
│   ├── step2_extract_signatures.py  → NMF de novo extraction
│   ├── step3_cosmic_assignment.py   → COSMIC v3.4 NNLS fitting
│   ├── step4_build_features.py      → Multi-omics feature matrix
│   ├── step5_process_methylation.py → Methylation processing (RAM-optimized)
│   ├── step6_classify.py            → Hybrid Biological AI classifier
│   ├── step7_visualize.py           → Publication figures
│   └── run_pipeline.sh              → Master pipeline runner
│
├── 04_intermediate_data/         # Intermediate pipeline outputs
│   └── sbs96_matrix.csv             (434 × 96 mutation channels)
│
├── 07_figures/                   # Publication-ready figures
│   ├── shap_bar.png / shap_summary.png
│   ├── tsne_umap.png
│   ├── signature_subtype_heatmap.png
│   ├── survival_curves.png
│   └── tmb_distribution.png
│
├── 08_clinical_data/             # Clinical & molecular subtype data
│   ├── clinical_data.csv            (GDC + cBioPortal merged)
│   ├── cbio_clinical.csv
│   └── gdc_clinical.csv
│
├── 09_documentation/             # Detailed documentation
│   ├── Step_by_Step_Documentation.txt
│   ├── pipeline_README.txt
│   └── README_Hinglish.txt
│
├── requirements.txt
├── LICENSE                       # MIT License
└── README.md                     # This file
```

---

## ⚙️ Installation

### Prerequisites
- Python 3.8+
- ~8 GB RAM (for methylation processing)
- ~500 MB disk space (excluding raw data)

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

3. **Methylation Data** (optional): Download Illumina 450K SeSAMe Level 3 beta-value files from GDC for TCGA-STAD

---

## 🚀 Quick Start

### Run Full Pipeline

```bash
cd AI-Driven-Multi-Omics-Gastric-Cancer-Subtyping
bash 03_pipeline_scripts/run_pipeline.sh
```

### Run Individual Steps

```bash
# Step 1: Process MAF files → SBS96 matrix + gene features
python 03_pipeline_scripts/step1_process_maf.py

# Step 2: Extract de novo signatures via NMF
python 03_pipeline_scripts/step2_extract_signatures.py

# Step 3: Map to COSMIC v3.4 via NNLS
python 03_pipeline_scripts/step3_cosmic_assignment.py

# Step 4: Build multi-omics feature matrix
python 03_pipeline_scripts/step4_build_features.py

# Step 5: Process DNA methylation data
python 03_pipeline_scripts/step5_process_methylation.py

# Step 6: Hybrid Biological AI classification
python 03_pipeline_scripts/step6_classify.py

# Step 7: Generate publication figures
python 03_pipeline_scripts/step7_visualize.py
```

### Run Specific Steps Only

```bash
# Run only steps 6-7 (classification + visualization)
bash 03_pipeline_scripts/run_pipeline.sh 6 7
```

---

## 🧠 Methodology

### What's New vs. the Original Pipeline?

| Feature | [Original](https://github.com/Ahsansayz/Machine-Learning-Mutational-Signatures) | **This Version** |
|---|---|---|
| Input Data | MAF only | **MAF + Methylation + Clinical** |
| Features | 179 (signatures + clinical) | **279 (multi-omics)** |
| Gene Mutations | ❌ | ✅ CDH1, RHOA, TP53, PIK3CA + 16 more |
| Methylation | ❌ | ✅ Promoter CpGs, CIMP score, PCA |
| Classification | Standard XGBoost | **Hybrid Biological AI** |
| POLE Detection | ML only (poor recall) | **SBS10 + TMB biological rules (100%)** |
| GS Detection | ML only | **CDH1/RHOA/TP53 biology-assisted** |
| Best Accuracy | 81.1% | **87.2%** |
| MSI Recall | 100% | **100%** |

### Hybrid Biological AI System

The classifier combines three complementary strategies:

1. **Biological Rule Engine** — Uses domain knowledge to detect subtypes with strong genomic markers:
   - **POLE**: SBS10a/b/c/d signature burden + elevated TMB with MSI safeguard
   - **GS**: CDH1 truncating mutations + RHOA mutations + TP53 wild-type

2. **SMOTE-XGBoost ML Model** — Handles the remaining classification with class-imbalance correction via SMOTE oversampling

3. **Per-Class Threshold Optimization** — Optimizes decision boundaries independently for each subtype to maximize per-class F1 scores

### Multi-Omics Feature Engineering

The pipeline constructs **279 features** across four data modalities:

| Category | Features | Description |
|---|---|---|
| COSMIC Signatures | 83 raw + 83 proportional | NNLS-fitted signature activities |
| Engineered | 13 | TMB, MSI burden, signature diversity, etc. |
| Gene Mutations | 53 | Binary status + count for 20 key genes |
| Methylation | 44 | PCA components, CIMP score, promoter CpGs |
| Clinical | 3 | Age, stage, vital status |

### DNA Methylation Processing

The methylation module is **RAM-optimized** with a two-pass architecture:
- **Pass 1** (~50 MB RAM): Extracts targeted CpG probes (MLH1, CDH1, CDKN2A, MGMT, BRCA1, RUNX3) and computes summary statistics
- **Pass 2** (~500 MB RAM): PCA on top 5,000 most variable probes for unsupervised feature extraction

---

## 📊 Visualizations

The pipeline generates 6 publication-ready figures:

| Figure | Description |
|---|---|
| `shap_summary.png` | SHAP beeswarm plot showing feature impact per subtype |
| `shap_bar.png` | Mean absolute SHAP values (global importance) |
| `tsne_umap.png` | t-SNE & UMAP embeddings colored by molecular subtype |
| `signature_subtype_heatmap.png` | Mean COSMIC signature contribution per subtype |
| `survival_curves.png` | Kaplan-Meier overall survival curves with log-rank test |
| `tmb_distribution.png` | Tumor mutation burden box plots by subtype |

---

## 📚 References

- **TCGA-STAD**: The Cancer Genome Atlas - Stomach Adenocarcinoma ([Nature 2014](https://doi.org/10.1038/nature13480))
- **COSMIC v3.4**: Catalogue of Somatic Mutations in Cancer ([Tate et al., 2019](https://doi.org/10.1093/nar/gky1015))
- **NMF for Signatures**: Alexandrov et al., *Nature* 2013 ([doi:10.1038/nature12477](https://doi.org/10.1038/nature12477))
- **XGBoost**: Chen & Guestrin, *KDD* 2016 ([doi:10.1145/2939672.2939785](https://doi.org/10.1145/2939672.2939785))
- **SHAP**: Lundberg & Lee, *NeurIPS* 2017

---

## 📄 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- [TCGA Research Network](https://www.cancer.gov/tcga) for making genomic data publicly available
- [GDC Data Portal](https://portal.gdc.cancer.gov/) for data access infrastructure
- [cBioPortal](https://www.cbioportal.org/) for clinical data integration
- [COSMIC](https://cancer.sanger.ac.uk/signatures/) for the mutational signature reference database

---

<p align="center">
  <sub>Built with ❤️ for precision oncology research</sub>
</p>
