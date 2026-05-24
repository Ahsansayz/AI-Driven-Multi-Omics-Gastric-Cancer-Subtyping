================================================================================
  SIGPRO — AI-Driven Mutational Signature Analysis of Gastric Cancer
  Using TCGA-STAD Next-Generation Sequencing Data
================================================================================

PROJECT SUMMARY (English):
  This project analyzes 434 TCGA-STAD (Stomach Adenocarcinoma) whole-exome
  sequencing samples to extract mutational signatures, assign COSMIC v3.4
  references, and classify molecular subtypes (CIN/MSI/GS/EBV/POLE) using
  machine learning. Best result: XGBoost with 81.1% accuracy, MSI subtype
  classified at 100%.

PROJECT SUMMARY (Hinglish):
  Is project mein humne 434 stomach cancer patients ka mutation data liya
  TCGA se. Usse 96 tarah ke mutation patterns nikale. Phir NMF se 9
  signatures extract kiye, COSMIC database se match kiya, clinical subtypes
  download kiye, aur 5 ML models train karke 81% accuracy se cancer subtypes
  predict kiye. MSI subtype 100% sahi classify hua!


FOLDER STRUCTURE:
═════════════════

  sigpro/
  │
  ├── README.txt                    ← Ye file (master overview)
  │
  ├── 01_raw_data/                  ← INPUT: 434 MAF files from GDC
  │   ├── maf_files/                   (434 .maf files)
  │   ├── MANIFEST.txt                 (GDC download manifest)
  │   └── README.txt
  │
  ├── 02_pipeline_scripts/          ← ALL 10 Python/Bash scripts
  │   ├── step1_build_sbs96_matrix.py
  │   ├── step2_extract_signatures.py
  │   ├── step3_cosmic_assignment.py
  │   ├── step4_get_clinical_data.py
  │   ├── step5_feature_engineering.py
  │   ├── step6_ml_classification.py
  │   ├── step7_visualize.py
  │   ├── step8_generate_report.py
  │   ├── run_pipeline.sh              (master runner)
  │   ├── generate_paper.py            (DOCX paper generator)
  │   └── README.txt
  │
  ├── 03_sbs96_matrix/              ← Step 1 OUTPUT: Mutation matrix
  │   ├── sbs96_matrix.csv             (431 samples × 96 channels)
  │   └── README.txt
  │
  ├── 04_denovo_signatures/         ← Step 2 OUTPUT: NMF extraction
  │   ├── signature_profiles.csv       (9 signatures × 96 channels)
  │   ├── signature_activities.csv     (431 samples × 9 activities)
  │   ├── rank_selection_metrics.csv
  │   ├── rank_selection.png
  │   ├── signature_profiles_96ch.png
  │   ├── signature_activities_heatmap.png
  │   └── README.txt
  │
  ├── 05_cosmic_assignment/         ← Step 3 OUTPUT: COSMIC mapping
  │   ├── cosmic_activities.csv        (431 × 83 active sigs)
  │   ├── cosmic_cosine_similarities.csv
  │   ├── top_signatures_per_sample.csv
  │   ├── COSMIC_v3.4_SBS_GRCh38.txt  (reference file)
  │   ├── cosmic_signature_landscape.png
  │   ├── cosine_similarity_distribution.png
  │   └── README.txt
  │
  ├── 06_clinical_data/             ← Step 4 OUTPUT: Patient data
  │   ├── clinical_data.csv            (merged final file)
  │   ├── gdc_clinical.csv             (from GDC API)
  │   ├── cbio_clinical.csv            (from cBioPortal)
  │   └── README.txt
  │
  ├── 07_feature_matrix/            ← Step 5 OUTPUT: ML-ready data
  │   ├── ml_features.csv             (375 samples × 179 features)
  │   ├── ml_labels.csv               (375 subtype labels)
  │   └── README.txt
  │
  ├── 08_ml_results/                ← Step 6 OUTPUT: ML models
  │   ├── model_comparison.csv
  │   ├── predictions.csv
  │   ├── feature_importance.csv
  │   ├── best_model.pkl               (trained XGBoost model)
  │   ├── best_confusion_matrix.png
  │   ├── confusion_matrices_all.png
  │   ├── roc_curves.png
  │   ├── model_comparison_chart.png
  │   ├── feature_importance.png
  │   └── README.txt
  │
  ├── 09_figures/                   ← Step 7 OUTPUT: Publication plots
  │   ├── shap_summary.png
  │   ├── shap_bar.png
  │   ├── tsne_umap.png
  │   ├── signature_subtype_heatmap.png
  │   ├── survival_curves.png
  │   ├── tmb_distribution.png
  │   └── README.txt
  │
  ├── 10_reports/                   ← Step 8 OUTPUT: Final reports
  │   ├── report.html                  (interactive HTML dashboard)
  │   ├── Research_Paper_(...).docx    (full research paper)
  │   └── README.txt
  │
  ├── 11_documentation/             ← Step-by-step guides
  │   ├── Step_by_Step_Documentation.txt  (English, every command)
  │   ├── README_Hinglish.txt             (Hinglish explainer)
  │   └── README.txt
  │
  └── 12_old_exploratory/           ← Old/initial scripts (before rewrite)
      ├── build_matrix.py, extract_signatures.py, clustering.py, etc.
      └── README.txt


PIPELINE FLOW:
══════════════

  MAF files (434)
       │
       ▼
  ┌─ Step 1 ─┐   96-channel SBS trinucleotide matrix
  └──────────┘   (431 samples × 96 mutation types)
       │
       ▼
  ┌─ Step 2 ─┐   De novo NMF extraction → 9 signatures
  └──────────┘
       │
       ▼
  ┌─ Step 3 ─┐   COSMIC v3.4 NNLS fitting → 83 active signatures
  └──────────┘   (mean cosine similarity = 0.933)
       │
       ▼
  ┌─ Step 4 ─┐   Clinical data (GDC + cBioPortal) + TMB
  └──────────┘   → 375 samples with official TCGA subtypes
       │
       ▼
  ┌─ Step 5 ─┐   Feature engineering → 179 features
  └──────────┘   (83 raw + 83 proportions + 13 engineered)
       │
       ▼
  ┌─ Step 6 ─┐   ML Classification (RF, XGBoost, SVM, MLP, GB)
  └──────────┘   🏆 XGBoost: 81.1% acc, F1=0.765, AUC=0.870
       │
       ▼
  ┌─ Step 7 ─┐   Visualization (SHAP, t-SNE, UMAP, Survival)
  └──────────┘
       │
       ▼
  ┌─ Step 8 ─┐   HTML Report + DOCX Research Paper
  └──────────┘


KEY RESULTS:
════════════

  • 431 TCGA-STAD samples processed
  • 154,534 somatic mutations captured in 96-channel SBS matrix
  • 9 de novo signatures extracted (NMF, stability=1.0)
  • 83 active COSMIC v3.4 signatures (cosine similarity=0.933)
  • 375 samples with official molecular subtype labels
  • Best model: XGBoost (Accuracy=81.1%, F1=0.765, AUC=0.870)
  • MSI subtype: 100% correctly classified (73/73)
  • Top features: SBS54, SBS15, TMB, MSI_sig_burden_prop


HOW TO RE-RUN:
══════════════

  # 1. Activate environment
  source /home/param/miniconda3/etc/profile.d/conda.sh
  conda activate base

  # 2. Go to scripts folder
  cd /home/param/Downloads/sigpro/02_pipeline_scripts

  # 3. Run full pipeline
  bash run_pipeline.sh

  # Or individual steps
  python step1_build_sbs96_matrix.py
  python step6_ml_classification.py   # etc.


DEPENDENCIES:
═════════════

  pip install pandas numpy scikit-learn xgboost matplotlib seaborn \
    tqdm shap imbalanced-learn lifelines umap-learn SigProfilerAssignment

================================================================================
