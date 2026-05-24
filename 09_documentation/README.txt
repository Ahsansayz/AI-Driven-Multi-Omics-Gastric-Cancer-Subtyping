================================================================================
              DOCUMENTATION — TABLE OF CONTENTS
================================================================================
  Project: AI-Driven Mutational Signature Analysis of Gastric Cancer
  Location: /home/param/Downloads/Project/documentation/
================================================================================


📄 MAIN DOCUMENTATION FILE:
─────────────────────────────
  Step_by_Step_Documentation.txt
      ↳ Complete step-by-step guide of EVERYTHING done in this project
      ↳ Every command, every output, every result — in sequence
      ↳ Covers: setup, dependencies, data download, all 8 pipeline steps,
        research paper generation, and file inventory


📂 scripts/ — FULL SOURCE CODE OF EACH SCRIPT:
──────────────────────────────────────────────────
  scripts/step1_build_sbs96_matrix.py.txt     → 96-ch SBS matrix builder
  scripts/step2_extract_signatures.py.txt     → NMF de novo extraction
  scripts/step3_cosmic_assignment.py.txt      → COSMIC v3.4 NNLS fitting
  scripts/step4_get_clinical_data.py.txt      → GDC + cBioPortal download
  scripts/step5_feature_engineering.py.txt    → Feature matrix construction
  scripts/step6_ml_classification.py.txt      → ML training & evaluation
  scripts/step7_visualize.py.txt              → SHAP, t-SNE, survival plots
  scripts/step8_generate_report.py.txt        → HTML report generator
  scripts/run_pipeline.sh.txt                 → Master runner script
  scripts/generate_paper.py.txt               → DOCX research paper generator


📂 figures/ — ALL PUBLICATION FIGURES:
────────────────────────────────────────
  figures/rank_selection.png                → NMF rank selection plot
  figures/signature_profiles_96ch.png       → 9 de novo signature profiles
  figures/cosmic_signature_landscape.png    → COSMIC signature landscape
  figures/cosine_similarity_distribution.png → Fitting quality histogram
  figures/best_confusion_matrix.png         → XGBoost confusion matrix
  figures/roc_curves.png                    → ROC curves (all subtypes)
  figures/model_comparison_chart.png        → 5-model performance comparison
  figures/feature_importance.png            → Top 20 feature importances
  figures/shap_bar.png                      → SHAP feature importance
  figures/tsne_umap.png                     → t-SNE + UMAP visualization
  figures/signature_subtype_heatmap.png     → Signature-subtype associations
  figures/survival_curves.png               → Kaplan-Meier OS curves
  figures/tmb_distribution.png              → TMB box plots by subtype


🔑 KEY RESULTS SUMMARY:
─────────────────────────
  ✅ 431 TCGA-STAD samples processed
  ✅ 154,534 somatic mutations in 96-channel SBS matrix
  ✅ 9 de novo signatures extracted
  ✅ 83 active COSMIC v3.4 signatures mapped (cos_sim = 0.933)
  ✅ 375 samples with official TCGA molecular subtypes
  ✅ Best model: XGBoost (81.1% accuracy, F1=0.765, AUC=0.870)
  ✅ MSI subtype: 100% correct (73/73)
  ✅ Full research paper: output/Research_Paper_(...).docx


================================================================================
