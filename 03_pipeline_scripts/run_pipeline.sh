#!/bin/bash
# ============================================================
# TCGA-STAD Gastric Cancer Mutational Signature Pipeline
# ============================================================
# Full pipeline: MAF → Signatures → Features → Classification
#
# Usage: cd Dessertation_final && bash 03_pipeline_scripts/run_pipeline.sh
# ============================================================

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "============================================================"
echo "TCGA-STAD Gastric Cancer Mutational Signature Pipeline"
echo "============================================================"
echo "Project: $PROJECT_DIR"
echo ""

# ============================================================
# Setup: Create symlinks so scripts can find data
# ============================================================
echo "--- Setting up project structure ---"
cd "$PROJECT_DIR"

# Create symlinks to expected paths
[ ! -e maf_files ] && ln -sf 01_raw_data/maf_files maf_files
[ ! -e cosmic ]    && ln -sf 01_raw_data/cosmic_reference/cosmic cosmic 2>/dev/null || true
[ ! -d output ]    && mkdir -p output
[ ! -d data ]      && mkdir -p data

# Copy pre-existing output data to output/ if available
for f in 04_intermediate_data/sbs96_matrix.csv; do
    [ -f "$f" ] && [ ! -f output/$(basename "$f") ] && cp "$f" output/
done

# Copy clinical data if available
if [ -d 08_clinical_data ]; then
    cp 08_clinical_data/*.csv data/ 2>/dev/null || true
fi

# Copy signature/cosmic results if available
if [ -d 04_intermediate_data/signatures ]; then
    mkdir -p output/signatures
    cp 04_intermediate_data/signatures/* output/signatures/ 2>/dev/null || true
fi
if [ -d 04_intermediate_data/cosmic_assignment ]; then
    mkdir -p output/cosmic_assignment
    cp 04_intermediate_data/cosmic_assignment/* output/cosmic_assignment/ 2>/dev/null || true
fi

# Copy methylation data if available
if [ -d 02_methylation_data ]; then
    mkdir -p output/methylation
    cp 02_methylation_data/*.csv output/methylation/ 2>/dev/null || true
    cp 02_methylation_data/*.json output/methylation/ 2>/dev/null || true
fi

# Copy pre-computed feature matrices if available
for f in 05_feature_matrices/ml_features.csv 05_feature_matrices/ml_labels.csv \
         05_feature_matrices/ml_features_multiomics.csv \
         05_feature_matrices/ml_features_multiomics_v2.csv \
         05_feature_matrices/gene_mutation_features.csv; do
    [ -f "$f" ] && [ ! -f output/$(basename "$f") ] && cp "$f" output/
done

echo "  ✅ Project structure ready"
echo ""

# ============================================================
# Parse arguments
# ============================================================
START_STEP=${1:-1}
END_STEP=${2:-7}

echo "Running steps $START_STEP to $END_STEP"
echo ""

# ============================================================
# Step 1: Process MAF files (SBS96 + Gene Features)
# ============================================================
if [ "$START_STEP" -le 1 ] && [ "$END_STEP" -ge 1 ]; then
    echo "============================================"
    echo "STEP 1: Process MAF Files"
    echo "============================================"
    python3 "$SCRIPT_DIR/step1_process_maf.py"
    echo ""
fi

# ============================================================
# Step 2: Extract De Novo Signatures (NMF)
# ============================================================
if [ "$START_STEP" -le 2 ] && [ "$END_STEP" -ge 2 ]; then
    echo "============================================"
    echo "STEP 2: Extract Signatures (NMF)"
    echo "============================================"
    python3 "$SCRIPT_DIR/step2_extract_signatures.py"
    echo ""
fi

# ============================================================
# Step 3: COSMIC Assignment
# ============================================================
if [ "$START_STEP" -le 3 ] && [ "$END_STEP" -ge 3 ]; then
    echo "============================================"
    echo "STEP 3: COSMIC Assignment"
    echo "============================================"
    python3 "$SCRIPT_DIR/step3_cosmic_assignment.py"
    echo ""
fi

# ============================================================
# Step 4: Build Feature Matrix (Clinical + Signatures + Multi-omics)
# ============================================================
if [ "$START_STEP" -le 4 ] && [ "$END_STEP" -ge 4 ]; then
    echo "============================================"
    echo "STEP 4: Build Feature Matrix"
    echo "============================================"
    python3 "$SCRIPT_DIR/step4_build_features.py"
    echo ""
fi

# ============================================================
# Step 5: Process Methylation (if raw data available)
# ============================================================
if [ "$START_STEP" -le 5 ] && [ "$END_STEP" -ge 5 ]; then
    echo "============================================"
    echo "STEP 5: Process Methylation"
    echo "============================================"
    if [ -f output/methylation/methylation_features.csv ]; then
        echo "  Methylation features already exist — skipping"
        echo "  (Delete output/methylation/ to re-process)"
    else
        python3 "$SCRIPT_DIR/step5_process_methylation.py"
    fi
    echo ""
fi

# ============================================================
# Step 6: ML Classification (Hybrid Biological AI)
# ============================================================
if [ "$START_STEP" -le 6 ] && [ "$END_STEP" -ge 6 ]; then
    echo "============================================"
    echo "STEP 6: ML Classification"
    echo "============================================"
    python3 "$SCRIPT_DIR/step6_classify.py"
    echo ""
fi

# ============================================================
# Step 7: Visualization
# ============================================================
if [ "$START_STEP" -le 7 ] && [ "$END_STEP" -ge 7 ]; then
    echo "============================================"
    echo "STEP 7: Visualization"
    echo "============================================"
    python3 "$SCRIPT_DIR/step7_visualize.py"
    echo ""
fi

# ============================================================
# Copy results back to organized folders
# ============================================================
echo "============================================"
echo "Organizing results..."
echo "============================================"

# Copy feature matrices
for f in output/ml_features.csv output/ml_labels.csv output/ml_features_multiomics.csv \
         output/ml_features_multiomics_v2.csv output/gene_mutation_features.csv; do
    [ -f "$f" ] && cp "$f" 05_feature_matrices/ 2>/dev/null
done

# Copy ML results
[ -d output/ml_results_hybrid ] && cp -r output/ml_results_hybrid 06_ml_results/stage4_hybrid_final 2>/dev/null || true

# Copy figures
[ -d output/figures ] && cp output/figures/* 07_figures/ 2>/dev/null || true

echo "  ✅ Results organized"
echo ""
echo "============================================================"
echo "✅ Pipeline complete!"
echo "============================================================"
echo ""
echo "Key outputs:"
echo "  Feature matrix:  output/ml_features_multiomics_v2.csv"
echo "  ML results:      output/ml_results_hybrid/"
echo "  Figures:         07_figures/"
