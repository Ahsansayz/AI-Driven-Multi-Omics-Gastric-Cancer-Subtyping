#!/usr/bin/env python3
"""
DNA Methylation Data Processing for Multi-Omics Integration
=============================================================
RAM-OPTIMIZED VERSION: Uses chunked processing to stay under 8GB RAM.

Processes GDC SeSAMe Level 3 beta-value methylation files and extracts
features for integration with the mutational signature pipeline.

Input:  DNA methylation beta-value files from GDC (457 samples)
        Located in: ../DNA methylation/<uuid>/<file>.methylation_array.sesame.level3betas.txt

Output: output/methylation/methylation_features.csv
        output/methylation/methylation_targeted.csv
        output/methylation/methylation_pca.csv
        output/methylation/methylation_summary.csv
"""

import os
import sys
import json
import glob
import gc
import requests
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

try:
    from sklearn.decomposition import IncrementalPCA
    from sklearn.preprocessing import StandardScaler
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    print("⚠️ scikit-learn required. Install: pip install scikit-learn")

from tqdm import tqdm

# ============================================================
# Configuration
# ============================================================
METHYLATION_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "DNA methylation")
OUTPUT_DIR = "output/methylation"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# PCA settings
TOP_VARIABLE_PROBES = 5000
N_PCA_COMPONENTS = 30

# Quality filters
MAX_NA_FRACTION = 0.20

# ============================================================
# Biologically important CpG probes for gastric cancer subtypes
# ============================================================
TARGETED_PROBES = {
    'MLH1_promoter': [
        'cg00893636', 'cg03588108', 'cg11600697', 'cg13846866',
        'cg21490561', 'cg23658326', 'cg10990328', 'cg02279071',
    ],
    'CDH1_promoter': [
        'cg00419022', 'cg07420274', 'cg13924996', 'cg16218540',
        'cg19427610', 'cg24859243',
    ],
    'CDKN2A_promoter': [
        'cg04026675', 'cg06839912', 'cg08516516', 'cg12840719',
        'cg14069088', 'cg26349275',
    ],
    'MGMT_promoter': [
        'cg12434587', 'cg12981137',
    ],
    'BRCA1_promoter': [
        'cg08993267', 'cg19088651', 'cg04658354',
    ],
    'TP73_promoter': [
        'cg14094988', 'cg21163590',
    ],
    'RUNX3_promoter': [
        'cg02494853', 'cg15688816', 'cg18301702',
    ],
}

CIMP_PANEL_PROBES = [
    'cg00339556', 'cg01009664', 'cg01192891', 'cg01707559',
    'cg02941816', 'cg03217795', 'cg04523975', 'cg05254747',
    'cg06493994', 'cg07553761', 'cg08993267', 'cg10636246',
    'cg11284293', 'cg12434587', 'cg13846866', 'cg14094988',
    'cg15688816', 'cg16218540', 'cg17576918', 'cg18301702',
]

# Collect ALL probes we need for targeted features (so we can extract them in pass 1)
ALL_TARGETED_CG = set()
for probes in TARGETED_PROBES.values():
    ALL_TARGETED_CG.update(probes)
ALL_TARGETED_CG.update(CIMP_PANEL_PROBES)


def map_uuids_to_barcodes(methylation_dir):
    """Map GDC folder UUIDs to TCGA sample barcodes via GDC API.
    
    The folder name in a GDC download IS the GDC file UUID.
    We query the GDC /files endpoint with these folder UUIDs to get
    the associated TCGA sample barcodes.
    """

    cache_file = os.path.join(OUTPUT_DIR, "uuid_barcode_map.json")

    if os.path.exists(cache_file):
        print("  📋 Loading cached UUID→barcode mapping...")
        with open(cache_file) as f:
            mapping = json.load(f)
            if len(mapping) > 0:
                print(f"  Loaded {len(mapping)} cached mappings")
                return mapping
            else:
                print("  Cache was empty, re-querying...")

    print("  🌐 Querying GDC API for UUID→barcode mapping...")
    print("     (This will be cached for future runs)")

    folder_uuids = [d for d in os.listdir(methylation_dir)
                    if os.path.isdir(os.path.join(methylation_dir, d))
                    and not d.startswith('.')]

    uuid_to_barcode = {}
    batch_size = 50

    for i in tqdm(range(0, len(folder_uuids), batch_size), desc="  Querying GDC"):
        batch = folder_uuids[i:i + batch_size]

        payload = {
            "filters": {
                "op": "in",
                "content": {
                    "field": "file_id",
                    "value": batch
                }
            },
            "fields": "file_id,cases.samples.submitter_id,cases.submitter_id",
            "format": "JSON",
            "size": str(len(batch))
        }

        try:
            response = requests.post(
                "https://api.gdc.cancer.gov/files",
                json=payload,
                headers={"Content-Type": "application/json"},
                timeout=30
            )
            data = response.json()

            for hit in data.get('data', {}).get('hits', []):
                file_id = hit.get('file_id', '')
                cases = hit.get('cases', [])
                if cases:
                    samples = cases[0].get('samples', [])
                    if samples:
                        barcode = samples[0].get('submitter_id', '')
                        if barcode:
                            uuid_to_barcode[file_id] = barcode
                    elif 'submitter_id' in cases[0]:
                        barcode = cases[0]['submitter_id']
                        if barcode:
                            uuid_to_barcode[file_id] = barcode
        except Exception as e:
            print(f"  ⚠️ API batch failed: {e}")
            continue

    print(f"  Mapped {len(uuid_to_barcode)} / {len(folder_uuids)} folders to TCGA barcodes")

    with open(cache_file, 'w') as f:
        json.dump(uuid_to_barcode, f, indent=2)

    return uuid_to_barcode


# ============================================================
# PASS 1: Lightweight scan — only targeted probes + stats
# This uses ~50 MB RAM instead of 8 GB
# ============================================================
def pass1_targeted_and_stats(methylation_dir, uuid_to_barcode):
    """Extract targeted probe values and per-sample summary stats in one pass.
    
    Instead of loading the full 457×486K matrix, we read each sample file
    one at a time, extract only the ~50 targeted probes + compute summary
    stats, then discard the rest. Peak RAM: ~50 MB.
    """
    
    print("\n  ═══════════════════════════════════════════════════")
    print("  PASS 1: Targeted features + summary stats (low RAM)")
    print("  ═══════════════════════════════════════════════════")
    
    folders = [d for d in os.listdir(methylation_dir)
               if os.path.isdir(os.path.join(methylation_dir, d))
               and not d.startswith('.')]
    
    targeted_rows = {}   # barcode → {probe: beta}
    summary_rows = {}    # barcode → {stat: value}
    probe_sums = {}      # probe_id → running sum (for variance calc later)
    probe_sq_sums = {}   # probe_id → running sum of squares
    probe_counts = {}    # probe_id → count of non-NA values
    sample_count = 0
    failed = 0
    
    for folder_uuid in tqdm(folders, desc="  Pass 1 scanning"):
        barcode = uuid_to_barcode.get(folder_uuid)
        if not barcode:
            failed += 1
            continue
        
        folder_path = os.path.join(methylation_dir, folder_uuid)
        txt_files = glob.glob(os.path.join(folder_path, "*.level3betas.txt"))
        if not txt_files:
            failed += 1
            continue
        
        try:
            # Read the full file for this ONE sample
            df = pd.read_csv(txt_files[0], sep='\t', header=None,
                            names=['probe', 'beta'], na_values=['NA', 'NaN', ''],
                            dtype={'probe': str, 'beta': np.float32})
            df = df.set_index('probe')['beta']
            
            # --- Extract targeted probes ---
            targeted_vals = {}
            for probe_id in ALL_TARGETED_CG:
                if probe_id in df.index:
                    val = df[probe_id]
                    if not np.isnan(val):
                        targeted_vals[probe_id] = float(val)
            targeted_rows[barcode] = targeted_vals
            
            # --- Compute per-sample summary stats ---
            valid = df.dropna()
            stats = {}
            stats['global_mean_beta'] = float(valid.mean())
            stats['global_var_beta'] = float(valid.var())
            stats['frac_hypermethylated'] = float((valid > 0.7).mean())
            stats['frac_hypomethylated'] = float((valid < 0.3).mean())
            
            # Shannon entropy
            bins = np.histogram(valid.values, bins=10, range=(0, 1))[0]
            probs = bins / bins.sum()
            probs = probs[probs > 0]
            stats['methylation_entropy'] = float(-np.sum(probs * np.log2(probs)))
            
            # CIMP probes
            cimp_vals = [df[p] for p in CIMP_PANEL_PROBES if p in df.index and not np.isnan(df[p])]
            if cimp_vals:
                stats['cimp_score'] = float(np.mean(cimp_vals))
                stats['cimp_high'] = int(stats['cimp_score'] > 0.3)
            
            summary_rows[barcode] = stats
            
            # --- Accumulate running stats for variance calculation (for PCA probe selection) ---
            for probe_id, val in df.items():
                if not np.isnan(val):
                    probe_sums[probe_id] = probe_sums.get(probe_id, 0.0) + val
                    probe_sq_sums[probe_id] = probe_sq_sums.get(probe_id, 0.0) + val * val
                    probe_counts[probe_id] = probe_counts.get(probe_id, 0) + 1
            
            sample_count += 1
            
            # Free memory
            del df, valid
            
        except Exception as e:
            failed += 1
            continue
    
    if failed > 0:
        print(f"  ⚠️ Failed/unmapped: {failed} samples")
    print(f"  ✅ Loaded {sample_count} samples successfully")
    
    # --- Build targeted features DataFrame ---
    targeted_df = pd.DataFrame.from_dict(targeted_rows, orient='index')
    
    # Compute gene-level features
    targeted_features = pd.DataFrame(index=targeted_df.index)
    for gene_name, probes in TARGETED_PROBES.items():
        available = [p for p in probes if p in targeted_df.columns]
        if available:
            targeted_features[f'{gene_name}_mean'] = targeted_df[available].mean(axis=1)
            targeted_features[f'{gene_name}_max'] = targeted_df[available].max(axis=1)
            print(f"    {gene_name}: {len(available)}/{len(probes)} probes found")
        else:
            print(f"    ⚠️ {gene_name}: 0/{len(probes)} probes found")
    
    # --- Build summary features DataFrame ---
    summary_df = pd.DataFrame.from_dict(summary_rows, orient='index')
    
    # --- Compute probe variances from running stats ---
    print("\n  📊 Computing probe variances for PCA probe selection...")
    probe_variances = {}
    for probe_id in probe_sums:
        n = probe_counts[probe_id]
        if n > 1 and n >= sample_count * (1 - MAX_NA_FRACTION):
            mean = probe_sums[probe_id] / n
            var = (probe_sq_sums[probe_id] / n) - (mean * mean)
            if var > 0.001:  # variance filter
                probe_variances[probe_id] = var
    
    # Filter out sex/SNP probes
    probe_variances = {p: v for p, v in probe_variances.items() 
                       if not p.startswith(('ch.', 'rs'))}
    
    # Select top variable probes for PCA
    sorted_probes = sorted(probe_variances.items(), key=lambda x: x[1], reverse=True)
    top_pca_probes = [p for p, v in sorted_probes[:TOP_VARIABLE_PROBES]]
    print(f"  Selected {len(top_pca_probes)} high-variance probes for PCA")
    
    # Cleanup running stats
    del probe_sums, probe_sq_sums, probe_counts, probe_variances
    gc.collect()
    
    return targeted_features, summary_df, top_pca_probes, sample_count


# ============================================================
# PASS 2: Load ONLY top-variable probes for PCA (chunked)
# This uses ~500 MB RAM instead of 8 GB
# ============================================================
def pass2_pca_features(methylation_dir, uuid_to_barcode, top_pca_probes):
    """Load only the top variable probes and compute PCA incrementally.
    
    Instead of loading 486K probes, we load only 5000 probes per sample.
    Peak RAM: ~500 MB for the 457×5000 matrix.
    """
    
    print("\n  ═══════════════════════════════════════════════════")
    print(f"  PASS 2: PCA on {len(top_pca_probes)} selected probes")
    print("  ═══════════════════════════════════════════════════")
    
    if not HAS_SKLEARN:
        print("  ⚠️ scikit-learn not available, skipping PCA")
        return pd.DataFrame()
    
    top_set = set(top_pca_probes)
    
    folders = [d for d in os.listdir(methylation_dir)
               if os.path.isdir(os.path.join(methylation_dir, d))
               and not d.startswith('.')]
    
    pca_rows = {}
    
    for folder_uuid in tqdm(folders, desc="  Pass 2 loading"):
        barcode = uuid_to_barcode.get(folder_uuid)
        if not barcode:
            continue
        
        folder_path = os.path.join(methylation_dir, folder_uuid)
        txt_files = glob.glob(os.path.join(folder_path, "*.level3betas.txt"))
        if not txt_files:
            continue
        
        try:
            df = pd.read_csv(txt_files[0], sep='\t', header=None,
                            names=['probe', 'beta'], na_values=['NA', 'NaN', ''],
                            dtype={'probe': str, 'beta': np.float32})
            df = df.set_index('probe')['beta']
            
            # Extract ONLY the top variable probes
            row = {p: df.get(p, np.nan) for p in top_pca_probes}
            pca_rows[barcode] = row
            
            del df
            
        except Exception:
            continue
    
    # Build small matrix (457 × 5000) — ~9 MB
    print(f"\n  🔗 Building PCA matrix: {len(pca_rows)} × {len(top_pca_probes)}")
    pca_matrix = pd.DataFrame.from_dict(pca_rows, orient='index')
    pca_matrix = pca_matrix[top_pca_probes]  # ensure column order
    
    del pca_rows
    gc.collect()
    
    # Impute NaN with column median
    na_count = pca_matrix.isna().sum().sum()
    if na_count > 0:
        print(f"  🔧 Imputing {na_count} missing values (median)...")
        pca_matrix = pca_matrix.fillna(pca_matrix.median())
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(pca_matrix.values)
    X_scaled = np.nan_to_num(X_scaled, nan=0.0)
    
    # PCA
    n_comp = min(N_PCA_COMPONENTS, X_scaled.shape[0], X_scaled.shape[1])
    print(f"  🧮 Fitting PCA with {n_comp} components...")
    
    from sklearn.decomposition import PCA
    pca = PCA(n_components=n_comp, random_state=42)
    X_pca = pca.fit_transform(X_scaled)
    
    pca_cols = [f'meth_PC{i+1}' for i in range(n_comp)]
    pca_df = pd.DataFrame(X_pca, index=pca_matrix.index, columns=pca_cols)
    
    ev = pca.explained_variance_ratio_
    cum = np.cumsum(ev)
    print(f"    Variance explained: PC1={ev[0]:.3f}, PC1-5={cum[4]:.3f}, Total={cum[-1]:.3f}")
    
    del pca_matrix, X_scaled
    gc.collect()
    
    return pca_df


def align_with_mutation_data(methylation_features):
    """Align methylation features with existing mutation-based feature matrix."""

    print("\n  🔗 Aligning with mutation-based features...")

    ml_features_file = "output/ml_features.csv"
    if not os.path.exists(ml_features_file):
        print("  ⚠️ ml_features.csv not found — saving all methylation samples")
        return methylation_features, None

    mutation_features = pd.read_csv(ml_features_file, index_col=0)
    print(f"    Mutation features: {mutation_features.shape[0]} samples")
    print(f"    Methylation features: {methylation_features.shape[0]} samples")

    # Try direct match first
    common = mutation_features.index.intersection(methylation_features.index)

    if len(common) == 0:
        print("    Direct match failed, trying patient ID matching...")

        # Build maps: patient_id → full barcode
        mut_pid_map = {}
        for idx in mutation_features.index:
            pid = str(idx)[:15]
            if pid not in mut_pid_map:
                mut_pid_map[pid] = idx

        meth_pid_map = {}
        for idx in methylation_features.index:
            pid = str(idx)[:15]
            if pid not in meth_pid_map:
                meth_pid_map[pid] = idx

        # Try 12-char patient ID if 15-char doesn't work
        if not set(mut_pid_map.keys()).intersection(set(meth_pid_map.keys())):
            mut_pid_map = {}
            for idx in mutation_features.index:
                pid = str(idx)[:12]
                if pid not in mut_pid_map:
                    mut_pid_map[pid] = idx

            meth_pid_map = {}
            for idx in methylation_features.index:
                pid = str(idx)[:12]
                if pid not in meth_pid_map:
                    meth_pid_map[pid] = idx

        common_pids = set(mut_pid_map.keys()).intersection(set(meth_pid_map.keys()))
        print(f"    Patient ID matches: {len(common_pids)}")

        if common_pids:
            rename_map = {meth_pid_map[pid]: mut_pid_map[pid] for pid in common_pids}
            methylation_aligned = methylation_features.loc[
                [meth_pid_map[pid] for pid in common_pids]
            ].rename(index=rename_map)

            return methylation_aligned, mutation_features
    else:
        print(f"    Direct matches: {len(common)}")
        return methylation_features.loc[common], mutation_features

    return methylation_features, mutation_features


def main():
    print("=" * 60)
    print("DNA Methylation Processing (RAM-Optimized)")
    print("=" * 60)

    # Check input directory
    if not os.path.exists(METHYLATION_DIR):
        print(f"❌ Methylation directory not found: {METHYLATION_DIR}")
        print("   Please check the path.")
        sys.exit(1)

    folders = [d for d in os.listdir(METHYLATION_DIR)
               if os.path.isdir(os.path.join(METHYLATION_DIR, d))
               and not d.startswith('.')]
    print(f"\n  Found {len(folders)} sample folders")

    # Step 1: Map UUIDs to TCGA barcodes
    uuid_to_barcode = map_uuids_to_barcodes(METHYLATION_DIR)

    # Step 2 (PASS 1): Extract targeted + summary features (low RAM)
    targeted_features, summary_features, top_pca_probes, n_samples = \
        pass1_targeted_and_stats(METHYLATION_DIR, uuid_to_barcode)
    
    gc.collect()
    print(f"\n  📊 Pass 1 complete: {n_samples} samples processed")
    print(f"     Targeted features: {targeted_features.shape[1]}")
    print(f"     Summary features: {summary_features.shape[1]}")

    # Step 3 (PASS 2): PCA on top variable probes only (medium RAM)
    pca_features = pass2_pca_features(METHYLATION_DIR, uuid_to_barcode, top_pca_probes)
    
    gc.collect()

    # Step 4: Combine all methylation features
    print("\n  🔗 Combining all methylation features...")
    
    # Align indices across all three DataFrames
    common_idx = targeted_features.index.intersection(summary_features.index)
    if not pca_features.empty:
        common_idx = common_idx.intersection(pca_features.index)
    
    all_meth_features = pd.concat([
        targeted_features.loc[common_idx],
        summary_features.loc[common_idx],
        pca_features.loc[common_idx] if not pca_features.empty else pd.DataFrame(index=common_idx),
    ], axis=1)
    
    # Fill any remaining NaN
    all_meth_features = all_meth_features.fillna(all_meth_features.median())
    
    print(f"    Total methylation features: {all_meth_features.shape[1]}")
    print(f"    Samples: {all_meth_features.shape[0]}")

    # Step 5: Align with mutation data
    meth_aligned, mutation_features = align_with_mutation_data(all_meth_features)

    # Step 6: Save outputs
    print("\n  💾 Saving outputs...")

    all_meth_features.to_csv(os.path.join(OUTPUT_DIR, "methylation_features.csv"))
    print(f"    Saved: methylation_features.csv ({all_meth_features.shape})")

    targeted_features.to_csv(os.path.join(OUTPUT_DIR, "methylation_targeted.csv"))
    summary_features.to_csv(os.path.join(OUTPUT_DIR, "methylation_summary.csv"))
    if not pca_features.empty:
        pca_features.to_csv(os.path.join(OUTPUT_DIR, "methylation_pca.csv"))

    meth_aligned.to_csv(os.path.join(OUTPUT_DIR, "methylation_features_aligned.csv"))
    print(f"    Saved: methylation_features_aligned.csv ({meth_aligned.shape})")

    # Create combined multi-omics matrix
    if mutation_features is not None:
        common = mutation_features.index.intersection(meth_aligned.index)
        if len(common) > 0:
            multiomics = pd.concat([
                mutation_features.loc[common],
                meth_aligned.loc[common]
            ], axis=1)
            multiomics.to_csv("output/ml_features_multiomics.csv")
            print(f"\n  🧬 Multi-omics feature matrix saved!")
            print(f"     Shape: {multiomics.shape[0]} samples × {multiomics.shape[1]} features")
            print(f"     Mutation features: {mutation_features.shape[1]}")
            print(f"     Methylation features: {meth_aligned.shape[1]}")
        else:
            print("\n  ⚠️ No overlapping samples between mutation and methylation data")

    # Summary
    print(f"\n{'=' * 60}")
    print(f"📊 Methylation Processing Summary")
    print(f"{'=' * 60}")
    print(f"  Samples loaded:       {n_samples}")
    print(f"  Targeted features:    {targeted_features.shape[1]}")
    print(f"  Summary features:     {summary_features.shape[1]}")
    print(f"  PCA features:         {pca_features.shape[1] if not pca_features.empty else 0}")
    print(f"  Total meth features:  {all_meth_features.shape[1]}")
    print(f"  Aligned with mutations: {meth_aligned.shape[0]} samples")

    print(f"\n✅ Methylation processing complete!")
    print(f"📁 All outputs saved in: {OUTPUT_DIR}/")
    print(f"\n🚀 Next step: Run ml_classification.py with --multiomics flag")
    print(f"   or use output/ml_features_multiomics.csv as the feature matrix")


if __name__ == "__main__":
    main()
