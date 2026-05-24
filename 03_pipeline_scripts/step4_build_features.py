#!/usr/bin/env python3
"""
Step 4: Build Feature Matrix (Clinical + Signature + Methylation)
=================================================================
MERGED: get_clinical_data.py + feature_engineering.py

1. Downloads clinical data & molecular subtypes from GDC/cBioPortal
2. Calculates TMB from MAF files
3. Engineers features from COSMIC signature activities
4. Integrates methylation features (if available)
5. Merges gene-level mutation features (from step1)
6. Produces final multi-omics feature matrix

Input:  output/cosmic_assignment/cosmic_activities.csv (from Step 3)
        output/gene_mutation_features.csv (from Step 1)
        output/methylation/ (from Step 5, optional)
        maf_files/*.maf (for barcodes + TMB)
Output: data/clinical_data.csv
        output/ml_features.csv
        output/ml_features_multiomics.csv (if methylation available)
        output/ml_features_multiomics_v2.csv (if gene features available)
        output/ml_labels.csv
"""

import os, sys, json, ssl
import pandas as pd
import urllib.request
import urllib.parse
import numpy as np
import warnings
warnings.filterwarnings('ignore')

OUTPUT_DIR = "output"
CLINICAL_DIR = "data"
MAF_DIR = "maf_files"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(CLINICAL_DIR, exist_ok=True)


# ============================================================
# PART A: Clinical Data & TMB
# ============================================================

def extract_barcodes_from_mafs():
    """Extract TCGA barcodes from MAF files."""
    print("  Extracting TCGA barcodes from MAF files...")
    barcodes = {}
    maf_files = sorted([f for f in os.listdir(MAF_DIR) if f.endswith('.maf')])
    
    for filename in maf_files:
        filepath = os.path.join(MAF_DIR, filename)
        try:
            df = pd.read_csv(filepath, sep='\t', comment='#',
                             usecols=['Tumor_Sample_Barcode'], nrows=1)
            barcode = df['Tumor_Sample_Barcode'].iloc[0]
            barcodes[filename] = {
                'full_barcode': barcode,
                'patient_id': barcode[:12],
                'sample_id': barcode[:15]
            }
        except Exception as e:
            print(f"    Warning: Error reading {filename}: {e}")
    
    print(f"  Found {len(barcodes)} unique samples")
    return barcodes


def download_clinical_from_gdc(patient_ids):
    """Download clinical data from GDC API."""
    print("\n  Downloading clinical data from GDC API...")
    url = "https://api.gdc.cancer.gov/cases"
    all_clinical = []
    batch_size = 100
    patient_list = list(set(patient_ids))
    
    fields = ",".join([
        "submitter_id", "demographic.gender", "demographic.race",
        "demographic.vital_status", "demographic.days_to_death",
        "diagnoses.age_at_diagnosis", "diagnoses.tumor_stage",
        "diagnoses.primary_diagnosis", "diagnoses.days_to_last_follow_up",
        "project.project_id"
    ])
    
    for i in range(0, len(patient_list), batch_size):
        batch = patient_list[i:i+batch_size]
        filters = json.dumps({"op": "in", "content": {"field": "submitter_id", "value": batch}})
        
        # Use POST to avoid URL length/encoding issues
        try:
            post_data = urllib.parse.urlencode({
                "filters": filters,
                "fields": fields,
                "size": str(batch_size),
                "format": "json"
            }).encode('utf-8')
            
            req = urllib.request.Request(url, data=post_data, method='POST')
            req.add_header('Content-Type', 'application/x-www-form-urlencoded')
            
            with urllib.request.urlopen(req, timeout=30) as response:
                data = json.loads(response.read().decode())
            
            for hit in data.get('data', {}).get('hits', []):
                clinical = {'patient_id': hit.get('submitter_id', '')}
                demo = hit.get('demographic', {})
                clinical['gender'] = demo.get('gender', '')
                clinical['vital_status'] = demo.get('vital_status', '')
                clinical['days_to_death'] = demo.get('days_to_death', '')
                
                diags = hit.get('diagnoses', [{}])
                if diags:
                    diag = diags[0]
                    clinical['age_at_diagnosis'] = diag.get('age_at_diagnosis', '')
                    clinical['tumor_stage'] = diag.get('tumor_stage', '')
                    clinical['days_to_last_follow_up'] = diag.get('days_to_last_follow_up', '')
                
                all_clinical.append(clinical)
        except Exception as e:
            print(f"    Warning: GDC API batch {i//batch_size + 1} error: {e}")
    
    print(f"  Downloaded clinical data for {len(all_clinical)} patients")
    return all_clinical


def download_subtypes():
    """Download molecular subtypes from cBioPortal."""
    print("\n  Downloading molecular subtypes from cBioPortal...")
    
    # Create SSL context that handles certificate issues
    ctx = ssl.create_default_context()
    try:
        # Try normal SSL first
        url = "https://www.cbioportal.org/api/clinical-data/fetch?clinicalDataType=PATIENT&projection=SUMMARY"
        body = {"attributeIds": ["SUBTYPE"],
                "studyViewFilter": {"studyIds": ["stad_tcga_pan_can_atlas_2018"]}}
        
        req = urllib.request.Request(url, data=json.dumps(body).encode('utf-8'), method='POST')
        req.add_header('Content-Type', 'application/json')
        req.add_header('Accept', 'application/json')
        
        try:
            response = urllib.request.urlopen(req, timeout=30, context=ctx)
        except ssl.SSLError:
            # Fallback: skip SSL verification for cBioPortal
            ctx_unverified = ssl.create_default_context()
            ctx_unverified.check_hostname = False
            ctx_unverified.verify_mode = ssl.CERT_NONE
            response = urllib.request.urlopen(req, timeout=30, context=ctx_unverified)
            print("  (Used fallback SSL context)")
        
        data = json.loads(response.read().decode())
        response.close()
        
        subtypes = {item.get('patientId', ''): item.get('value', '') for item in data}
        print(f"  Got subtypes for {len(subtypes)} patients")
        return subtypes
    except Exception as e:
        print(f"  Warning: Could not download subtypes: {e}")
        return {}


def generate_tmb():
    """Calculate Tumor Mutation Burden for each sample."""
    print("\n  Calculating TMB...")
    tmb = {}
    exome_size_mb = 30.0
    
    for filename in os.listdir(MAF_DIR):
        if not filename.endswith('.maf'):
            continue
        filepath = os.path.join(MAF_DIR, filename)
        try:
            df = pd.read_csv(filepath, sep='\t', comment='#',
                             usecols=['Tumor_Sample_Barcode', 'Variant_Classification', 'Variant_Type'])
            barcode = df['Tumor_Sample_Barcode'].iloc[0]
            nonsyn_classes = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del',
                              'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins',
                              'Splice_Site', 'Nonstop_Mutation', 'Translation_Start_Site']
            nonsyn_count = df[df['Variant_Classification'].isin(nonsyn_classes)].shape[0]
            tmb[barcode] = {'total_mutations': df.shape[0], 'nonsynonymous_mutations': nonsyn_count,
                            'TMB': nonsyn_count / exome_size_mb}
        except:
            pass
    
    print(f"  Calculated TMB for {len(tmb)} samples")
    return tmb


# ============================================================
# PART B: Feature Engineering
# ============================================================

def load_signature_data():
    """Load signature activities — prefer COSMIC, fall back to de novo."""
    cosmic_file = "output/cosmic_assignment/cosmic_activities.csv"
    denovo_file = "output/signatures/signature_activities.csv"
    
    if os.path.exists(cosmic_file):
        print("  Loading COSMIC signature activities...")
        activities = pd.read_csv(cosmic_file, index_col=0)
        source = "COSMIC"
    elif os.path.exists(denovo_file):
        print("  Loading de novo signature activities...")
        activities = pd.read_csv(denovo_file, index_col=0)
        source = "de_novo"
    else:
        print("  No signature activities found! Run Steps 2-3 first.")
        sys.exit(1)
    
    # Remove zero-variance signatures
    nonzero_sigs = activities.columns[activities.var() > 0]
    activities = activities[nonzero_sigs]
    print(f"  Loaded {source} activities: {activities.shape[0]} samples x {activities.shape[1]} signatures")
    return activities, source


def engineer_features(activities, clinical_df=None, source="COSMIC"):
    """Build the feature matrix with engineered features."""
    features = activities.copy()
    total_activity = features.sum(axis=1)
    proportions = features.div(total_activity + 1e-10, axis=0)
    proportions.columns = [f"{c}_prop" for c in proportions.columns]
    
    eng = pd.DataFrame(index=features.index)
    
    if source == "COSMIC":
        # MSI signature burden
        msi_sigs = [s for s in ['SBS6', 'SBS14', 'SBS15', 'SBS20', 'SBS21', 'SBS26', 'SBS44'] if s in features.columns]
        if msi_sigs:
            eng['MSI_sig_burden'] = features[msi_sigs].sum(axis=1)
            eng['MSI_sig_burden_prop'] = proportions[[f"{s}_prop" for s in msi_sigs if f"{s}_prop" in proportions.columns]].sum(axis=1)
        
        # APOBEC burden
        apobec_sigs = [s for s in ['SBS2', 'SBS13'] if s in features.columns]
        if apobec_sigs:
            eng['APOBEC_burden'] = features[apobec_sigs].sum(axis=1)
            eng['APOBEC_burden_prop'] = proportions[[f"{s}_prop" for s in apobec_sigs if f"{s}_prop" in proportions.columns]].sum(axis=1)
        
        # Clock-like signatures
        clock_sigs = [s for s in ['SBS1', 'SBS5'] if s in features.columns]
        if clock_sigs:
            eng['clock_burden'] = features[clock_sigs].sum(axis=1)
        
        if 'SBS3' in features.columns:
            eng['HRD_burden'] = features['SBS3']
        
        sbs17_sigs = [s for s in ['SBS17a', 'SBS17b'] if s in features.columns]
        if sbs17_sigs:
            eng['SBS17_burden'] = features[sbs17_sigs].sum(axis=1)
        
        if 'SBS1' in features.columns and 'SBS5' in features.columns:
            eng['SBS1_SBS5_ratio'] = features['SBS1'] / (features['SBS5'] + 1e-10)
    
    eng['dominant_sig_prop'] = proportions.max(axis=1)
    eng['n_active_sigs'] = (proportions > 0.05).sum(axis=1)
    eng['log_total_activity'] = np.log1p(total_activity)
    
    feature_matrix = pd.concat([features, proportions, eng], axis=1)
    
    # Add clinical features (TMB, age, gender)
    if clinical_df is not None and len(clinical_df) > 0:
        if 'full_barcode' in clinical_df.columns:
            clinical_indexed = clinical_df.drop_duplicates('full_barcode').set_index('full_barcode')
        elif 'patient_id' in clinical_df.columns:
            clinical_indexed = clinical_df.drop_duplicates('patient_id').set_index('patient_id')
        else:
            clinical_indexed = pd.DataFrame()
        
        if not clinical_indexed.empty:
            for col_name, src_col in [('TMB', 'TMB'), ('age_years', 'age_years'), ('gender_male', 'gender')]:
                if src_col not in clinical_indexed.columns:
                    continue
                val_map = {}
                for sample in feature_matrix.index:
                    for key in [sample, sample[:12]]:
                        if key in clinical_indexed.index:
                            val = clinical_indexed.loc[key, src_col]
                            if col_name == 'gender_male':
                                val_map[sample] = 1 if str(val).lower() == 'male' else 0
                            else:
                                val_map[sample] = val
                            break
                if val_map:
                    feature_matrix[col_name] = feature_matrix.index.map(val_map)
                    if col_name == 'TMB':
                        feature_matrix['log_TMB'] = np.log1p(feature_matrix['TMB'].fillna(0))
    
    return feature_matrix.fillna(0)


# ============================================================
# PART C: Multi-omics Integration
# ============================================================

def integrate_methylation(feature_matrix):
    """Integrate methylation features if available."""
    meth_path = "output/methylation/methylation_features_aligned.csv"
    if not os.path.exists(meth_path):
        meth_path = "output/methylation/methylation_features.csv"
    
    if not os.path.exists(meth_path):
        print("  No methylation features found — skipping integration")
        return None
    
    print("  Integrating methylation features...")
    meth_df = pd.read_csv(meth_path, index_col=0)
    
    # Align by patient ID
    feat_pid = {str(idx)[:16]: idx for idx in feature_matrix.index}
    meth_pid = {str(idx)[:16]: idx for idx in meth_df.index}
    common_pids = set(feat_pid.keys()) & set(meth_pid.keys())
    
    if not common_pids:
        # Try 12-char patient ID
        feat_pid = {str(idx)[:12]: idx for idx in feature_matrix.index}
        meth_pid = {str(idx)[:12]: idx for idx in meth_df.index}
        common_pids = set(feat_pid.keys()) & set(meth_pid.keys())
    
    if common_pids:
        meth_rename = {meth_pid[pid]: feat_pid[pid] for pid in common_pids}
        meth_aligned = meth_df.loc[list(meth_rename.keys())].rename(index=meth_rename)
        
        feat_common = feature_matrix.loc[list(meth_rename.values())]
        multiomics = pd.concat([feat_common, meth_aligned.loc[feat_common.index]], axis=1)
        
        print(f"  Multi-omics: {multiomics.shape[0]} samples x {multiomics.shape[1]} features")
        return multiomics
    
    return None


def integrate_gene_features(feature_matrix):
    """Integrate gene-level mutation features from step1."""
    gene_path = os.path.join(OUTPUT_DIR, "gene_mutation_features.csv")
    if not os.path.exists(gene_path):
        print("  No gene mutation features found — skipping")
        return None
    
    print("  Integrating gene mutation features...")
    gene_df = pd.read_csv(gene_path, index_col=0)
    
    # Align
    feat_pid = {str(idx)[:16]: idx for idx in feature_matrix.index}
    gene_pid = {str(idx)[:16]: idx for idx in gene_df.index}
    common_pids = set(feat_pid.keys()) & set(gene_pid.keys())
    
    if common_pids:
        rename_map = {gene_pid[pid]: feat_pid[pid] for pid in common_pids}
        gene_aligned = gene_df.loc[list(rename_map.keys())].rename(index=rename_map)
        
        existing_cols = set(feature_matrix.columns)
        new_cols = [c for c in gene_aligned.columns if c not in existing_cols]
        
        enhanced = pd.concat([
            feature_matrix.loc[list(rename_map.values())],
            gene_aligned[new_cols].loc[list(rename_map.values())]
        ], axis=1)
        
        print(f"  Enhanced: {enhanced.shape[0]} samples x {enhanced.shape[1]} features (+{len(new_cols)} gene-level)")
        return enhanced
    
    return None


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 60)
    print("Step 4: Build Feature Matrix")
    print("=" * 60)
    
    # --- A: Clinical data ---
    print("\n--- A: Clinical Data ---")
    barcodes = extract_barcodes_from_mafs()
    barcode_df = pd.DataFrame.from_dict(barcodes, orient='index')
    patient_ids = barcode_df['patient_id'].unique().tolist()
    
    clinical_data = download_clinical_from_gdc(patient_ids)
    subtypes = download_subtypes()
    tmb_data = generate_tmb()
    
    # Build master clinical table
    master = barcode_df.reset_index().rename(columns={'index': 'maf_file'})
    clinical_df = pd.DataFrame(clinical_data) if clinical_data else pd.DataFrame()
    if not clinical_df.empty:
        master = master.merge(clinical_df, on='patient_id', how='left')
    if subtypes:
        master['molecular_subtype'] = master['patient_id'].map(subtypes)
    
    tmb_df = pd.DataFrame.from_dict(tmb_data, orient='index')
    if not tmb_df.empty:
        master = master.merge(tmb_df, left_on='full_barcode', right_index=True, how='left')
    
    if 'age_at_diagnosis' in master.columns:
        master['age_years'] = pd.to_numeric(master['age_at_diagnosis'], errors='coerce') / 365.25
    
    clinical_path = os.path.join(CLINICAL_DIR, "clinical_data.csv")
    master.to_csv(clinical_path, index=False)
    print(f"  Saved: {clinical_path}")
    
    # --- B: Feature engineering ---
    print("\n--- B: Feature Engineering ---")
    activities, source = load_signature_data()
    feature_matrix = engineer_features(activities, master, source)
    
    # Extract labels
    labels = None
    if 'molecular_subtype' in master.columns:
        label_rows = master.dropna(subset=['molecular_subtype'])[['full_barcode', 'molecular_subtype']]
        if len(label_rows) > 0:
            labels = label_rows.drop_duplicates('full_barcode').set_index('full_barcode')
    
    # Fallback: use pre-existing labels if API download failed
    existing_labels_path = os.path.join(OUTPUT_DIR, "ml_labels.csv")
    if (labels is None or len(labels) == 0) and os.path.exists(existing_labels_path):
        print("  API labels unavailable — using pre-existing ml_labels.csv")
        labels = pd.read_csv(existing_labels_path, index_col=0)
        print(f"  Loaded {len(labels)} pre-existing labels")
    
    # Align features & labels
    if labels is not None and len(labels) > 0:
        common = feature_matrix.index.intersection(labels.index)
        if len(common) == 0:
            # Patient ID matching fallback
            feat_pids = {idx: idx[:12] for idx in feature_matrix.index}
            label_pids = {idx: idx[:12] for idx in labels.index}
            pid_to_label = {label_pids[idx]: labels.loc[idx, 'molecular_subtype'] for idx in labels.index}
            matched = {idx: pid_to_label[pid] for idx, pid in feat_pids.items() if pid in pid_to_label}
            if matched:
                labels = pd.DataFrame.from_dict(matched, orient='index', columns=['molecular_subtype'])
                common = feature_matrix.index.intersection(labels.index)
        
        if len(common) > 0:
            feature_matrix.loc[common].to_csv(os.path.join(OUTPUT_DIR, "ml_features.csv"))
            labels.loc[common].to_csv(os.path.join(OUTPUT_DIR, "ml_labels.csv"))
            print(f"  Saved: ml_features.csv ({len(common)} labeled samples x {feature_matrix.shape[1]} features)")
            print(f"  Saved: ml_labels.csv")
        else:
            feature_matrix.to_csv(os.path.join(OUTPUT_DIR, "ml_features.csv"))
            print(f"  Warning: No labels could be matched to features")
    else:
        feature_matrix.to_csv(os.path.join(OUTPUT_DIR, "ml_features.csv"))
        print(f"  Saved: ml_features.csv (no labels available)")
    
    # --- C: Multi-omics ---
    print("\n--- C: Multi-omics Integration ---")
    multiomics = integrate_methylation(feature_matrix)
    if multiomics is not None:
        multiomics.to_csv(os.path.join(OUTPUT_DIR, "ml_features_multiomics.csv"))
        print(f"  Saved: ml_features_multiomics.csv")
        
        # Add gene features on top
        enhanced = integrate_gene_features(multiomics)
        if enhanced is not None:
            enhanced.to_csv(os.path.join(OUTPUT_DIR, "ml_features_multiomics_v2.csv"))
            print(f"  Saved: ml_features_multiomics_v2.csv")
    else:
        # Try gene features on base matrix
        enhanced = integrate_gene_features(feature_matrix)
        if enhanced is not None:
            enhanced.to_csv(os.path.join(OUTPUT_DIR, "ml_features_multiomics_v2.csv"))
    
    # Summary
    print(f"\n{'=' * 60}")
    print("Step 4 Summary")
    print("=" * 60)
    print(f"  Samples: {feature_matrix.shape[0]}")
    print(f"  Base features: {feature_matrix.shape[1]}")
    if multiomics is not None:
        print(f"  Multi-omics features: {multiomics.shape[1]}")
    if 'molecular_subtype' in master.columns:
        print(f"\n  Subtype Distribution:")
        for sub, cnt in master['molecular_subtype'].value_counts().items():
            print(f"    {sub}: {cnt}")
    
    print(f"\nStep 4 complete!")


if __name__ == "__main__":
    main()
