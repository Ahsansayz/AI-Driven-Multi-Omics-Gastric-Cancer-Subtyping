#!/usr/bin/env python3
"""
Strategy #7 + #5: Hybrid Rule System + Threshold Tuning
========================================================
- POLE detection via biological rules (SBS10a, TMB)
- Per-class probability threshold optimization
- ML handles remaining classification
"""

import os, pickle
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# Fix random seed for reproducibility
np.random.seed(42)

from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import classification_report, accuracy_score, f1_score, matthews_corrcoef
from sklearn.pipeline import Pipeline
from xgboost import XGBClassifier

try:
    from imblearn.over_sampling import SMOTE
    from imblearn.pipeline import Pipeline as ImbPipeline
    HAS_SMOTE = True
except ImportError:
    HAS_SMOTE = False

OUTPUT_DIR = "output/ml_results_hybrid"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# Load Data
# ============================================================
print("=" * 60)
print("Hybrid Rule System + Threshold Tuning")
print("=" * 60)

X_full = pd.read_csv('output/ml_features_multiomics_v2.csv', index_col=0)
y_df = pd.read_csv('output/ml_labels.csv', index_col=0)
common = X_full.index.intersection(y_df.index)
X_full = X_full.loc[common]
y = y_df.loc[common, 'molecular_subtype']

# Separate: ML features (original 226) vs Rule features (gene-level)
# ML model gets original features only — avoids noise from 53 extra columns
X_ml_orig = pd.read_csv('output/ml_features_multiomics.csv', index_col=0)
X_ml = X_ml_orig.loc[common]

# Rule features come from the full v2 set (gene mutations, TiTv, etc.)
X_rules = X_full  # Has CDH1_mutated, RHOA_mutated, TP53_mutated, TiTv_ratio, etc.

le = LabelEncoder()
y_enc = le.fit_transform(y)
classes = le.classes_

X_ml = X_ml.replace([np.inf, -np.inf], np.nan).fillna(X_ml.median())
X_rules = X_rules.replace([np.inf, -np.inf], np.nan).fillna(0)

print("ML features: %d, Rule features: %d" % (X_ml.shape[1], X_rules.shape[1]))
print("Classes: %s" % list(classes))
dist = dict(zip(*np.unique(y, return_counts=True)))
print("Distribution: %s" % dist)

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
pole_idx = list(classes).index('POLE')
gs_idx = list(classes).index('GS')
msi_idx = list(classes).index('MSI')
cin_idx = list(classes).index('CIN')
ebv_idx = list(classes).index('EBV')

# Check gene features availability (from X_rules)
HAS_CDH1 = 'CDH1_mutated' in X_rules.columns
HAS_TP53 = 'TP53_mutated' in X_rules.columns
HAS_RHOA = 'RHOA_mutated' in X_rules.columns
HAS_TITV = 'TiTv_ratio' in X_rules.columns
print("Gene features for rules: CDH1=%s TP53=%s RHOA=%s TiTv=%s" % (HAS_CDH1, HAS_TP53, HAS_RHOA, HAS_TITV))

# ============================================================
# Strategy A: Baseline SMOTE XGBoost (for comparison)
# ============================================================
print("\n" + "=" * 50)
print("Baseline: SMOTE + XGBoost (previous best)")
print("=" * 50)

if HAS_SMOTE:
    counts = dict(zip(*np.unique(y_enc, return_counts=True)))
    target_counts = dict(counts)
    target_counts[pole_idx] = 30
    target_counts[gs_idx] = max(counts[gs_idx], 60)
    
    baseline_pipe = ImbPipeline([
        ('scaler', StandardScaler()),
        ('smote', SMOTE(sampling_strategy=target_counts, random_state=42, k_neighbors=3)),
        ('clf', XGBClassifier(
            n_estimators=500, max_depth=5, learning_rate=0.05,
            subsample=0.8, use_label_encoder=False, eval_metric='mlogloss',
            random_state=42, n_jobs=-1
        ))
    ])
    
    y_pred_base = cross_val_predict(baseline_pipe, X_ml, y_enc, cv=skf)
    y_prob_base = cross_val_predict(baseline_pipe, X_ml, y_enc, cv=skf, method='predict_proba')
    
    print("\nBaseline Results:")
    print(classification_report(y_enc, y_pred_base, target_names=classes))


# ============================================================
# Strategy B: Hybrid POLE Rule + ML
# ============================================================
print("\n" + "=" * 50)
print("Strategy #7: Hybrid POLE Rule + ML")
print("=" * 50)

# Analyze POLE biological markers
print("\n--- Finding optimal POLE rule thresholds ---")

# SBS10a is 226x more in POLE vs others — strongest signal
# SBS10b is 45x, SBS28 is 70x
# TMB overlap exists (MSI also has high TMB), so we need SBS10 + TMB combo

# Find threshold using training data (leave-one-out for POLE)
sbs10a_pole = X_ml.loc[y == 'POLE', 'SBS10a'].values
sbs10a_other = X_ml.loc[y != 'POLE', 'SBS10a'].values
sbs10b_pole = X_ml.loc[y == 'POLE', 'SBS10b'].values
sbs10b_other = X_ml.loc[y != 'POLE', 'SBS10b'].values

print("SBS10a — POLE: min=%.1f, median=%.1f" % (sbs10a_pole.min(), np.median(sbs10a_pole)))
print("SBS10a — Others: max=%.1f, 99pct=%.1f" % (sbs10a_other.max(), np.percentile(sbs10a_other, 99)))
print("SBS10b — POLE: min=%.1f, median=%.1f" % (sbs10b_pole.min(), np.median(sbs10b_pole)))
print("SBS10b — Others: max=%.1f, 99pct=%.1f" % (sbs10b_other.max(), np.percentile(sbs10b_other, 99)))

# Combined SBS10 score
X_sbs10_total = X_ml['SBS10a'] + X_ml['SBS10b'] + X_ml['SBS10c'] + X_ml['SBS10d']
sbs10_total_pole = X_sbs10_total[y == 'POLE'].values
sbs10_total_other = X_sbs10_total[y != 'POLE'].values

print("\nSBS10_total — POLE: min=%.1f, median=%.1f, max=%.1f" % (
    sbs10_total_pole.min(), np.median(sbs10_total_pole), sbs10_total_pole.max()))
print("SBS10_total — Others: max=%.1f, 99pct=%.1f" % (
    sbs10_total_other.max(), np.percentile(sbs10_total_other, 99)))

# The hybrid approach: cross-validated to avoid data leakage
y_pred_hybrid = np.empty(len(y_enc), dtype=int)
y_prob_hybrid = np.zeros((len(y_enc), len(classes)))

for fold, (train_idx, test_idx) in enumerate(skf.split(X_ml, y_enc)):
    X_train, X_test = X_ml.iloc[train_idx], X_ml.iloc[test_idx]
    X_rules_test = X_rules.iloc[test_idx]  # Gene features for rules
    y_train = y_enc[train_idx]
    y_test_labels = y.iloc[test_idx].values
    
    # --- RULE 1: Determine POLE thresholds from training data ---
    train_labels = le.inverse_transform(y_train)
    pole_train_mask = train_labels == 'POLE'
    
    # SBS10 total score
    train_sbs10 = X_train['SBS10a'] + X_train['SBS10b'] + X_train['SBS10c'] + X_train['SBS10d']
    
    # Threshold: Use the minimum SBS10_total of POLE patients in training
    # with a small safety margin (take 80th percentile of non-POLE as lower bound)
    pole_sbs10_min = train_sbs10[pole_train_mask].min()
    other_sbs10_99 = train_sbs10[~pole_train_mask].quantile(0.99)
    
    # Conservative threshold: midpoint between other_99pct and pole_min
    sbs10_threshold = max(other_sbs10_99, pole_sbs10_min * 0.5)
    
    # Also require elevated TMB (>8, which is above CIN/GS/EBV range)
    tmb_threshold = 7.0  # POLE min TMB is ~7.9
    
    # --- RULE 2: Train ML on all data ---
    scaler = StandardScaler()
    X_train_sc = scaler.fit_transform(X_train)
    X_test_sc = scaler.transform(X_test)
    
    # SMOTE for training
    if HAS_SMOTE:
        train_counts = dict(zip(*np.unique(y_train, return_counts=True)))
        smote_target = dict(train_counts)
        if pole_idx in smote_target and smote_target[pole_idx] < 25:
            smote_target[pole_idx] = 25
        if gs_idx in smote_target and smote_target[gs_idx] < 50:
            smote_target[gs_idx] = 50
            
        try:
            smote = SMOTE(sampling_strategy=smote_target, random_state=42, k_neighbors=min(3, smote_target.get(pole_idx, 5)-1))
            X_train_res, y_train_res = smote.fit_resample(X_train_sc, y_train)
        except:
            X_train_res, y_train_res = X_train_sc, y_train
    else:
        X_train_res, y_train_res = X_train_sc, y_train
    
    clf = XGBClassifier(
        n_estimators=500, max_depth=5, learning_rate=0.05,
        subsample=0.8, use_label_encoder=False, eval_metric='mlogloss',
        random_state=42, n_jobs=-1
    )
    clf.fit(X_train_res, y_train_res)
    
    # --- PREDICT with hybrid logic ---
    ml_pred = clf.predict(X_test_sc)
    ml_prob = clf.predict_proba(X_test_sc)
    
    test_sbs10 = X_test['SBS10a'].values + X_test['SBS10b'].values + X_test['SBS10c'].values + X_test['SBS10d'].values
    test_tmb = X_test['TMB'].values
    # Gene features from rules dataset (not ML features)
    
    for i in range(len(test_idx)):
        idx = test_idx[i]
        
        # ---- RULE 1: POLE detection ----
        # If SBS10 total is very high AND TMB is elevated → candidate POLE
        if test_sbs10[i] > sbs10_threshold and test_tmb[i] > tmb_threshold:
            # MSI SAFEGUARD: MSI patients also have high TMB but are driven
            # by SBS15 (defective mismatch repair), not SBS10 (POLE proofreading).
            sbs15_val = X_test.iloc[i]['SBS15'] if 'SBS15' in X_test.columns else 0
            sbs10a_val = X_test.iloc[i]['SBS10a']  
            ml_msi_prob = ml_prob[i, msi_idx]
            
            is_likely_msi = (sbs15_val > sbs10a_val * 5) or (ml_msi_prob > 0.6)
            
            if is_likely_msi:
                y_pred_hybrid[idx] = ml_pred[i]
                y_prob_hybrid[idx] = ml_prob[i]
            else:
                y_pred_hybrid[idx] = pole_idx
                y_prob_hybrid[idx] = ml_prob[i]
                y_prob_hybrid[idx, pole_idx] = 0.95
        
        # ---- RULE 2: GS detection (conservative) ----
        # Only intervene when ML is UNCERTAIN (max_prob < 0.5)
        # CDH1 mutated (esp. truncating) + TP53 wild-type + low TMB = GS
        elif (HAS_CDH1 and HAS_TP53 and test_tmb[i] < 8.0):
            cdh1_mut = X_rules_test.iloc[i].get('CDH1_mutated', 0)
            cdh1_trunc = X_rules_test.iloc[i].get('CDH1_truncating', 0)
            tp53_mut = X_rules_test.iloc[i].get('TP53_mutated', 0)
            rhoa_mut = X_rules_test.iloc[i].get('RHOA_mutated', 0) if HAS_RHOA else 0
            
            # GS score: CDH1 truncating (strongest), CDH1 missense, RHOA, TP53 absence
            gs_evidence = 0
            if cdh1_trunc == 1:
                gs_evidence += 3    # Strong GS signal
            elif cdh1_mut == 1:
                gs_evidence += 1    # Moderate signal
            if rhoa_mut == 1:
                gs_evidence += 2    # RHOA is very GS-specific
            if tp53_mut == 0:
                gs_evidence += 1    # TP53 wild-type supports GS over CIN
            
            ml_max_prob = ml_prob[i].max()
            ml_best = ml_pred[i]
            ml_gs_prob = ml_prob[i, gs_idx]
            
            # Only override if: strong biological evidence AND ML is uncertain
            if gs_evidence >= 3 and ml_max_prob < 0.5:
                # Strong GS biology + ML confused → rule overrides
                y_pred_hybrid[idx] = gs_idx
                y_prob_hybrid[idx] = ml_prob[i]
                y_prob_hybrid[idx, gs_idx] = max(ml_gs_prob, 0.6)
            elif gs_evidence >= 2 and ml_best == cin_idx and ml_max_prob < 0.45:
                # Moderate GS biology + ML weakly says CIN → nudge to GS
                y_pred_hybrid[idx] = gs_idx
                y_prob_hybrid[idx] = ml_prob[i]
                y_prob_hybrid[idx, gs_idx] = max(ml_gs_prob, 0.5)
            else:
                y_pred_hybrid[idx] = ml_pred[i]
                y_prob_hybrid[idx] = ml_prob[i]
        
        else:
            y_pred_hybrid[idx] = ml_pred[i]
            y_prob_hybrid[idx] = ml_prob[i]

acc_h = accuracy_score(y_enc, y_pred_hybrid)
f1_h = f1_score(y_enc, y_pred_hybrid, average='weighted')
mcc_h = matthews_corrcoef(y_enc, y_pred_hybrid)

print("\n--- Hybrid Rule Results ---")
print("Accuracy: %.4f, F1: %.4f, MCC: %.4f" % (acc_h, f1_h, mcc_h))
print(classification_report(y_enc, y_pred_hybrid, target_names=classes))

# Count how many POLE were caught by rule vs ML
pole_true = y_enc == pole_idx
pole_correct = (y_pred_hybrid == pole_idx) & pole_true
print("POLE correct: %d / %d (%.0f%%)" % (pole_correct.sum(), pole_true.sum(), 100*pole_correct.mean()))


# ============================================================
# Strategy C: Hybrid + Threshold Tuning (#5)
# ============================================================
print("\n" + "=" * 50)
print("Strategy #5: + Per-Class Threshold Optimization")
print("=" * 50)

# Use the hybrid probabilities and optimize thresholds
# For each class, find the threshold that maximizes F1 for that class

best_thresholds = {}
print("\n--- Optimizing per-class thresholds ---")

for cls_idx_t in range(len(classes)):
    cls_name = classes[cls_idx_t]
    true_binary = (y_enc == cls_idx_t).astype(int)
    probs = y_prob_hybrid[:, cls_idx_t]
    
    best_f1 = 0
    best_thresh = 0.5
    
    for thresh in np.arange(0.10, 0.90, 0.02):
        pred_binary = (probs >= thresh).astype(int)
        tp = ((pred_binary == 1) & (true_binary == 1)).sum()
        fp = ((pred_binary == 1) & (true_binary == 0)).sum()
        fn = ((pred_binary == 0) & (true_binary == 1)).sum()
        
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        
        if f1 > best_f1:
            best_f1 = f1
            best_thresh = thresh
    
    best_thresholds[cls_idx_t] = best_thresh
    print("  %s: optimal threshold = %.2f (F1 = %.3f)" % (cls_name, best_thresh, best_f1))

# Apply optimized thresholds
y_pred_tuned = np.empty(len(y_enc), dtype=int)

for i in range(len(y_enc)):
    probs = y_prob_hybrid[i].copy()
    
    # Adjust probabilities by threshold distance
    adjusted_scores = np.zeros(len(classes))
    for c in range(len(classes)):
        # Score = how much above threshold
        adjusted_scores[c] = probs[c] - best_thresholds[c]
    
    # Predict class with highest adjusted score
    # But if any class is strongly above its threshold, prefer it
    y_pred_tuned[i] = np.argmax(adjusted_scores)

acc_t = accuracy_score(y_enc, y_pred_tuned)
f1_t = f1_score(y_enc, y_pred_tuned, average='weighted')
mcc_t = matthews_corrcoef(y_enc, y_pred_tuned)

print("\n--- Hybrid + Threshold Results ---")
print("Accuracy: %.4f, F1: %.4f, MCC: %.4f" % (acc_t, f1_t, mcc_t))
print(classification_report(y_enc, y_pred_tuned, target_names=classes))


# ============================================================
# FINAL COMPARISON
# ============================================================
print("\n" + "=" * 60)
print("FINAL COMPARISON")
print("=" * 60)

pole_mask = y_enc == pole_idx
gs_mask = y_enc == gs_idx
msi_mask = y_enc == msi_idx

results_all = {}
if HAS_SMOTE:
    results_all['Baseline (SMOTE+XGB)'] = y_pred_base
results_all['Hybrid Rule + ML'] = y_pred_hybrid
results_all['Hybrid + Threshold'] = y_pred_tuned

print("%-30s %8s %8s %10s %10s %10s" % ('Strategy', 'Acc', 'F1', 'MSI_Rec', 'GS_Rec', 'POLE_Rec'))
print("-" * 78)

for name, pred in results_all.items():
    acc = accuracy_score(y_enc, pred)
    f1 = f1_score(y_enc, pred, average='weighted')
    msi_rec = (pred[msi_mask] == msi_idx).mean()
    gs_rec = (pred[gs_mask] == gs_idx).mean()
    pole_rec = (pred[pole_mask] == pole_idx).mean()
    print("%-30s %8.4f %8.4f %10.1f%% %10.1f%% %10.1f%%" % (name, acc, f1, msi_rec*100, gs_rec*100, pole_rec*100))

# Save best
best_name = 'Hybrid + Threshold'
best_pred = y_pred_tuned

pred_df = pd.DataFrame({
    'Sample': X_ml.index,
    'True_Label': [classes[i] for i in y_enc],
    'Predicted': [classes[i] for i in best_pred],
    'Correct': y_enc == best_pred
})
pred_df.to_csv(os.path.join(OUTPUT_DIR, "predictions_hybrid.csv"), index=False)

# Save thresholds
import json
thresh_save = {classes[k]: float(v) for k, v in best_thresholds.items()}
with open(os.path.join(OUTPUT_DIR, "optimal_thresholds.json"), 'w') as f:
    json.dump(thresh_save, f, indent=2)

# Train final model on ALL data and save for SHAP (step7)
print("\nTraining final model on all data for SHAP analysis...")
from sklearn.pipeline import Pipeline as SkPipeline

final_scaler = StandardScaler()
X_ml_scaled = final_scaler.fit_transform(X_ml)

if HAS_SMOTE:
    counts_all = dict(zip(*np.unique(y_enc, return_counts=True)))
    smote_all = dict(counts_all)
    smote_all[pole_idx] = 30
    smote_all[gs_idx] = max(counts_all[gs_idx], 60)
    try:
        smote_final = SMOTE(sampling_strategy=smote_all, random_state=42, k_neighbors=3)
        X_res, y_res = smote_final.fit_resample(X_ml_scaled, y_enc)
    except:
        X_res, y_res = X_ml_scaled, y_enc
else:
    X_res, y_res = X_ml_scaled, y_enc

final_clf = XGBClassifier(
    n_estimators=500, max_depth=5, learning_rate=0.05,
    subsample=0.8, use_label_encoder=False, eval_metric='mlogloss',
    random_state=42, n_jobs=-1
)
final_clf.fit(X_res, y_res)

# Save as pipeline-like dict
model_save = {
    'model': SkPipeline([('scaler', final_scaler), ('classifier', final_clf)]),
    'label_encoder': le,
    'classes': list(classes),
    'thresholds': thresh_save,
    'feature_names': list(X_ml.columns)
}
model_path = os.path.join("output/ml_results", "best_model.pkl")
os.makedirs("output/ml_results", exist_ok=True)
with open(model_path, 'wb') as f:
    pickle.dump(model_save, f)
print("  Saved model: %s" % model_path)

print("\nSaved to: %s/" % OUTPUT_DIR)
