#!/usr/bin/env python3
"""
Step 1: Process MAF Files — SBS96 Matrix + Gene-Level Features
==============================================================
MERGED: build_sbs96_matrix.py + extract_gene_features.py

Reads 434 GDC masked somatic mutation MAF files in a SINGLE PASS and
extracts BOTH the 96-channel SBS trinucleotide matrix AND gene-level
mutation features (CDH1, RHOA, TP53, Ti/Tv ratio, etc.)

Input:  maf_files/*.maf (434 GDC masked somatic mutation MAF files)
Output: output/sbs96_matrix.csv          (434 samples × 96 mutation types)
        output/gene_mutation_features.csv (434 samples × 53 gene features)
"""

import os, sys
import pandas as pd
import numpy as np
from collections import Counter, OrderedDict
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# CONFIGURATION
# ============================================================
MAF_DIR = "maf_files"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# SBS96 Constants
# ============================================================
PURINES = {'A', 'G'}
PYRIMIDINES = {'C', 'T'}
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
BASES = ['A', 'C', 'G', 'T']
SUB_TYPES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']

def get_96_channels():
    """Generate the 96 SBS channel labels in standard order."""
    channels = []
    for sub in SUB_TYPES:
        ref = sub[0]
        for five_prime in BASES:
            for three_prime in BASES:
                label = f"{five_prime}[{sub}]{three_prime}"
                channels.append(label)
    return channels

SBS96_CHANNELS = get_96_channels()

# ============================================================
# Gene Feature Constants
# ============================================================
GS_MARKER_GENES = ['CDH1', 'RHOA', 'ARID1A', 'CDH2', 'CLDN18']
CIN_MARKER_GENES = ['TP53', 'ERBB2', 'VEGFA', 'CCNE1', 'MYC']
OTHER_GENES = ['PIK3CA', 'KRAS', 'PTEN', 'APC', 'SMAD4', 
               'CTNNB1', 'FBXW7', 'ERBB3', 'TGFBR2', 'RNF43']
ALL_GENES = GS_MARKER_GENES + CIN_MARKER_GENES + OTHER_GENES

SILENT_TYPES = {'Silent', 'Intron', "3'UTR", "5'UTR", "3'Flank", "5'Flank", 'IGR', 'RNA'}
TRUNCATING_TYPES = {'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins',
                    'Splice_Site', 'Splice_Region'}
TRANSITIONS = {'AG', 'GA', 'CT', 'TC'}


# ============================================================
# CORE: Single-pass MAF processing
# ============================================================
def classify_mutation(ref, alt, context):
    """Classify a single SNP into one of 96 SBS channels."""
    if ref not in 'ACGT' or alt not in 'ACGT' or ref == alt:
        return None
    if context is None or len(context) < 3:
        return None
    
    mid = len(context) // 2
    five_prime = context[mid - 1].upper()
    three_prime = context[mid + 1].upper()
    
    if five_prime not in 'ACGT' or three_prime not in 'ACGT':
        return None
    
    if ref in PURINES:
        ref = COMPLEMENT[ref]
        alt = COMPLEMENT[alt]
        five_prime, three_prime = COMPLEMENT[three_prime], COMPLEMENT[five_prime]
    
    sub = f"{ref}>{alt}"
    return f"{five_prime}[{sub}]{three_prime}"


def process_single_maf(filepath):
    """
    Process a single MAF file and extract BOTH SBS96 counts AND gene features.
    Returns: (sample_barcode, sbs96_counts, gene_features)
    """
    try:
        cols_needed = ['Hugo_Symbol', 'Variant_Classification', 'Variant_Type',
                       'Tumor_Sample_Barcode', 'Reference_Allele', 'Tumor_Seq_Allele2',
                       'CONTEXT']
        df = pd.read_csv(filepath, sep='\t', comment='#',
                         usecols=lambda c: c in cols_needed,
                         low_memory=False)
    except Exception as e:
        return None, Counter(), {}
    
    if len(df) == 0:
        return None, Counter(), {}
    
    barcode = df['Tumor_Sample_Barcode'].iloc[0]
    
    # ---- PART A: SBS96 Channel Counts ----
    sbs96_counts = Counter()
    snps = df[df['Variant_Type'] == 'SNP']
    
    if 'CONTEXT' in df.columns:
        for _, row in snps.iterrows():
            ref = str(row['Reference_Allele']).upper()
            alt = str(row['Tumor_Seq_Allele2']).upper()
            context = str(row.get('CONTEXT', ''))
            channel = classify_mutation(ref, alt, context)
            if channel and channel in SBS96_CHANNELS:
                sbs96_counts[channel] += 1
    
    # ---- PART B: Gene-Level Features ----
    features = {}
    
    # Non-silent mutations
    nonsil = df[~df['Variant_Classification'].isin(SILENT_TYPES)]
    mutated_genes = set(nonsil['Hugo_Symbol'].unique())
    
    # Binary mutation status for each gene
    for gene in ALL_GENES:
        features[gene + '_mutated'] = 1 if gene in mutated_genes else 0
    
    # Mutation counts per gene
    gene_counts = nonsil[nonsil['Hugo_Symbol'].isin(ALL_GENES)]['Hugo_Symbol'].value_counts()
    for gene in ALL_GENES:
        features[gene + '_mut_count'] = int(gene_counts.get(gene, 0))
    
    # CDH1 truncating mutations (critical for GS subtype)
    cdh1_muts = nonsil[nonsil['Hugo_Symbol'] == 'CDH1']
    features['CDH1_truncating'] = 1 if any(cdh1_muts['Variant_Classification'].isin(TRUNCATING_TYPES)) else 0
    features['CDH1_missense'] = 1 if 'Missense_Mutation' in cdh1_muts['Variant_Classification'].values else 0
    
    # Variant type counts
    snp_count = (df['Variant_Type'] == 'SNP').sum()
    ins_count = (df['Variant_Type'] == 'INS').sum()
    del_count = (df['Variant_Type'] == 'DEL').sum()
    features['total_snps'] = int(snp_count)
    features['total_indels'] = int(ins_count + del_count)
    features['indel_fraction'] = float(features['total_indels'] / max(1, snp_count + ins_count + del_count))
    
    # Ti/Tv ratio
    ti_count = tv_count = 0
    if len(snps) > 0:
        for _, row in snps.iterrows():
            ref = str(row['Reference_Allele']).upper()
            alt = str(row['Tumor_Seq_Allele2']).upper()
            if ref in 'ACGT' and alt in 'ACGT' and ref != alt:
                if ref + alt in TRANSITIONS:
                    ti_count += 1
                else:
                    tv_count += 1
    
    features['Ti_count'] = ti_count
    features['Tv_count'] = tv_count
    features['TiTv_ratio'] = float(ti_count / max(1, tv_count))
    
    # Variant classification distribution
    vc_counts = nonsil['Variant_Classification'].value_counts()
    features['missense_count'] = int(vc_counts.get('Missense_Mutation', 0))
    features['nonsense_count'] = int(vc_counts.get('Nonsense_Mutation', 0))
    features['frameshift_count'] = int(vc_counts.get('Frame_Shift_Del', 0) + vc_counts.get('Frame_Shift_Ins', 0))
    features['splice_count'] = int(vc_counts.get('Splice_Site', 0) + vc_counts.get('Splice_Region', 0))
    
    total_nonsil = len(nonsil)
    truncating = features['nonsense_count'] + features['frameshift_count'] + features['splice_count']
    features['truncating_fraction'] = float(truncating / max(1, total_nonsil))
    
    return barcode, sbs96_counts, features


# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 60)
    print("Step 1: Process MAF Files (SBS96 + Gene Features)")
    print("=" * 60)
    
    maf_files = sorted([f for f in os.listdir(MAF_DIR) if f.endswith('.maf')])
    print(f"\nFound {len(maf_files)} MAF files")
    
    if len(maf_files) == 0:
        print("No MAF files found!")
        sys.exit(1)
    
    # Single pass over all MAF files
    sbs96_results = {}
    gene_results = {}
    skipped = 0
    
    for filename in tqdm(maf_files, desc="Processing MAF files"):
        filepath = os.path.join(MAF_DIR, filename)
        barcode, sbs96_counts, gene_features = process_single_maf(filepath)
        
        if barcode is None:
            skipped += 1
            continue
        
        sbs96_results[barcode] = sbs96_counts
        gene_results[barcode] = gene_features
    
    print(f"\nProcessed {len(sbs96_results)} samples, skipped {skipped}")
    
    # ---- Build SBS96 Matrix ----
    print(f"\n{'=' * 40}")
    print("Building SBS96 Matrix...")
    
    sbs96_matrix = pd.DataFrame.from_dict(sbs96_results, orient='index')
    for ch in SBS96_CHANNELS:
        if ch not in sbs96_matrix.columns:
            sbs96_matrix[ch] = 0
    sbs96_matrix = sbs96_matrix[SBS96_CHANNELS].fillna(0).astype(int)
    sbs96_matrix.index.name = 'Sample'
    
    sbs96_matrix['Total'] = sbs96_matrix.sum(axis=1)
    sbs96_matrix = sbs96_matrix.sort_values('Total', ascending=False)
    total_col = sbs96_matrix.pop('Total')
    
    sbs96_path = os.path.join(OUTPUT_DIR, "sbs96_matrix.csv")
    sbs96_matrix.to_csv(sbs96_path)
    
    print(f"  Samples:          {sbs96_matrix.shape[0]}")
    print(f"  Mutation types:   {sbs96_matrix.shape[1]} (should be 96)")
    print(f"  Total mutations:  {sbs96_matrix.values.sum():,}")
    print(f"  Median per sample: {int(total_col.median()):,}")
    print(f"  Saved: {sbs96_path}")
    
    # ---- Build Gene Feature Matrix ----
    print(f"\n{'=' * 40}")
    print("Building Gene Feature Matrix...")
    
    gene_df = pd.DataFrame.from_dict(gene_results, orient='index')
    gene_df.index.name = 'Sample'
    
    gene_path = os.path.join(OUTPUT_DIR, "gene_mutation_features.csv")
    gene_df.to_csv(gene_path)
    
    print(f"  Samples:  {gene_df.shape[0]}")
    print(f"  Features: {gene_df.shape[1]}")
    print(f"  Saved: {gene_path}")
    
    # ---- Validation: Gene Mutation Rates by Subtype ----
    labels_path = os.path.join(OUTPUT_DIR, "ml_labels.csv")
    if os.path.exists(labels_path):
        y_df = pd.read_csv(labels_path, index_col=0)
        common = gene_df.index.intersection(y_df.index)
        
        if len(common) > 0:
            gene_sub = gene_df.loc[common]
            y_sub = y_df.loc[common, 'molecular_subtype']
            
            print("\n=== Gene Mutation Rates by Subtype ===")
            for gene in GS_MARKER_GENES + CIN_MARKER_GENES[:3]:
                col = gene + '_mutated'
                if col in gene_sub.columns:
                    print(f"\n  {gene}:")
                    for sub in ['CIN', 'EBV', 'GS', 'MSI', 'POLE']:
                        mask = y_sub == sub
                        if mask.sum() > 0:
                            rate = gene_sub.loc[mask, col].mean() * 100
                            print(f"    {sub}: {rate:.1f}% ({int(gene_sub.loc[mask, col].sum())}/{mask.sum()})")
            
            print("\n=== Ti/Tv Ratio by Subtype ===")
            for sub in ['CIN', 'EBV', 'GS', 'MSI', 'POLE']:
                mask = y_sub == sub
                if mask.sum() > 0:
                    ratio = gene_sub.loc[mask, 'TiTv_ratio'].median()
                    print(f"  {sub}: median Ti/Tv = {ratio:.2f}")
    
    print(f"\n{'=' * 60}")
    print("Step 1 complete!")
    print("=" * 60)
    
    return sbs96_matrix, gene_df


if __name__ == "__main__":
    main()
