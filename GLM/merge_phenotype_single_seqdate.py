#!/usr/bin/env python3
"""
merge_phenotype_single_seqdate.py
Aggregate embryo-level data to patient-level with single seq_date covariate
NO PATERNAL AGE
FIXED: Use correct PC file format (indivID column)
"""
import pandas as pd
import numpy as np

# Cluster path root -- update to match your own environment.
PROJECT_ROOT = "/path/to/project"

# Input files
EMBRYO_FILE = f"{PROJECT_ROOT}/imputation/outputs/ligate/maffilt/dataset_ONE_with_seqdate.csv"
PC_FILE = f"{PROJECT_ROOT}/ancestry_inference/scripts/laser.SeqPC.coord"
OUTPUT_FILE = f"{PROJECT_ROOT}/imputation/outputs/ligate/maffilt/merged_phenotype_single_seqdate.csv"

print("="*80)
print("MERGING PHENOTYPE WITH SINGLE SEQ DATE COVARIATE")
print("="*80)

# Load embryo-level data
print(f"\nLoading embryo data: {EMBRYO_FILE}")
embryo_df = pd.read_csv(EMBRYO_FILE)
print(f"  {len(embryo_df)} embryos")

# Ensure patient_ID is string type
embryo_df['patient_ID'] = embryo_df['patient_ID'].astype(str)

# Create binary aneuploidy indicator
embryo_df['is_aneuploid'] = (embryo_df['CCS_interpretation'] == 'Abnormal').astype(int)

# Aggregate to patient level (NO PATERNAL AGE)
print("\nAggregating to patient level...")
patient_df = embryo_df.groupby('patient_ID').agg({
    'is_aneuploid': 'sum',
    'embryo_ID': 'count',
    'maternal_age': 'mean',
    'seq_date_numeric': 'first',
    'sequencing_date': 'first',
    'file_name': 'first',
    'path': 'first'
}).reset_index()

# Rename columns
patient_df.rename(columns={
    'is_aneuploid': 'aneuploid_egg_num',
    'embryo_ID': 'embryo_count',
    'maternal_age': 'mean_maternal_age',
    'file_name': 'first_file_name',
    'path': 'first_path'
}, inplace=True)

# Calculate euploid count and aneuploidy rate
patient_df['euploid_egg_num'] = patient_df['embryo_count'] - patient_df['aneuploid_egg_num']
patient_df['aneuploidy_rate'] = patient_df['aneuploid_egg_num'] / patient_df['embryo_count']

print(f"  {len(patient_df)} patients")

# Load PCs (CORRECT FORMAT)
print(f"\nLoading PCs: {PC_FILE}")
pc_data = pd.read_csv(PC_FILE, sep=r'\s+')  # Use whitespace separator
print(f"  {len(pc_data)} samples in PC file")
print(f"  PC file columns: {list(pc_data.columns)}")

# Extract sample ID from indivID column
# Format: "1468_merged.sorted" -> "1468"
pc_data['sample_id'] = pc_data['indivID'].str.replace('_merged.sorted', '', regex=False)
pc_data['sample_id'] = pc_data['sample_id'].astype(str)

# Select PC columns (first 8)
pc_cols = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8']
pc_subset = pc_data[['sample_id'] + pc_cols].copy()

print(f"  Extracted {len(pc_subset)} patient IDs from PC data")
print(f"  Sample PC IDs: {pc_subset['sample_id'].head(10).tolist()}")

# Debug: Check patient_ID formats
print(f"\nDebug - patient_df patient_IDs (first 10):")
print(patient_df['patient_ID'].head(10).tolist())
print(f"\nDebug - pc_subset sample_ids (first 10):")
print(pc_subset['sample_id'].head(10).tolist())

# Merge
print("\nMerging...")
merged_df = patient_df.merge(pc_subset, 
                             left_on='patient_ID', 
                             right_on='sample_id', 
                             how='inner')

# Drop redundant sample_id column
if 'sample_id' in merged_df.columns:
    merged_df = merged_df.drop('sample_id', axis=1)

print(f"  {len(merged_df)} patients after merge")

# Check what was lost
n_lost = len(patient_df) - len(merged_df)
if n_lost > 0:
    print(f"  WARNING: Lost {n_lost} patients in merge")
    
    # Show examples of IDs that didn't match
    pheno_ids = set(patient_df['patient_ID'])
    pc_ids = set(pc_subset['sample_id'])
    
    missing_in_pc = pheno_ids - pc_ids
    missing_in_pheno = pc_ids - pheno_ids
    
    if missing_in_pc:
        print(f"\n  {len(missing_in_pc)} patient IDs in phenotype not found in PC data")
        print(f"  Examples: {list(missing_in_pc)[:5]}")
    
    if missing_in_pheno:
        print(f"\n  {len(missing_in_pheno)} IDs in PC data not found in phenotype")
        print(f"  Examples: {list(missing_in_pheno)[:5]}")

# Select final columns (NO PATERNAL AGE)
final_columns = [
    'patient_ID',
    'embryo_count',
    'aneuploid_egg_num',
    'euploid_egg_num',
    'mean_maternal_age',
    'aneuploidy_rate',
    'sequencing_date',
    'seq_date_numeric',
    'first_file_name',
    'first_path',
    'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8'
]

output_df = merged_df[final_columns]

# Save
output_df.to_csv(OUTPUT_FILE, index=False)

print(f"\nSaved: {OUTPUT_FILE}")
print(f"\nFinal dataset:")
print(f"  Patients: {len(output_df)}")
print(f"  Columns: {len(output_df.columns)}")

print("\nColumn names:")
print(list(output_df.columns))

if len(output_df) > 0:
    print("\nSequencing date distribution:")
    print(output_df.groupby('sequencing_date').size().sort_index())
    
    print("\nSeq date numeric stats:")
    print(output_df['seq_date_numeric'].describe())
    
    print("\nSample of data:")
    print(output_df.head())
else:
    print("\nERROR: No data in final output!")

print("\n" + "="*80)
print("COMPLETE")
print("="*80)