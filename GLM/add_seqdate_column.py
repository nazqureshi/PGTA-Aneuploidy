#!/usr/bin/env python3
"""
add_seqdate_column.py
Add sequencing date as a single numeric column to existing phenotype file
FIXED: Correctly extract month from YYMMDD format
"""
import pandas as pd
import numpy as np
import re

# Cluster path root -- update to match your own environment.
PROJECT_ROOT = "/path/to/project"

# Input/output files
INPUT_FILE = f"{PROJECT_ROOT}/imputation/outputs/ligate/maffilt/dataset_ONE.csv"
OUTPUT_FILE = f"{PROJECT_ROOT}/imputation/outputs/ligate/maffilt/dataset_ONE_with_seqdate.csv"

print("="*80)
print("ADDING SEQUENCING DATE COLUMN")
print("="*80)

# Load data
print(f"\nLoading: {INPUT_FILE}")
df = pd.read_csv(INPUT_FILE)
print(f"  {len(df)} rows (embryos)")

# Extract sequencing date from path
def extract_seq_date(path):
    """
    Extract YYYY_MM from path like:
    'K:/MiSeq4/MiSeqOutput/2021/October_2021/211004_...'
                            ^^^^             ^^^^^^
                            YYYY           YYMMDD
    
    The date code format is YYMMDD where:
    - YY = 2-digit year (21 = 2021)
    - MM = 2-digit month (10 = October)
    - DD = 2-digit day (04 = 4th)
    
    Returns: '2021_10' (year_month)
    """
    if pd.isna(path):
        return None
    
    # Match: YYYY/Month_YYYY/YYMMDD_...
    match = re.search(r'(\d{4})/\w+_\d{4}/\d{2}(\d{2})\d{2}_', path)
    #                   ^^^^  capture year    ^^  capture month (middle 2 digits of YYMMDD)
    
    if match:
        year = match.group(1)      # e.g., "2021"
        month = match.group(2)     # e.g., "10" (from 211004)
        return f"{year}_{month}"
    
    return None

print("\nExtracting sequencing dates from paths...")
df['sequencing_date'] = df['path'].apply(extract_seq_date)

# Check how many were extracted
n_extracted = df['sequencing_date'].notna().sum()
print(f"  Successfully extracted: {n_extracted} / {len(df)} ({n_extracted/len(df)*100:.1f}%)")

# Show unique dates
print(f"\nUnique sequencing dates found:")
date_counts = df['sequencing_date'].value_counts().sort_index()
print(date_counts)

# Convert to numeric (months since January 2021)
def date_to_numeric(date_str):
    """
    Convert YYYY_MM to months since 2021-01
    Example: 2021_10 -> 10 (October 2021 = 10 months since Jan 2021)
             2022_03 -> 15 (March 2022 = 15 months since Jan 2021)
    """
    if pd.isna(date_str):
        return np.nan
    
    year, month = date_str.split('_')
    # Months since January 2021 (baseline = 1)
    months_since_start = (int(year) - 2021) * 12 + int(month)
    return months_since_start

print("\nConverting to numeric (months since Jan 2021)...")
df['seq_date_numeric'] = df['sequencing_date'].apply(date_to_numeric)

# Show distribution
print(f"\nSequencing date numeric distribution:")
print(f"  Min: {df['seq_date_numeric'].min()} ({df[df['seq_date_numeric'] == df['seq_date_numeric'].min()]['sequencing_date'].iloc[0]})")
print(f"  Max: {df['seq_date_numeric'].max()} ({df[df['seq_date_numeric'] == df['seq_date_numeric'].max()]['sequencing_date'].iloc[0]})")
print(f"  Mean: {df['seq_date_numeric'].mean():.2f}")
print(f"  Median: {df['seq_date_numeric'].median():.2f}")
print(f"  Std: {df['seq_date_numeric'].std():.2f}")

# Show example rows
print("\nExample rows with new columns:")
print(df[['patient_ID', 'embryo_ID', 'path', 'sequencing_date', 'seq_date_numeric']].head(10).to_string(index=False))

# Verify the fix worked
print("\nVerification - checking October 2021 samples:")
october_samples = df[df['path'].str.contains('October_2021', na=False)].head(3)
print(october_samples[['path', 'sequencing_date', 'seq_date_numeric']].to_string(index=False))

# Save
df.to_csv(OUTPUT_FILE, index=False)
print(f"\nSaved: {OUTPUT_FILE}")

print("\n" + "="*80)
print("COMPLETE")
print("="*80)
print(f"\nNew columns added:")
print(f"  - sequencing_date (YYYY_MM format)")
print(f"  - seq_date_numeric (months since Jan 2021)")