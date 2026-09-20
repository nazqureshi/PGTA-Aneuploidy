#!/usr/bin/env python3
"""
step3_generate_gt_meta_no_totalembryos.py
Generate gt_meta files WITHOUT total_embryos covariate
Assumes data is already MAF/INFO filtered
"""
import pandas as pd
import numpy as np
import os
import subprocess
import sys
import gc
import time

def process_chromosome_step3(chromosome):
    """Process chromosome data for Step 3"""
    print(f"Processing chromosome {chromosome} - Step 3", flush=True)
    print(f"Working directory: {os.getcwd()}", flush=True)
    
    # File paths (assumes pre-filtered data)
    gt_file = f"gt_5616pts_ligated_merged_chr{chromosome}_maffilt.txt"
    output_file = f"gt_meta_no_totalembryos_chr{chromosome}.csv"
    meta_file = "merged_phenotype_single_seqdate.csv"
    bcf_file = f"5616pts_ligated_merged_chr{chromosome}_maffilt.bcf"
    
    # Check genotype file
    if not os.path.exists(gt_file):
        print(f"ERROR: Genotype file not found: {gt_file}", flush=True)
        return False
    
    file_size_gb = os.path.getsize(gt_file) / (1024**3)
    print(f"Found genotype file: {gt_file} ({file_size_gb:.1f} GB)", flush=True)
    
    # Load metadata
    print("Loading metadata...", flush=True)
    try:
        meta_seq = pd.read_csv(meta_file)
        meta_seq["patient_ID"] = meta_seq["patient_ID"].astype(str)
        meta_ids = list(meta_seq["patient_ID"])
        print(f"   Loaded {len(meta_seq)} samples from metadata", flush=True)
        print(f"   Metadata columns: {list(meta_seq.columns)}", flush=True)
        
        # Verify required columns (NO total_embryos, NO paternal_age)
        required_cols = [
            'patient_ID', 
            'embryo_count',
            'aneuploid_egg_num',
            'euploid_egg_num',
            'mean_maternal_age',
            'seq_date_numeric',
            'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8'
        ]
        
        missing_cols = [col for col in required_cols if col not in meta_seq.columns]
        if missing_cols:
            print(f"ERROR: Missing required columns: {missing_cols}", flush=True)
            print(f"Available columns: {list(meta_seq.columns)}", flush=True)
            return False
        
        print(f"   ✓ All required columns present", flush=True)
        print(f"   ✓ NO total_embryos covariate", flush=True)
        print(f"   ✓ NO paternal_age covariate", flush=True)
        
    except Exception as e:
        print(f"ERROR loading metadata: {e}", flush=True)
        import traceback
        traceback.print_exc()
        return False
    
    # Get sample names from BCF
    print("Getting sample order from BCF...", flush=True)
    try:
        result = subprocess.run(['bcftools', 'query', '-l', bcf_file], 
                               capture_output=True, text=True, check=True)
        sample_paths = result.stdout.strip().split('\n')
        
        sample_names = []
        for path in sample_paths:
            filename = path.split('/')[-1]
            sample_id = str(filename.split('_')[0])
            sample_names.append(sample_id)
        
        print(f"   Found {len(sample_names)} samples in BCF", flush=True)
        
    except Exception as e:
        print(f"ERROR getting sample names: {e}", flush=True)
        return False
    
    # Process genotype file
    print("Processing genotype matrix...", flush=True)
    
    try:
        # Collect variant IDs
        print("   Collecting variant IDs...", flush=True)
        variant_ids = []
        line_count = 0
        
        with open(gt_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t', 1)
                if len(parts) > 0:
                    variant_ids.append(parts[0])
                    line_count += 1
                    if line_count % 100000 == 0:
                        print(f"      Read {line_count:,} variants", flush=True)
        
        print(f"   Total variants: {len(variant_ids):,}", flush=True)
        
        # Prepare ordered metadata
        print("   Preparing ordered metadata...", flush=True)
        sample_to_idx = {name: idx for idx, name in enumerate(sample_names)}
        
        common_samples = [s for s in sample_names if s in meta_ids]
        common_indices = [sample_to_idx[s] for s in common_samples]
        
        print(f"   Found {len(common_samples)} overlapping samples", flush=True)
        
        if len(common_samples) == 0:
            print("ERROR: No overlapping samples found!", flush=True)
            return False
        
        # Order metadata
        meta_ordered = []
        for sample_id in common_samples:
            sample_data = meta_seq[meta_seq["patient_ID"] == sample_id]
            if len(sample_data) > 0:
                meta_ordered.append(sample_data.iloc[0])
        
        metadata_df = pd.DataFrame(meta_ordered)
        metadata_df.reset_index(drop=True, inplace=True)
        
        # Select columns (NO total_embryos)
        metadata_columns_ordered = [
            'patient_ID',
            'embryo_count',
            'aneuploid_egg_num',
            'euploid_egg_num',
            'mean_maternal_age',
            'aneuploidy_rate',
            'seq_date_numeric',
            'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8'
        ]
        
        metadata_df = metadata_df[metadata_columns_ordered]
        
        print(f"   Metadata columns in output: {list(metadata_df.columns)}", flush=True)
        print(f"   Number of metadata columns: {len(metadata_df.columns)}", flush=True)
        
        # Write header
        print("   Writing output file header...", flush=True)
        header = list(metadata_df.columns) + variant_ids
        
        with open(output_file, 'w') as f:
            f.write(','.join(map(str, header)) + '\n')
        
        # Process in batches
        print("   Processing genotype data in batches...", flush=True)
        
        samples_per_batch = 100
        num_sample_batches = (len(common_samples) + samples_per_batch - 1) // samples_per_batch
        
        for batch_idx in range(num_sample_batches):
            start_idx = batch_idx * samples_per_batch
            end_idx = min((batch_idx + 1) * samples_per_batch, len(common_samples))
            batch_samples = common_samples[start_idx:end_idx]
            batch_indices = common_indices[start_idx:end_idx]
            
            print(f"   Processing samples {start_idx+1}-{end_idx} of {len(common_samples)}...", flush=True)
            
            batch_data = np.zeros((len(batch_samples), len(variant_ids)), dtype=np.float32)
            
            with open(gt_file, 'r') as f:
                var_idx = 0
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) > 1:
                        for i, sample_idx in enumerate(batch_indices):
                            if sample_idx + 1 < len(parts):
                                try:
                                    batch_data[i, var_idx] = float(parts[sample_idx + 1])
                                except (ValueError, IndexError):
                                    batch_data[i, var_idx] = np.nan
                    var_idx += 1
                    
                    if var_idx % 100000 == 0:
                        print(f"      Processed {var_idx:,} variants for current batch", flush=True)
            
            print(f"   Writing batch {batch_idx+1}/{num_sample_batches} to file...", flush=True)
            
            with open(output_file, 'a') as f:
                for i in range(len(batch_samples)):
                    meta_values = metadata_df.iloc[start_idx + i].values
                    row_data = list(meta_values) + list(batch_data[i, :])
                    f.write(','.join(map(str, row_data)) + '\n')
            
            del batch_data
            gc.collect()
            
            print(f"   Batch {batch_idx+1}/{num_sample_batches} complete", flush=True)
        
        print(f"Successfully created: {output_file}", flush=True)
        
    except Exception as e:
        print(f"ERROR during processing: {e}", flush=True)
        import traceback
        traceback.print_exc()
        return False
    
    # Verification
    print("\nVerification...", flush=True)
    try:
        file_size_gb = os.path.getsize(output_file) / (1024**3)
        
        with open(output_file, 'r') as f:
            header_line = f.readline().strip()
        
        num_cols = len(header_line.split(','))
        
        print(f"   Output file: {output_file}", flush=True)
        print(f"   File size: {file_size_gb:.1f} GB", flush=True)
        print(f"   Number of columns: {num_cols:,}", flush=True)
        print(f"   Expected: {len(metadata_df.columns)} + {len(variant_ids):,} = {len(metadata_df.columns) + len(variant_ids):,}", flush=True)
        
        if num_cols == len(metadata_df.columns) + len(variant_ids):
            print("   ✓ Column count matches", flush=True)
        
        header_cols = header_line.split(',')
        metadata_cols_in_output = header_cols[:len(metadata_df.columns)]
        print(f"   Metadata columns: {metadata_cols_in_output}", flush=True)
        
        if 'total_embryos' not in metadata_cols_in_output:
            print("   ✓ total_embryos correctly excluded", flush=True)
            
        print(f"SUCCESS: Chromosome {chromosome} processed!", flush=True)
        return True
        
    except Exception as e:
        print(f"ERROR in verification: {e}", flush=True)
        return False

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 step3_generate_gt_meta_no_totalembryos.py <chromosome_number>")
        sys.exit(1)
    
    try:
        chromosome = int(sys.argv[1])
        if chromosome < 1 or chromosome > 22:
            print("Chromosome must be between 1 and 22")
            sys.exit(1)
    except ValueError:
        print("Chromosome must be an integer")
        sys.exit(1)
    
    # Cluster path root -- update to match your own environment.
    work_dir = "/path/to/project/imputation/outputs/ligate/maffilt"
    os.chdir(work_dir)
    
    print("=" * 60)
    print(f"STEP 3: Generate GT Meta File (No total_embryos)")
    print(f"Chromosome: {chromosome}")
    print("=" * 60)
    
    start_time = time.time()
    success = process_chromosome_step3(chromosome)
    elapsed_time = time.time() - start_time
    
    print(f"\nProcessing time: {elapsed_time/60:.1f} minutes", flush=True)
    
    if success:
        print(f"SUCCESS: Completed for chromosome {chromosome}", flush=True)
        sys.exit(0)
    else:
        print(f"FAILED: chromosome {chromosome}", flush=True)
        sys.exit(1)

if __name__ == "__main__":
    main()