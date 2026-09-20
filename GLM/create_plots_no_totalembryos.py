#!/usr/bin/env python3
"""
create_plots_no_totalembryos.py
Generate comprehensive plots for no total_embryos model
Pre-filtered data (MAF≥0.05, INFO≥0.20)
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import os
import gc

# ============================================================================
# CONFIGURATION
# ============================================================================

# Cluster path root -- update to match your own environment.
PROJECT_ROOT = "/path/to/project/imputation/outputs/ligate/maffilt/asso_no_totalembryos"

GWAS_FILE = f"{PROJECT_ROOT}/all_chromosomes_combined.csv"
OUTPUT_DIR = f"{PROJECT_ROOT}/plots"

# ============================================================================
# FUNCTIONS
# ============================================================================

def parse_variant_id(variant_id):
    """Extract chromosome and position from variant ID"""
    try:
        if ':' in variant_id:
            parts = variant_id.split(':')
        elif '_' in variant_id:
            parts = variant_id.split('_')
        else:
            return None, None
        
        chrom = parts[0].replace('chr', '').replace('Chr', '')
        
        try:
            chrom = int(chrom)
        except ValueError:
            return None, None
        
        pos = int(parts[1])
        return chrom, pos
        
    except (IndexError, ValueError):
        return None, None

def load_gwas_with_positions(gwas_file):
    """Load GWAS results and parse chr/pos"""
    print("\nLoading GWAS results...")
    
    chunks = []
    chunk_size = 1000000
    total_rows = 0
    
    for chunk in pd.read_csv(gwas_file, chunksize=chunk_size):
        total_rows += len(chunk)
        
        chunk['variant'] = chunk['variant'].str.strip()
        
        # Parse chr/pos
        parsed = chunk['variant'].apply(parse_variant_id)
        chunk['chromosome'] = parsed.apply(lambda x: x[0])
        chunk['position'] = parsed.apply(lambda x: x[1])
        
        # Remove unparseable
        chunk = chunk.dropna(subset=['chromosome', 'position'])
        chunk['chromosome'] = chunk['chromosome'].astype(int)
        chunk['position'] = chunk['position'].astype(int)
        
        # Keep only autosomes
        chunk = chunk[chunk['chromosome'].between(1, 22)]
        
        if len(chunk) > 0:
            chunks.append(chunk)
        
        if total_rows % 10000000 == 0:
            print(f"  Processed {total_rows:,} rows...")
        
        del chunk
        gc.collect()
    
    print(f"  Combining chunks...")
    df = pd.concat(chunks, ignore_index=True)
    
    # Calculate -log10(p) if not present
    if 'neg_log10_p' not in df.columns:
        df['neg_log10_p'] = -np.log10(df['p'])
    
    print(f"  Final dataset: {len(df):,} variants")
    
    return df

def calculate_lambda(gwas_df):
    """Calculate genomic inflation factor"""
    t_values = gwas_df['t'].dropna().values
    
    if len(t_values) == 0:
        return np.nan
    
    chi2_values = t_values ** 2
    median_chi2 = np.median(chi2_values)
    lambda_gc = median_chi2 / 0.4549364
    
    return lambda_gc

def create_qq_plot(gwas_df, lambda_val, output_file):
    """Create QQ plot"""
    print("\nCreating QQ plot...")
    
    p_values = gwas_df['p'].dropna().values
    
    if len(p_values) == 0:
        print("  No p-values available")
        return
    
    n = len(p_values)
    expected_p = np.arange(1, n + 1) / (n + 1)
    observed_p = np.sort(p_values)
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    ax.scatter(-np.log10(expected_p), -np.log10(observed_p), 
              s=3, alpha=0.5, c='darkblue')
    
    max_val = max(-np.log10(expected_p[0]), -np.log10(observed_p[0]))
    ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.8, linewidth=2,
            label='Expected under null')
    
    # 95% CI
    try:
        alpha = 0.05
        lower_bound = -np.log10(np.maximum(1e-10, 
            stats.beta.ppf(alpha/2, np.arange(1, n+1), np.arange(n, 0, -1))))
        upper_bound = -np.log10(np.maximum(1e-10, 
            stats.beta.ppf(1-alpha/2, np.arange(1, n+1), np.arange(n, 0, -1))))
        
        ax.fill_between(-np.log10(expected_p), lower_bound, upper_bound, 
                        alpha=0.2, color='gray', label='95% CI')
    except:
        pass
    
    ax.set_xlabel('Expected -log₁₀(p-value)', fontsize=16)
    ax.set_ylabel('Observed -log₁₀(p-value)', fontsize=16)
    ax.set_title(f'QQ Plot (No total_embryos covariate)\n'
                 f'MAF ≥ 0.05, INFO ≥ 0.20\n'
                 f'n = {n:,} variants, λ = {lambda_val:.4f}',
                 fontsize=18, fontweight='bold')
    
    ax.text(0.05, 0.95, f'λ = {lambda_val:.4f}\nn = {n:,}', 
            transform=ax.transAxes, fontsize=14,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.grid(True, alpha=0.3)
    ax.legend(loc='lower right', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")

def create_manhattan_plot(gwas_df, lambda_val, output_file):
    """Create genome-wide Manhattan plot"""
    print("\nCreating Manhattan plot (all chromosomes)...")
    
    df = gwas_df.copy()
    
    # Calculate cumulative positions
    chr_lengths = df.groupby('chromosome')['position'].max()
    chr_starts = {}
    cumulative_pos = 0
    
    for chr_num in range(1, 23):
        chr_starts[chr_num] = cumulative_pos
        if chr_num in chr_lengths.index:
            cumulative_pos += chr_lengths[chr_num] + 5000000
    
    df['cum_pos'] = df.apply(lambda row: chr_starts.get(row['chromosome'], 0) + row['position'], axis=1)
    
    # Calculate chromosome centers
    chr_centers = {}
    for chr_num in range(1, 23):
        chr_data = df[df['chromosome'] == chr_num]
        if len(chr_data) > 0:
            chr_centers[chr_num] = chr_data['cum_pos'].mean()
    
    fig, ax = plt.subplots(figsize=(20, 10))
    
    colors = ['#1b9e77', '#d95f02']
    
    for chr_num in range(1, 23):
        chr_data = df[df['chromosome'] == chr_num]
        if len(chr_data) > 0:
            color = colors[chr_num % 2]
            
            # Non-significant
            non_sig = chr_data[chr_data['p'] >= 1e-5]
            if len(non_sig) > 0:
                ax.scatter(non_sig['cum_pos'], non_sig['neg_log10_p'], 
                          c=color, s=1, alpha=0.5)
            
            # Suggestive
            suggestive = chr_data[(chr_data['p'] < 1e-5) & (chr_data['p'] >= 5e-8)]
            if len(suggestive) > 0:
                ax.scatter(suggestive['cum_pos'], suggestive['neg_log10_p'], 
                          c='orange', s=20, alpha=0.8, edgecolors='black', linewidth=0.5)
            
            # Genome-wide significant
            significant = chr_data[chr_data['p'] < 5e-8]
            if len(significant) > 0:
                ax.scatter(significant['cum_pos'], significant['neg_log10_p'], 
                          c='red', s=40, alpha=1.0, edgecolors='black', linewidth=1,
                          marker='D', label='p < 5e-8' if chr_num == 1 else "")
    
    # Significance lines
    ax.axhline(y=-np.log10(5e-8), color='red', linestyle='--', linewidth=1.5, 
               label='Genome-wide (5e-8)', alpha=0.7)
    ax.axhline(y=-np.log10(1e-5), color='orange', linestyle='--', linewidth=1.0, 
               label='Suggestive (1e-5)', alpha=0.7)
    
    ax.set_xlabel('Chromosome', fontsize=16, fontweight='bold')
    ax.set_ylabel('-log₁₀(p-value)', fontsize=16, fontweight='bold')
    ax.set_title(f'Manhattan Plot (No total_embryos covariate)\n'
                 f'MAF ≥ 0.05, INFO ≥ 0.20 | n = {len(df):,} variants | λ = {lambda_val:.4f}',
                 fontsize=18, fontweight='bold', pad=20)
    
    if chr_centers:
        ax.set_xticks(list(chr_centers.values()))
        ax.set_xticklabels(list(chr_centers.keys()), fontsize=12)
    
    max_y = max(12, df['neg_log10_p'].max() + 0.5)
    ax.set_ylim(0, max_y)
    
    ax.legend(loc='upper right', fontsize=11)
    ax.grid(True, axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")

def create_per_chromosome_plots(gwas_df, output_dir):
    """Create individual Manhattan plot for each chromosome"""
    print("\nCreating per-chromosome Manhattan plots...")
    
    chr_dir = os.path.join(output_dir, 'per_chromosome')
    os.makedirs(chr_dir, exist_ok=True)
    
    for chr_num in range(1, 23):
        chr_data = gwas_df[gwas_df['chromosome'] == chr_num].copy()
        
        if len(chr_data) == 0:
            continue
        
        fig, ax = plt.subplots(figsize=(14, 6))
        
        # All variants
        ax.scatter(chr_data['position']/1e6, chr_data['neg_log10_p'], 
                  s=10, alpha=0.6, c='darkblue', label='All variants')
        
        # Significant
        sig_data = chr_data[chr_data['p'] < 5e-8]
        if len(sig_data) > 0:
            ax.scatter(sig_data['position']/1e6, sig_data['neg_log10_p'], 
                      s=100, c='red', marker='D', edgecolors='black', linewidth=1,
                      label=f'{len(sig_data)} genome-wide sig.')
        
        # Suggestive
        sug_data = chr_data[(chr_data['p'] >= 5e-8) & (chr_data['p'] < 1e-5)]
        if len(sug_data) > 0:
            ax.scatter(sug_data['position']/1e6, sug_data['neg_log10_p'], 
                      s=50, c='orange', marker='o', edgecolors='black', linewidth=0.5,
                      alpha=0.8, label=f'{len(sug_data)} suggestive')
        
        # Significance lines
        ax.axhline(y=-np.log10(5e-8), color='red', linestyle='--', alpha=0.5, label='p = 5e-8')
        ax.axhline(y=-np.log10(1e-5), color='orange', linestyle='--', alpha=0.3, label='p = 1e-5')
        
        ax.set_xlabel(f'Position on Chromosome {chr_num} (Mb)', fontsize=13, fontweight='bold')
        ax.set_ylabel('-log₁₀(p-value)', fontsize=13, fontweight='bold')
        ax.set_title(f'Chromosome {chr_num} (No total_embryos)\n{len(chr_data):,} variants',
                    fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=10)
        
        plt.tight_layout()
        chr_file = os.path.join(chr_dir, f'manhattan_chr{chr_num}.png')
        plt.savefig(chr_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  Chr {chr_num}: {len(chr_data):,} variants, {len(sig_data)} sig., {len(sug_data)} suggestive")
    
    print(f"  Saved to: {chr_dir}/")

def create_effect_size_plots(gwas_df, lambda_val, output_file):
    """Create effect size distribution plots"""
    print("\nCreating effect size plots...")
    
    sig_results = gwas_df[gwas_df['p'] < 1e-5].copy()
    
    if len(sig_results) == 0:
        print("  No significant results for effect size plot")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Plot 1: Effect size vs p-value
    ax1.scatter(sig_results['beta'], -np.log10(sig_results['p']), 
               alpha=0.6, s=30, c='darkblue')
    ax1.axhline(y=-np.log10(5e-8), color='red', linestyle='--', alpha=0.5, label='p = 5e-8')
    ax1.axvline(x=0, color='gray', linestyle=':', alpha=0.5)
    ax1.set_xlabel('Effect Size (β)', fontsize=13)
    ax1.set_ylabel('-log₁₀(p-value)', fontsize=13)
    ax1.set_title('Effect Size vs Significance', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Effect size distribution
    ax2.hist(sig_results['beta'], bins=30, alpha=0.7, color='darkblue', edgecolor='black')
    ax2.axvline(x=0, color='red', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Effect Size (β)', fontsize=13)
    ax2.set_ylabel('Count', fontsize=13)
    ax2.set_title(f'Effect Size Distribution (p < 1e-5)\nn = {len(sig_results)} variants', 
                 fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    mean_beta = sig_results['beta'].mean()
    median_beta = sig_results['beta'].median()
    ax2.text(0.02, 0.98, f'Mean β: {mean_beta:.3f}\nMedian β: {median_beta:.3f}',
             transform=ax2.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    fig.suptitle(f'Effect Sizes (No total_embryos | λ = {lambda_val:.4f})',
                fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")

def create_top_hits_table(gwas_df, output_file):
    """Create table of top hits"""
    print("\nCreating top hits table...")
    
    top_hits = gwas_df.nsmallest(100, 'p')[
        ['variant', 'chromosome', 'position', 'beta', 'se', 't', 'p', 'n']
    ].copy()
    
    top_hits.to_csv(output_file, index=False)
    print(f"  Saved: {output_file}")
    
    return top_hits

# ============================================================================
# MAIN
# ============================================================================

def main():
    print("="*80)
    print("GWAS PLOTS - NO TOTAL_EMBRYOS COVARIATE")
    print("="*80)
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load GWAS
    gwas_df = load_gwas_with_positions(GWAS_FILE)
    
    if gwas_df is None or len(gwas_df) == 0:
        print("ERROR: No data to plot!")
        return
    
    # Calculate lambda
    lambda_val = calculate_lambda(gwas_df)
    print(f"\nLambda (λ) = {lambda_val:.4f}")
    
    # Create plots
    print("\nCreating plots...")
    print("-" * 80)
    
    create_qq_plot(gwas_df, lambda_val, 
                   os.path.join(OUTPUT_DIR, 'qq_plot.png'))
    
    create_manhattan_plot(gwas_df, lambda_val,
                         os.path.join(OUTPUT_DIR, 'manhattan_all_chromosomes.png'))
    
    create_per_chromosome_plots(gwas_df, OUTPUT_DIR)
    
    create_effect_size_plots(gwas_df, lambda_val,
                            os.path.join(OUTPUT_DIR, 'effect_size_plots.png'))
    
    # Create tables
    top_hits = create_top_hits_table(gwas_df,
                                     os.path.join(OUTPUT_DIR, 'top_100_hits.csv'))
    
    # Summary
    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)
    
    sig_5e8 = len(gwas_df[gwas_df['p'] < 5e-8])
    sig_1e5 = len(gwas_df[(gwas_df['p'] >= 5e-8) & (gwas_df['p'] < 1e-5)])
    
    print(f"\nSummary:")
    print(f"  Model: NO total_embryos covariate")
    print(f"  Covariates: mat_age + PC1-8 + seq_date + geno")
    print(f"  Total variants: {len(gwas_df):,}")
    print(f"  Lambda: {lambda_val:.4f}")
    print(f"  Genome-wide significant (p < 5e-8): {sig_5e8}")
    print(f"  Suggestive (5e-8 ≤ p < 1e-5): {sig_1e5}")
    
    if len(top_hits) > 0:
        print(f"\nTop 5 hits:")
        print(top_hits.head()[['variant', 'chromosome', 'position', 'p', 'beta']].to_string(index=False))
    
    print(f"\nAll outputs saved to: {OUTPUT_DIR}/")

if __name__ == "__main__":
    main()