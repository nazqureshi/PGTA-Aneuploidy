#!/bin/bash

# Script to generate 22 individual SLURM scripts for chromosome-specific filtering
# Usage: bash generate_af_filt_scripts.sh

# Cluster path root -- update to match your own environment.
PROJECT_ROOT="/path/to/project"

# Create directory for individual scripts if it doesn't exist
mkdir -p chr_scripts

# Generate individual SLURM scripts for each chromosome
for chr in {1..22}; do
    script_name="chr_scripts/af_filt_chr${chr}.slurm"
    
    cat > "$script_name" << EOF
#!/bin/bash
#SBATCH --partition=main-redhat
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --export=ALL
#SBATCH --mem=96GB
#SBATCH --job-name=step_0_af_filt_chr${chr}
#SBATCH --output=${PROJECT_ROOT}/imputation/outputs/ligate/maffilt/logs/af_filt_chr${chr}.%j.out
#SBATCH --error=${PROJECT_ROOT}/imputation/outputs/ligate/maffilt/logs/af_filt_chr${chr}.%j.err
#SBATCH --requeue

# Load modules
module use /projects/community/redhat/modulefiles
module load bcftools/1.23.1

cd ${PROJECT_ROOT}/imputation/outputs/ligate

echo "Starting filtering for chromosome ${chr} at \$(date)"

# Filter variants by imputation score >= 0.2 AND minor allele frequency >= 0.0.0 Since now testing different filters to reduce p-val inflation
bcftools view -i 'INFO>=0.2' --threads 32 \\
    5616pts_ligated_merged_chr${chr}.bcf -Ou | \\
bcftools view -q 0.05:minor --threads 32 \\
    -Ou -o maffilt/5616pts_ligated_merged_chr${chr}_maffilt.bcf

echo "Filtering completed for chromosome ${chr}, extracting genotype dosages..."

# Extract genotype dosage
bcftools query -f '%ID [\\t%DS]\\n' \\
    maffilt/5616pts_ligated_merged_chr${chr}_maffilt.bcf \\
    > maffilt/gt_5616pts_ligated_merged_chr${chr}_maffilt.txt

echo "Chromosome ${chr} processing completed at \$(date)"
EOF

    echo "Generated script: $script_name"
done

echo ""
echo "All 22 scripts generated in chr_scripts/ directory"
echo ""
echo "To submit all jobs, run:"
echo "for script in chr_scripts/af_filt_chr*.slurm; do sbatch \$script; done"
echo ""
echo "To check job status:"
echo "squeue -u \$USER"
