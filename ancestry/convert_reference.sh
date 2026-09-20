#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=159GB
#SBATCH --time=48:00:00
#SBATCH --partition=genetics_1
#SBATCH --job-name=convert_ref
#SBATCH --output=../logs/convert_reference.%j.out
#SBATCH --error=../logs/convert_reference.%j.err

# Cluster path root -- update to match your own environment.
PROJECT_ROOT="/path/to/project"

# Define LASER directory
LASER_DIR="${PROJECT_ROOT}/ancestry_inference/scripts/LASER-2.04"

# Step C: Convert VCF to LASER format
if [ ! -f "../data/reference/1kg.geno" ]; then
  "$LASER_DIR/vcf2geno/vcf2geno" \
    --inVcf ../data/reference/1kg_maf10.vcf.gz \
    --out ../data/reference/1kg
fi

# Step D: Generate BED file
if [ ! -f "../data/reference/1kg.bed" ]; then
  awk '{if (NR > 1) {print "chr"$1, $2-1, $2;}}' ../data/reference/1kg.site > ../data/reference/1kg.bed
fi
