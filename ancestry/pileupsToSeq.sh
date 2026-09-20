#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100GB
#SBATCH --time=72:00:00
#SBATCH --partition=genetics_1
#SBATCH --job-name=pileupToSeq
#SBATCH --output=../logs/pileupToSeq.%j.out
#SBATCH --error=../logs/pileupToSeq.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu


# Load required modules
module use /projects/community/redhat/modulefiles
module load bcftools/1.23.1
module load samtools
module load python/2.7.12

# Cluster path roots -- update these to match your own environment.
PROJECT_ROOT="/path/to/project"
REF_DATA_ROOT="/path/to/reference_data"

python -B "${PROJECT_ROOT}/ancestry_inference/scripts/LASER-2.04/pileup2seq/pileup2seq.py" \
  -f "${REF_DATA_ROOT}/hg38_ref/hg38.fa" \
  -m "${PROJECT_ROOT}/ancestry_inference/data/reference/1kg.site" \
  -o ../data/seq/target_samples \
  $(cat "${PROJECT_ROOT}/ancestry_inference/data/pileup/pileup_list.txt")
