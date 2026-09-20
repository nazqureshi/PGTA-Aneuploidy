#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=150GB
#SBATCH --time=72:00:00
#SBATCH --partition=genetics_1
#SBATCH --job-name=PCA
#SBATCH --output=../logs/pCA.%j.out
#SBATCH --error=../logs/pCA.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu


# Load required modules, in this case laser is in the scripts directory, so need to put the path to it properly

#So, Define LASER directory
# Cluster path root -- update to match your own environment.
PROJECT_ROOT="/path/to/project"
LASER_DIR="${PROJECT_ROOT}/ancestry_inference/scripts/LASER-2.04"

${LASER_DIR}/laser \
  -g ../data/reference/1kg.geno \
  -o ../data/reference/laser_PCA \
  -k 8 \
  -pca 1 #explicitly run in PCA mode by enabling it, this is written at the end
