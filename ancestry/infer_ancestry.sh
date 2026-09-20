#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150GB
#SBATCH --time=112:00:00
#SBATCH --partition=genetics_1
#SBATCH --job-name=infer_ancestry
#SBATCH --export=ALL
#SBATCH --output=../logs/infer_ancestry.%j.out
#SBATCH --error=../logs/infer_ancestry.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu



# Load required modules, in this case laser is in the scripts directory, so need to put the path to it properly

#So, Define LASER directory
# Cluster path root -- update to match your own environment.
PROJECT_ROOT="/path/to/project"
LASER_DIR="${PROJECT_ROOT}/ancestry_inference/scripts/LASER-2.04"


# Step G: Infer ancestry
"$LASER_DIR/laser" \
  -c ../data/reference/laser_PCA.RefPC.coord\
  -g ../data/reference/1kg.geno \
  -s ../data/seq/target_samples.seq \
  -fmt 0 \
  -k 8 \
  -minc 0.01 \
  -seed 0 \
  -nt 24 \
  -r 5
