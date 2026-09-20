#!/bin/bash
#SBATCH --partition=main-redhat
#SBATCH --job-name=index_chr
#SBATCH --array=1-22%10
#SBATCH --output=logs/index_%A_%a.out
#SBATCH --error=logs/index_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --requeue
#SBATCH --cpus-per-task=1
#SBATCH --export=ALL

chr=$SLURM_ARRAY_TASK_ID

# Cluster path root -- update to match your own environment.
WORK_DIR="/path/to/project/imputation/outputs/ligate/maffilt"

module use /projects/community/redhat/modulefiles
module load bcftools/1.23.1

bcftools index "${WORK_DIR}/5616pts_ligated_merged_chr${chr}_maffilt.bcf"
