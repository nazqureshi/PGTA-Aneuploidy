#!/bin/bash
#SBATCH --job-name=glm_quasi
#SBATCH --partition=main-redhat
#SBATCH --array=6-21
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=72:00:00
#SBATCH --output=logs/R_glm_chr%a.out
#SBATCH --error=logs/R_glm_chr%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu
#SBATCH --requeue

# Create logs directory if it doesn't exist
mkdir -p logs

# Load R module (adjust version as needed for Amarel)
module use /projects/community/redhat/modulefiles
module load bcftools/1.23.1
module load R

export R_LIBS_USER="~/R/x86_64-pc-linux-gnu-library/4.1"

# Run the R script with chromosome number
Rscript glm_quasibinomial.R ${SLURM_ARRAY_TASK_ID}
