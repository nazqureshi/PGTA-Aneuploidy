#!/bin/bash
#SBATCH --partition=main-redhat
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --export=ALL
#SBATCH --mem=96GB
#SBATCH --job-name=modify_phenotype_file
#SBATCH --output=/path/to/project/imputation/outputs/ligate/maffilt/logs/step1_edit_pheno.%j.out
#SBATCH --error=/path/to/project/imputation/outputs/ligate/maffilt/logs/step1_edit_pheno.%j.err

# Cluster path root -- update to match your own environment.
WORK_DIR="/path/to/project/imputation/outputs/ligate/maffilt"

# Load Python (adjust module name as needed for Amarel)
module load python  # or whatever Python module is available

# activate python environment to use pandas
source ~/pandas_env/bin/activate

# Change to working directory
cd "${WORK_DIR}"

# Run the merge script
python add_seqdate_column.py


echo "Step 1 phenotype file modification completed"
