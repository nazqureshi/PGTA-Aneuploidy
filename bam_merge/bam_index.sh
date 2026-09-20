#!/bin/bash
#SBATCH --partition=main-redhat
#SBATCH --job-name=bam_index
#SBATCH --output=/path/to/project/bam_merge_sort/logs/bam_index.%A_%a.out
#SBATCH --error=/path/to/project/bam_merge_sort/logs/bam_index.%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --export=ALL
#SBATCH --array=1-281%50

set -euo pipefail

#NOTE: This job is being run on main, and as of May 2026 these nodes are not on redhat9 yet, so cannot load the modules via the redhat modules, need to run from amarel 1 or 2

#source /etc/profile.d/modules.sh
#module use /projects/community/redhat/modulefiles
module load samtools

# ==========================================================
# Paths
# ==========================================================
# Cluster path root -- update to match your own environment.
PROJECT_ROOT="/path/to/project"

OUTPUT_DIR="${PROJECT_ROOT}/bam_merge_sort/outputs"
FILE_LIST="${OUTPUT_DIR}/bam_files.txt"

# Each array task indexes 20 BAM files
FILES_PER_TASK=20

START_LINE=$(( (SLURM_ARRAY_TASK_ID - 1) * FILES_PER_TASK + 1 ))
END_LINE=$(( SLURM_ARRAY_TASK_ID * FILES_PER_TASK ))

echo "======================================="
echo "Started at: $(date)"
echo "Node: $(hostname)"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing lines ${START_LINE}-${END_LINE} from ${FILE_LIST}"
echo "======================================="

sed -n "${START_LINE},${END_LINE}p" "${FILE_LIST}" | while read -r BAM_FILE; do
    [ -z "${BAM_FILE}" ] && continue

    BAM_BASENAME=$(basename "${BAM_FILE}")
    INDEX_OUT="${OUTPUT_DIR}/${BAM_BASENAME}.bai"

    echo "---------------------------------------"
    echo "Input BAM: ${BAM_FILE}"
    echo "Output BAI: ${INDEX_OUT}"

    samtools index \
        -@ "${SLURM_CPUS_PER_TASK}" \
        "${BAM_FILE}" \
        "${INDEX_OUT}"
done

echo "======================================="
echo "Finished at: $(date)"
echo "======================================="
