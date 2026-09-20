#!/bin/bash
#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=merge_bam
#SBATCH --time=10:00:00
#SBATCH --export=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24GB

# ---------------------------------------------------------------------------
# merge_bam_array.sh
#
# SLURM job array script for merging BAM files by patient ID.
#
# Expects environment variables (passed via --export at submission):
#   BAM_LIST    — path to the job array input file (one patient per line)
#   OUTPUT_DIR  — directory for merged BAM output
#   LOG_DIR     — directory for log files
# ---------------------------------------------------------------------------

module load samtools
module load python
module load bcftools
source ~/pandas_env/bin/activate

# Read the line corresponding to this array task
PATIENT_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BAM_LIST}")
PATIENT_ID=$(echo "${PATIENT_LINE}" | cut -d',' -f1)
BAM_FILES=$(echo "${PATIENT_LINE}" | cut -d',' -f2-)

echo "Task ${SLURM_ARRAY_TASK_ID}: Patient ${PATIENT_ID}"
echo "BAM files: ${BAM_FILES}"

# --- Validation ---

if [ -z "${PATIENT_LINE}" ]; then
    echo "Error: No line found for task ${SLURM_ARRAY_TASK_ID} in ${BAM_LIST}" >&2
    exit 1
fi

if [ -z "${BAM_FILES}" ]; then
    echo "Error: No BAM files for patient ${PATIENT_ID}" >&2
    exit 1
fi

for file in ${BAM_FILES}; do
    if [ ! -f "${file}" ]; then
        echo "Error: BAM file does not exist: ${file}" >&2
        exit 1
    fi
done

# --- Merge ---

OUTPUT_BAM="${OUTPUT_DIR}/${PATIENT_ID}_merged.bam"

# Skip if output already exists (safe re-runs)
if [ -f "${OUTPUT_BAM}" ]; then
    echo "Output already exists, skipping: ${OUTPUT_BAM}"
    exit 0
fi

echo "Merging BAM files for patient ${PATIENT_ID} -> ${OUTPUT_BAM}"
samtools merge "${OUTPUT_BAM}" ${BAM_FILES}
MERGE_EXIT=$?

if [ ${MERGE_EXIT} -ne 0 ]; then
    echo "Error: samtools merge failed for patient ${PATIENT_ID} (exit code ${MERGE_EXIT})" >&2
    rm -f "${OUTPUT_BAM}"
    exit ${MERGE_EXIT}
fi

echo "Done merging for patient ${PATIENT_ID}"
