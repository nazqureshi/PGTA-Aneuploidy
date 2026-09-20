#!/bin/bash
#SBATCH --job-name=gls_compute_all_chr  # Job name
#SBATCH --time=240:00:00            
#SBATCH --partition=genetics_1
#SBATCH --mem=150G
#SBATCH --cpus-per-task=24             # Number of CPU cores per task
#SBATCH --export=ALL
#SBATCH --output=/path/to/project/gls/logs/slurm.%N.%j.out  
#SBATCH --error=/path/to/project/gls/logs/slurm.%N.%j.err   
#SBATCH --array=1-22%10                # Submit a job for each chromosome MAX 10 at a time
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu



set -euo pipefail                      # fail on errors 

# Reference path
REFGEN_PATH="/path/to/reference_data/hg38_ref"

# Sites path
SITES_PATH="/path/to/reference_data/glimpse/calgl"

# BAM list
BAM_LIST="/path/to/project/gls/scripts/bam_file_list.sorted.txt"

# Output path
OUT_PATH="/path/to/project/gls/outputs_calculated_gls"

LOG_DIR="/path/to/project/gls/logs"

# modules                        
module use /projects/community/redhat/modulefiles
module load bcftools/1.23.1
module load samtools

# ulimit for file descriptors -> avoids overloading the system with too many open files at once
ulimit -n 9999

# get the chromosome number from SLURM_ARRAY_TASK_ID
CHROM="${SLURM_ARRAY_TASK_ID}"          # CHANGED: quote variable assignment

# Reference genome for the current chromosome
REFGEN="${REFGEN_PATH}/chr${CHROM}.fa.gz"

# Sites VCF file for the current chromosome
VCF="${SITES_PATH}/1000GP.chr${CHROM}.cleaned.sites.vcf.gz"

# Sites TSV file for the current chromosome
TSV="${SITES_PATH}/1000GP.chr${CHROM}.cleaned.sites.tsv.gz"

# Output VCF file for the current chromosome
OUT="${OUT_PATH}/5616pts.chr${CHROM}.vcf.gz"

# Log output for the current chromosome
LOG="${LOG_DIR}/slurm.${SLURM_JOB_NODELIST}.${SLURM_JOB_ID}.chr${CHROM}.out"
ERR="${LOG_DIR}/slurm.${SLURM_JOB_NODELIST}.${SLURM_JOB_ID}.chr${CHROM}.err"


# troubleshooting/debugging errors with regards to running on cluster -> getting more verbose error logs 
exec > >(tee -a "$LOG") 2> >(tee -a "$ERR" >&2)  # CHANGED: send stdout/stderr to logs while preserving Slurm logs

echo "Starting chromosome ${CHROM} at $(date)"   # CHANGED
echo "Running on node(s): ${SLURM_JOB_NODELIST}" # CHANGED
echo "Using bcftools: $(which bcftools)"         # CHANGED
bcftools --version                              # CHANGED
echo "Using samtools: $(which samtools)"         # CHANGED
samtools --version | head -n 1                  # CHANGED

# Check if required files exist before running commands
if [ ! -f "$REFGEN" ]; then echo "Error: Missing reference genome $REFGEN"; exit 1; fi       # CHANGED
if [ ! -f "$VCF" ]; then echo "Error: Missing sites VCF $VCF"; exit 1; fi                    # CHANGED
if [ ! -f "$TSV" ]; then echo "Error: Missing sites TSV $TSV"; exit 1; fi                    # CHANGED
if [ ! -f "$BAM_LIST" ]; then echo "Error: Missing BAM list $BAM_LIST"; exit 1; fi           # CHANGED
if [ ! -x "./GLIMPSE_chunk" ]; then echo "Error: ./GLIMPSE_chunk missing or not executable"; exit 1; fi  # CHANGED

# Run Genotype Likelihood calculation for the current chromosome using bcftools
bcftools mpileup --threads 24 \
    -f ${REFGEN} \
    -I -E -a 'FORMAT/DP' \
    -T ${VCF} -r chr${CHROM} \
    -b ${BAM_LIST} -Ou | bcftools call --threads 24 -Aim -C alleles \
    -T ${TSV} -Oz -o ${OUT} >> "$LOG" 2>> "$ERR"

# Index the output VCF file
bcftools index --threads 24 -f ${OUT} >> "$LOG" 2>> "$ERR"

# Log the completion of the bcftools part
echo "Completed GL calculation for chromosome ${CHROM} at $(date),NEXT moving onto GLIMPSE CHUNKING" >> ${LOG}


# Check if GLIMPSE input file exists
GLIMPSE_INPUT="${SITES_PATH}/1000GP.chr${CHROM}.cleaned.sites.vcf.gz"
if [ ! -f "$GLIMPSE_INPUT" ]; then
    echo "Error: Missing GLIMPSE input $GLIMPSE_INPUT" | tee -a "$ERR"
    exit 1
fi



# GLIMPSE chunking step for the current chromosome
./GLIMPSE_chunk --threads 24 \
    --input ${SITES_PATH}/1000GP.chr${CHROM}.cleaned.sites.vcf.gz \
    --region chr${CHROM} \
    --window-mb 2 \
    --buffer-mb 0.2 \
    --output ${OUT_PATH}/chunks.chr${CHROM}.txt >> "$LOG" 2>> "$ERR"

# Log completion of the GLIMPSE chunking
echo "Completed GLIMPSE chunking for chromosome ${CHROM} at $(date)" >> ${LOG}
