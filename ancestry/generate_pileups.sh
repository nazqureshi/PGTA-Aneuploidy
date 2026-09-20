#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=150GB
#SBATCH --time=96:00:00
#SBATCH --partition=genetics_1
#SBATCH --job-name=generate_pileups
#SBATCH --output=../logs/generate_pileups.%j.out
#SBATCH --error=../logs/generate_pileups.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu

# Load required modules
module use /projects/community/redhat/modulefiles
module load bcftools/1.23.1
module load samtools

# Cluster path root -- update to match your own environment.
REF_DATA_ROOT="/path/to/reference_data"

# Step E: Convert BAM to pileup
N=8
pids=""
COUNTER=0

while IFS= read -r bam_file || [ -n "${bam_file}" ]; do
  ((COUNTER=COUNTER%N)); ((COUNTER++==0)) && wait $pids && pids=""
  SAMPID=$(basename "${bam_file::-4}")
  samtools mpileup -q 30 -Q 20 -f "${REF_DATA_ROOT}/hg38_ref/hg38.fa" -l ../data/reference/1kg.bed "${bam_file}" > ../data/pileup/${SAMPID}.pileup &
  pids="$pids $!"
done < ../data/bam/bam_list.txt
wait $pids

# List all pileup files
ls -A1 ../data/pileup/*.pileup > ../data/pileup/pileup_list.txt
