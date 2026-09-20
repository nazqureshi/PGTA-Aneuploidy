#!/bin/bash
#SBATCH --job-name=fix_chr21_overlap
#SBATCH --partition=main-redhat
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/path/to/project/imputation/outputs/imputed/logs/slurm.fix_chr21_overlap.%N.%j.out
#SBATCH --error=/path/to/project/imputation/outputs/imputed/logs/slurm.fix_chr21_overlap.%N.%j.err

# Load modules
module use /projects/community/redhat/modulefiles
module load bcftools/1.23.1

PROJECT_ROOT="/path/to/project"

cd "${PROJECT_ROOT}/imputation/outputs/imputed" || exit 1

for s in sample0 sample1; do
    echo "Processing ${s}.chr21.bcf"

    cp ${s}.chr21.bcf ${s}.chr21.ORIGINAL_OVERLAP.bcf
    cp ${s}.chr21.bcf.csi ${s}.chr21.ORIGINAL_OVERLAP.bcf.csi

    bcftools view \
      --exclude 'CHROM="chr21" && POS=10804701' \
      -O b \
      -o ${s}.chr21.tmp.bcf \
      ${s}.chr21.bcf

    bcftools index -f ${s}.chr21.tmp.bcf

    mv ${s}.chr21.tmp.bcf ${s}.chr21.bcf
    mv ${s}.chr21.tmp.bcf.csi ${s}.chr21.bcf.csi

    echo "Finished ${s}.chr21.bcf"
done
