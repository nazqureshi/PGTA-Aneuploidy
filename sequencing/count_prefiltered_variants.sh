#!/bin/bash

# Directories
base_dir="/path/to/project/imputation/outputs/ligate"
logs_dir="${base_dir}/logs"

for chr in {1..22}
do
  job_script="${base_dir}/impute_qual_chr${chr}.sh"
  
  cat > "$job_script" << EOF
#!/bin/bash
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --job-name=impute_qual_chr${chr}
#SBATCH --partition=main-redhat
#SBATCH --mem=48G
#SBATCH --cpus-per-task=1
#SBATCH --export=ALL
#SBATCH --output=${logs_dir}/impute_qual_chr${chr}.%N.%j.out
#SBATCH --error=${logs_dir}/impute_qual_chr${chr}.%N.%j.err

# Load modules
module use /projects/community/redhat/modulefiles
module load bcftools/1.23.1

# Extract ID and INFO for chromosome ${chr}
bcftools query -f '%ID\t%INFO\n' 5616pts_ligated_merged_chr${chr}.bcf > ${base_dir}/5616pts_PREFILTERED_impute_qual_chr${chr}.txt


# Count variant records (excluding headers)
num=\$(bcftools index -n 5616pts_ligated_merged_chr${chr}.bcf)
echo "chr${chr} has \${num} variants" >> ${base_dir}/variant_counts_PRE_filter_chr_ALL.txt

EOF

  sbatch "$job_script"
done
