#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=150GB
#SBATCH --time=96:00:00
#SBATCH --partition=genetics_1
#SBATCH --job-name=prepare_reference
#SBATCH --output=../logs/prepare_reference.%j.out
#SBATCH --error=../logs/prepare_reference.%j.err

# Load required modules
module use /projects/community/redhat/modulefiles
module load bcftools/1.23.1
module load samtools

# Set number of parallel jobs
N=8
declare -a pids  # Use an array to store PIDs

# Step A: Extract biallelic SNPs with MAF > 10% from 1000G VCFs

for i in {1..22}; do
  input_vcf="../data/reference/ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
  output_vcf="../data/reference/chr${i}_1kg_maf10.vcf.gz"

  if [ ! -f "$input_vcf" ]; then
    echo "WARNING: Missing input file $input_vcf" >&2
    continue
  fi

  ((i=i%N)); ((i++==0)) && wait "${pids[@]}" && unset pids
  bcftools view \
    "$input_vcf" \
    --exclude-types indels,mnps,ref,bnd,other \
    --min-alleles 2 \
    --max-alleles 2 \
    --min-af 0.1:minor \
    --phased \
    --exclude 'AN!=2*N_SAMPLES' \
    --output-file "$output_vcf" \
    --output-type z &

  pids+=($!)
done
wait "${pids[@]}"

# Step B: Merge filtered VCFs
rm -f ../data/reference/input_vcf_list.txt
for i in {1..22}; do
  if [ -f "../data/reference/chr${i}_1kg_maf10.vcf.gz" ]; then
    echo "../data/reference/chr${i}_1kg_maf10.vcf.gz" >> ../data/reference/input_vcf_list.txt
  else
    echo "WARNING: Missing chr${i}_1kg_maf10.vcf.gz" >&2
  fi
done

bcftools concat --threads $N --file-list ../data/reference/input_vcf_list.txt -o ../data/reference/1kg_maf10.vcf.gz -O z

rm -f ../data/reference/chr*.vcf.gz