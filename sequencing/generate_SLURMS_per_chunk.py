import os

# Cluster path roots -- update these to match your own environment.
PROJECT_ROOT = "/path/to/project"
REF_DATA_ROOT = "/path/to/reference_data"


def generate_slurm_scripts(chunk_file_path, chromosome_number, output_dir):
    # Make sure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Read the chunk lines
    with open(chunk_file_path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    for i, line in enumerate(lines):
        chunk_id = line.split()[0]
        script_filename = os.path.join(output_dir, f"impute_chr{chromosome_number}_{chunk_id}.slurm")

        slurm_script = f"""#!/bin/bash
#SBATCH --partition=main-redhat
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --requeue
#SBATCH --ntasks=1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --output={PROJECT_ROOT}/imputation/logs/chr{chromosome_number}_impute/slurm.impute.%N.%j.out
#SBATCH --error={PROJECT_ROOT}/imputation/logs/chr{chromosome_number}_impute/slurm.impute.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@institution.edu

# Load modules
module use /projects/community/redhat/modulefiles
module load bcftools/1.23.1


# Set chromosome ID
CHR_ID={chromosome_number}

# Set input files
VCF={PROJECT_ROOT}/gls/outputs_calculated_gls/5616pts.chr${{CHR_ID}}.vcf.gz
REF={REF_DATA_ROOT}/glimpse/calgl/1000GP.chr${{CHR_ID}}.cleaned.bcf
MAP={REF_DATA_ROOT}/genetic_maps.b38/chr${{CHR_ID}}.b38.gmap.gz

# ----------- MANUALLY SET THE LINE CONTENTS HERE -------------
LINE="{line}"
# -------------------------------------------------------------

# Parse the line
CHUNK_ID=$(echo $LINE | cut -d" " -f1)
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)

OUT={PROJECT_ROOT}/imputation/outputs/imputed/sample${{CHUNK_ID}}.chr${{CHR_ID}}.bcf

# Run GLIMPSE2_phase
{PROJECT_ROOT}/imputation/scripts/GLIMPSE2_phase --input-gl ${{VCF}} \\
                --reference ${{REF}} \\
                --map ${{MAP}} \\
                --input-region ${{IRG}} \\
                --output-region ${{ORG}} \\
                --output ${{OUT}} \\
                --threads 16

# Index the output
bcftools index -f ${{OUT}}
"""

        with open(script_filename, "w") as out_file:
            out_file.write(slurm_script)

        print(f"Generated: {script_filename}")

# Usage: one call per autosome (chr1-chr22), pointing at that chromosome's
# GLIMPSE chunk file produced by gl_calculation.sh.
for chrom in range(1, 23):
    chunk_file = f"{PROJECT_ROOT}/imputation/outputs/chunks/chunks.chr{chrom}.txt"
    generate_slurm_scripts(chunk_file, chrom, ".")
