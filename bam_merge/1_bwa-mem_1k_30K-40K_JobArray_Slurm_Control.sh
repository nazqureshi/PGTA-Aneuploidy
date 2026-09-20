#!/bin/bash

#SBATCH --partition=genetics_1       # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=bwa-mem-1k        # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=1           # Cores per task (>1 if multithread tasks)
#SBATCH --mem=4000                  # Real memory (RAM) required (MB)
#SBATCH --time=240:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=/path/to/project/logs/bwa-mem.batch.%A-%a.out     # STDOUT output file
#SBATCH --error=/path/to/project/logs/bwa-mem.batch.%A-%a.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

# Cluster path roots -- update these to match your own environment.
PROJECT_ROOT="/path/to/project"
SCRATCH_LOG_DIR="/path/to/scratch/logs"

#====Working folder=====#

cd "${PROJECT_ROOT}/fastq_batch2"


# Script Name: 1_bwa-mem_1k_30K-40K_JobArray_Slurm_Control.sh
# Description: This script process 1000 fastq files in one slurm job array, loops through 40,000 files in increments of 1,000.
# The program starts from file starting_file, does not cover the final batch of 41 files.

# Author: J.X.
# Date: 2024-12-26



# Define the total number of files, the increment, the starting file, and j for log handling
total_files=10000
increment=1000
starting_file=30000
j=0

# Loop through the files
for ((i=0; i<total_files; i+=increment)); do
    
    # Define the number of lines in the squeue output
    squeue_output_line=0

    # Loop until the number of squeue output lines equals 2 (one header line, one job for the current script)
    # If there are more than two lines (i.e., jobs running, waiting 30 seconds and check again)
    while [ "$squeue_output_line" -ne 2 ]; do
        echo "squeue has job runing. Waiting for 30 seconds..."
        sleep 30
        squeue_output_line=`squeue -u $USER | wc -l`
    done


    ##### Handle log files, starting from the 2nd loop #####
    if (( "$j" != 0 )); then

        echo "Collecting log file: $j"

        # Combine Error File
        grep -Ev "idx|process" ${SCRATCH_LOG_DIR}/bwa-mem.batch.*.err > "${PROJECT_ROOT}/logs/bwa-mem.batch30K_${j}_nolog.err"
        # Combine Out File
        grep -Ev "batch_num" ${SCRATCH_LOG_DIR}/bwa-mem.batch.*.out > "${PROJECT_ROOT}/logs/bwa-mem.batch30K_${j}_nolog.out"

        mv ${SCRATCH_LOG_DIR}/bwa-mem.batch.*.err ${SCRATCH_LOG_DIR}/ind_logs
        mv ${SCRATCH_LOG_DIR}/bwa-mem.batch.*.out ${SCRATCH_LOG_DIR}/ind_logs
    fi

    # Modify j for log collection
    j=$((j + 1))



    ##### Submit the next batch #####

    echo "Processing files from $((i + starting_file + 1)) to $((i + starting_file + increment))"
    echo "sbatch --export=ALL,BATCH=$((i + starting_file)) ${PROJECT_ROOT}/scripts/1.bwa-mem_Batch1k_test.slurm"
    sbatch --export=ALL,BATCH=$((i + starting_file)) "${PROJECT_ROOT}/scripts/1.bwa-mem_Batch1k_test.slurm"

done



# Define the number of lines in the squeue output
squeue_output_line=0

# Loop until the number of squeue output lines equals 2 (one header line, one job for the current script)
# If there are more than two lines (i.e., jobs running, waiting 30 seconds and check again)
while [ "$squeue_output_line" -ne 2 ]; do
    echo "squeue has job runing. Waiting for 30 seconds..."
    sleep 30
    squeue_output_line=`squeue -u $USER | wc -l`
done

# Process final batch of log files

echo "Collecting final batch log file: $j"

# Combine Error File
grep -Ev "idx|process" ${SCRATCH_LOG_DIR}/bwa-mem.batch.*.err > "${PROJECT_ROOT}/logs/bwa-mem.batch30K_${j}_nolog.err"
# Combine Out File
grep -Ev "batch_num" ${SCRATCH_LOG_DIR}/bwa-mem.batch.*.out > "${PROJECT_ROOT}/logs/bwa-mem.batch30K_${j}_nolog.out"

mv ${SCRATCH_LOG_DIR}/bwa-mem.batch.*.err ${SCRATCH_LOG_DIR}/ind_logs
mv ${SCRATCH_LOG_DIR}/bwa-mem.batch.*.out ${SCRATCH_LOG_DIR}/ind_logs
