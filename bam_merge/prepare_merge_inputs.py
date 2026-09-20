#!/usr/bin/env python3
"""
prepare_merge_inputs.py

Reads the filtered Dataset ONE CSV, groups BAM files by patient ID, and writes
chunked job array input files for SLURM submission.

The dataset CSV contains a 'file_name' column with .fastq.gz filenames.
These are converted to .bam filenames and resolved to absolute paths under BAM_DIR.
Patients are grouped and split into chunks of <= CHUNK_SIZE to respect Amarel's
max SLURM array size of 1001.

Output files (written to SCRIPT_DIR):
    DS_ONE_bam_file_list_1.txt
    DS_ONE_bam_file_list_2.txt
    ...

Each line format:
    patient_ID,/full/path/to/file1.bam /full/path/to/file2.bam ...
"""

import os
import pandas as pd

# ===========================================================================
# Configuration — all file paths defined here
# ===========================================================================

# Cluster path root -- update to match your own environment.
# (PGTA_ROOT is the top-level project directory; 1_Final_Pipeline_4_2026 and
# bam/ are sibling subdirectories under it.)
PGTA_ROOT   = "/path/to/project_root"
PROJECT_ROOT = f"{PGTA_ROOT}/1_Final_Pipeline_4_2026"

DATASET_CSV = f"{PROJECT_ROOT}/dataset_ONE.csv"
BAM_DIR     = f"{PGTA_ROOT}/bam"
SCRIPT_DIR  = f"{PROJECT_ROOT}/bam_merge_sort/scripts"

OUTPUT_PREFIX = "DS_ONE"
CHUNK_SIZE    = 500   # Max patients per SLURM array (Amarel limit is 1001)

# ===========================================================================
# Functions
# ===========================================================================

def prepare_bam_list(csv_path, bam_dir):
    """Read dataset CSV, convert filenames to BAM paths, group by patient."""
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} rows from {csv_path}")

    # The CSV has .fastq.gz filenames — convert to .bam and build full paths
    df["file_name"] = df["file_name"].str.replace(".fastq.gz", ".bam", regex=False)
    df["file_path"] = df["file_name"].apply(lambda x: os.path.join(bam_dir, x))

    grouped = df.groupby("patient_ID")["file_path"].apply(list)
    return grouped


def split_into_chunks(grouped, chunk_size):
    """Split grouped patients into chunks without breaking patient groups apart."""
    patient_ids = list(grouped.index)
    chunks = []
    for i in range(0, len(patient_ids), chunk_size):
        chunk_ids = patient_ids[i : i + chunk_size]
        chunks.append(grouped.loc[chunk_ids])
    return chunks


def write_chunk_file(chunk, output_path):
    """Write one chunk to a job array input file (one line per patient)."""
    with open(output_path, "w") as f:
        for patient_id, file_paths in chunk.items():
            cleaned = [fp.strip() for fp in file_paths if fp.strip()]
            f.write(f"{patient_id},{' '.join(cleaned)}\n")


# ===========================================================================
# Main
# ===========================================================================

def main():
    grouped = prepare_bam_list(DATASET_CSV, BAM_DIR)
    print(f"Total patients: {len(grouped)}")
    print(f"Total BAM files: {sum(len(v) for v in grouped)}")

    chunks = split_into_chunks(grouped, CHUNK_SIZE)
    print(f"Chunks needed: {len(chunks)}")

    os.makedirs(SCRIPT_DIR, exist_ok=True)
    for i, chunk in enumerate(chunks, start=1):
        output_path = os.path.join(SCRIPT_DIR, f"{OUTPUT_PREFIX}_bam_file_list_{i}.txt")
        write_chunk_file(chunk, output_path)
        print(f"  Wrote {output_path} ({len(chunk)} patients)")

    print("Done.")


if __name__ == "__main__":
    main()
