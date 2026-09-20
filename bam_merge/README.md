# BAM Merge Pipeline — Step 1

Merges per-embryo BAM files into per-patient merged BAMs using Dataset ONE filtering criteria.

## Directory Structure (on Amarel)

```
/path/to/project_root/
├── bam/                                    # Source BAM files (pre-mapped from fastq)
├── 1_Final_Pipeline_4_2026/
│   ├── dataset_ONE.csv                     # Filtered dataset from analysis notebooks
│   └── bam_merge_sort/                     # This pipeline step
│       ├── scripts/
│       │   ├── prepare_merge_inputs.py     # Generates chunked job array files from CSV
│       │   ├── merge_bam_array.sh          # SLURM job array template (samtools merge)
│       │   ├── submit_all.sh               # Runs prepare script + submits all chunks
│       │   └── DS_ONE_bam_file_list_*.txt  # Generated at runtime by prepare script
│       ├── outputs/                        # Merged BAMs: {patient_ID}_merged.bam
│       └── logs/                           # SLURM stdout/stderr logs
```

## What Each Script Does

- **prepare_merge_inputs.py** — Reads `dataset_ONE.csv`, converts `.fastq.gz` filenames to `.bam`,
  resolves full paths under the `bam/` directory, groups by patient ID, and splits into chunks
  of ≤1000 patients (Amarel array size limit). All paths are defined as variables at the top of the file.

- **merge_bam_array.sh** — SLURM job array script. Each task reads one patient line from a chunk file,
  validates all input BAMs exist, and runs `samtools merge`. Skips patients whose output already exists
  (safe for re-runs). Cleans up partial output on failure.

- **submit_all.sh** — Wrapper script: Runs the Python prep script, then submits one SLURM array per chunk.
  Limits to 100 concurrent tasks per array.

## Notes

- Amarel max array size is 1001 — the pipeline auto-chunks into groups of ≤1000 patients.
- Each task requests 4 CPUs, 24GB RAM, 10-hour wall time on `main` partition.
- Merged BAMs land in `outputs/` as `{patient_ID}_merged.bam`.
