# FINAL_gls_scripts

Current/active scripts for computing genotype likelihoods (GLs) from merged, sorted
patient BAMs, plus a compiled `GLIMPSE_chunk` binary used at the end of the
per-chromosome job. Runs on Amarel via SLURM.

## Scripts

- **`create_bam_file_list.sh`** — Builds the sorted-BAM file list consumed by
  `gl_calculation.sh`. Finds all `*_merged.sorted.bam` files under the BAM-merge
  output directory (`.../bam_merge_sort/outputs`) and writes their full paths, one
  per line, to `bam_file_list.sorted.txt`. Note in the script: "we need the SORTED &
  INDEXED bam files, for gls, not the regular filtered & merged files."

- **`gl_calculation.sh`** — Main GL-calculation SLURM job array (`--array=1-22%10`,
  one task per human autosome, max 10 concurrent). For each chromosome: runs
  `bcftools mpileup` (against the GRCh38 reference and 1000GP cleaned sites VCF, over
  all BAMs in `bam_file_list.sorted.txt`) piped into `bcftools call -Aim -C alleles`
  to produce `5616pts.chr<N>.vcf.gz`, indexes it, then runs `./GLIMPSE_chunk` on the
  1000GP sites VCF to produce `chunks.chr<N>.txt` for that chromosome (2Mb windows,
  0.2Mb buffer). Has file-existence checks before each step and tees stdout/stderr to
  per-chromosome log files. "5616pts" in the output filename reflects the patient
  count in this dataset.

- **`count_vcf_records_old_vs_new.sh`** — QC/validation script: for each chromosome,
  compares record counts between the old pipeline's GL output
  (`initial_pipeline/gls/calculated_gls/`) and this pipeline's output
  (`1_Final_Pipeline_4_2026/gls/outputs_calculated_gls/`) using `bcftools view -H |
  wc -l`, writing a MATCH/DIFFER/MISSING status table to a TSV.

- **`GLIMPSE_chunk`** — Compiled GLIMPSE2 binary used by `gl_calculation.sh` to chunk
  each chromosome's sites VCF into imputation-sized windows. Vendored executable, not
  source-controlled here beyond the binary itself.

## Data files

- **`bam_file_list.sorted.txt`** — Generated output of `create_bam_file_list.sh`:
  full paths to all merged+sorted BAMs used as GL input.
- **`old_bam_file_list.txt`** — An earlier/alternate version of the BAM list, kept
  for reference (likely from the prior pipeline run referenced in
  `count_vcf_records_old_vs_new.sh`).
