# 3_GLM_4_2026

Runs a per-variant genetic association test for aneuploidy rate across imputed genotypes on chromosomes 1-22, then builds QQ/Manhattan/effect-size plots and a top-hits table from the results.

## Model

For each variant, `glm_quasibinomial.R` fits:

```
glm(cbind(aneuploid_egg_num, embryo_count - aneuploid_egg_num) ~
      mat_age + PC1 + PC2 + ... + PC8 + seq_date + geno,
    family = quasibinomial)
```

- **Unit of analysis**: patient (not embryo). `aneuploid_egg_num` / `embryo_count` are patient-level counts of aneuploid vs. total tested embryos (aggregated from embryo-level `CCS_interpretation == "Abnormal"` calls).
- **Predictor of interest**: `geno`, the imputed dosage (`DS` field) of a variant for that patient.
- **Covariates**: maternal age, 8 ancestry principal components (from LASER, `laser.SeqPC.coord`), and a numeric sequencing-batch date (`seq_date_numeric`, months since Jan 2021).
- **Explicitly excluded** (per comments in the scripts): `total_embryos` and `paternal_age` are not included as covariates in this model variant — the filenames/output paths (`asso_no_totalembryos/`, `_no_totalembryos`) reflect that this is one specific covariate configuration; other configurations may exist elsewhere in the project but are not in this folder.
- Variants are dropped per-chromosome if they are monomorphic or have fewer than 50 non-missing genotypes.
- Output per chromosome: `beta`, `se`, `t`, `p`, and dispersion for the `geno` term, plus a separate summary of the non-genotype covariate effect sizes (`covariate_contributions_chr{N}.csv`).

Genotype input is imputed dosage data filtered for imputation INFO ≥ 0.2 and minor allele frequency ≥ 0.05 (5,616 patients, chromosomes 1-22).

## Pipeline order

The scripts in this folder are SLURM job scripts (`#SBATCH` headers) or the underlying R/Python they call, meant to run on the Amarel cluster. Inferred order:

1. **`step0_maf_impute_filter.sh`** — generates 22 per-chromosome SLURM scripts (`chr_scripts/af_filt_chr{N}.slurm`) that filter the ligated/imputed BCF by `INFO>=0.2` and minor allele frequency `>=0.05`, writing `maffilt/5616pts_ligated_merged_chr{N}_maffilt.bcf` and extracting a genotype-dosage text matrix (`gt_5616pts_ligated_merged_chr{N}_maffilt.txt`) via `bcftools query -f '%ID [\t%DS]\n'`.
2. **`index_glm.sh`** — SLURM array job (chr 1-22) that runs `bcftools index` on each filtered BCF; needed before later `bcftools query -l` calls.
3. **`count_post_filtered_variants.sh`** — QC utility; for each chromosome, extracts `%ID`/`%INFO` from the filtered BCF and appends the post-filter variant count (`bcftools index -n`) to `variant_counts_post_filter_chr_ALL.txt`. Not required for the GLM itself.
4. **`step1_single_seqDate_first.sh`** → **`add_seqdate_column.py`** — adds a `sequencing_date` (`YYYY_MM`) and `seq_date_numeric` (months since Jan 2021) column to the embryo-level phenotype table (`dataset_ONE.csv` from `1_Phenotype_4_2026`), parsed from the sequencer output file path. Output: `dataset_ONE_with_seqdate.csv`. (See `prepare_dataset_one/` for a locally-run notebook version of this same step, used to produce the copy of `dataset_ONE_with_seqdate.csv` committed in this repo.)
5. **`step1_single_seqDate_second_merge.sh`** → **`merge_phenotype_single_seqdate.py`** — aggregates the embryo-level table to one row per patient (`aneuploid_egg_num`, `embryo_count`, `mean_maternal_age`, `seq_date_numeric`, etc.), then inner-joins in the 8 ancestry PCs from `laser.SeqPC.coord`. Output: `merged_phenotype_single_seqdate.csv` — the patient-level covariate table used directly by `glm_quasibinomial.R`.
6. **`redone_submit_SINGLE_SEQ_all_chrs_step3.sh`** → **`step3_generate_gt_meta_single_seqdate.py`** — SLURM submission loop (chr 1-22) that builds a combined wide table of patient metadata + per-variant genotype dosages (`gt_meta_no_totalembryos_chr{N}.csv`). Note: this output does not appear to be read by `glm_quasibinomial.R` (which reads the metadata CSV and the `gt_*.txt` dosage file separately/directly) — it looks like a parallel/alternate data-prep path rather than a hard dependency of the GLM step; kept here for reference but its exact downstream use wasn't confirmed from the code.
7. **GLM step** — three SLURM wrapper scripts all invoke the same `glm_quasibinomial.R <chr>`, split across chromosomes by expected runtime/memory:
   - **`run_BIG_R.glm.sh`** — array 1-5 (largest chromosomes), 128 GB mem, 128 hr.
   - **`run_R_glm.sh`** — array 6-21, 96 GB mem, 72 hr.
   - **`run_TEST_R_glm.sh`** — single run hardcoded to chromosome 22.
   - **`glm_quasibinomial.R`** — the core script (see Model above). Reads `merged_phenotype_single_seqdate.csv` and the chromosome's `gt_*.txt` dosage file in chunks of 10,000 variants, fits the quasibinomial GLM per variant in parallel (`parallel::mclapply`), and writes `asso_no_totalembryos/glm_quasibinomial_chr{N}.csv` (sorted by p-value) and `asso_no_totalembryos/covariate_contributions_chr{N}.csv`.
8. *(not scripted in this folder)* the 22 per-chromosome result CSVs must be concatenated into `asso_no_totalembryos/all_chromosomes_combined.csv` before plotting — no combine script was found here, so this step was presumably done manually or with an ad hoc command.
9. **`run_create_plots_no_total_embryos.sh`** → **`create_plots_no_totalembryos.py`** — loads `all_chromosomes_combined.csv`, computes the genomic inflation factor λ (median χ² of t-values / 0.4549364), and produces a QQ plot, a genome-wide Manhattan plot, per-chromosome Manhattan plots, effect-size plots, and a top-100-hits table. Significance thresholds used throughout: suggestive p < 1e-5, genome-wide p < 5e-8. Outputs are mirrored in `plots/` (see `plots/README.md`).
