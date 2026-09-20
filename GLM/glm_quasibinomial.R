userlib <- Sys.getenv("R_LIBS_USER")
if (nzchar(userlib)) .libPaths(c(userlib, .libPaths()))
library(data.table)
library(parallel)

# Cluster path root -- update to match your own environment.
work_dir <- "/path/to/project/imputation/outputs/ligate/maffilt"
setwd(work_dir)

chr <- as.integer(commandArgs(trailingOnly = TRUE)[1])
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))

cat(sprintf("Using %d cores for chromosome %d\n", n_cores, chr))
cat("Model: NO total_embryos, NO paternal_age\n")
cat("Covariates: mat_age + PC1-8 + seq_date + geno\n")

# Load metadata
metadata <- fread("merged_phenotype_single_seqdate.csv")
metadata$patient_ID <- as.character(metadata$patient_ID)

# Get sample IDs
sample_file <- sprintf("sample_list_chr%d.txt", chr)
if (!file.exists(sample_file)) {
    system(sprintf("bcftools query -l 5616pts_ligated_merged_chr%d_maffilt.bcf > %s", chr, sample_file))
}

sample_ids <- fread(sample_file, header = FALSE)[[1]]
sample_ids <- sub(".*/(\\d+)_.*", "\\1", sample_ids)

meta_matched <- metadata[match(sample_ids, metadata$patient_ID), ]

# Function to process one variant (NO total_embryos)
process_variant <- function(variant_id, geno_row, meta_matched) {
    
    geno <- as.numeric(geno_row)
    if (length(unique(geno[!is.na(geno)])) <= 1 || sum(!is.na(geno)) < 50) return(NULL)
    
    df <- data.frame(
        y = meta_matched$aneuploid_egg_num,
        n = meta_matched$embryo_count,
        mat_age = meta_matched$mean_maternal_age,
        PC1 = meta_matched$PC1, 
        PC2 = meta_matched$PC2, 
        PC3 = meta_matched$PC3, 
        PC4 = meta_matched$PC4,
        PC5 = meta_matched$PC5, 
        PC6 = meta_matched$PC6, 
        PC7 = meta_matched$PC7, 
        PC8 = meta_matched$PC8,
        seq_date = meta_matched$seq_date_numeric,
        geno = geno
    )
    df <- df[!is.na(df$geno), ]
    
    # Model WITHOUT total_embryos and paternal_age
    model <- glm(cbind(y, n - y) ~ mat_age + 
                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                 seq_date + geno,
                 data = df, family = quasibinomial)
    
    coef_summary <- summary(model)$coefficients
    
    result <- list(
        main = data.frame(
            variant = variant_id,
            beta = coef_summary["geno", "Estimate"],
            se = coef_summary["geno", "Std. Error"],
            t = coef_summary["geno", "t value"],
            p = coef_summary["geno", "Pr(>|t|)"],
            dispersion = summary(model)$dispersion,
            n = nrow(df)
        ),
        covariates = {
            coef_df <- as.data.frame(coef_summary)
            coef_df$covariate <- rownames(coef_df)
            coef_df <- coef_df[!coef_df$covariate %in% c("(Intercept)", "geno"), ]
            coef_df$variant <- variant_id
            coef_df
        }
    )
    
    return(result)
}

# Process in chunks
gt_file <- sprintf("gt_5616pts_ligated_merged_chr%d_maffilt.txt", chr)
n_variants <- as.integer(system(sprintf("wc -l < %s", gt_file), intern = TRUE))
chunk_size <- 10000

cat(sprintf("Total variants: %d, processing in chunks of %d\n", n_variants, chunk_size))

all_results <- list()
all_covariates <- list()

for (chunk_start in seq(1, n_variants, chunk_size)) {
    
    skip_rows <- chunk_start - 1
    nrows_to_read <- min(chunk_size, n_variants - chunk_start + 1)
    
    cat(sprintf("Reading chunk: variants %d-%d (%d/%d)\n", 
                chunk_start, chunk_start + nrows_to_read - 1, 
                chunk_start, n_variants))
    
    gt_chunk <- fread(gt_file, skip = skip_rows, nrows = nrows_to_read, header = FALSE)
    variant_ids_chunk <- gt_chunk[[1]]
    gt_chunk <- gt_chunk[, -1]
    
    cat(sprintf("Processing %d variants with %d cores...\n", nrow(gt_chunk), n_cores))
    
    chunk_results <- mclapply(1:nrow(gt_chunk), function(j) {
        process_variant(variant_ids_chunk[j], gt_chunk[j, ], meta_matched)
    }, mc.cores = n_cores)
    
    chunk_results <- chunk_results[!sapply(chunk_results, is.null)]
    
    if (length(chunk_results) > 0) {
        all_results[[length(all_results) + 1]] <- rbindlist(lapply(chunk_results, function(x) x$main))
        all_covariates[[length(all_covariates) + 1]] <- rbindlist(lapply(chunk_results, function(x) x$covariates))
    }
    
    rm(gt_chunk, chunk_results)
    gc()
    
    cat(sprintf("Chunk complete. Total results so far: %d\n", 
                sum(sapply(all_results, nrow))))
}

# Combine all results
cat("Combining all results...\n")
results_df <- rbindlist(all_results)
covariate_all <- rbindlist(all_covariates)

results_df <- results_df[order(results_df$p), ]
fwrite(results_df, sprintf("asso_no_totalembryos/glm_quasibinomial_chr%d.csv", chr))

covariate_summary <- covariate_all[, .(
    mean_beta = mean(Estimate),
    median_beta = median(Estimate),
    sd_beta = sd(Estimate),
    mean_abs_beta = mean(abs(Estimate)),
    n_sig_0.05 = sum(`Pr(>|t|)` < 0.05),
    n_sig_1e5 = sum(`Pr(>|t|)` < 1e-5),
    n_sig_5e8 = sum(`Pr(>|t|)` < 5e-8),
    pct_sig_0.05 = 100 * sum(`Pr(>|t|)` < 0.05) / .N
), by = covariate]

covariate_summary <- covariate_summary[order(-mean_abs_beta)]
fwrite(covariate_summary, sprintf("asso_no_totalembryos/covariate_contributions_chr%d.csv", chr))

cat("Done!\n")
