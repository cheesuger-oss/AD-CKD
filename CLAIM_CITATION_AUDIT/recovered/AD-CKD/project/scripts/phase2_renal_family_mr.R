#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
lib <- file.path(root, "project", "Rlib-4.5")
renv_lib <- file.path(root, "renv", "library", "macos", "R-4.5", "aarch64-apple-darwin20")
.libPaths(c(lib, renv_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

ids <- c("serum_cystatin_c", "egfrcys", "bun_ea", "ckd_ea", "microalbuminuria", "rapid3", "ckdi25", "uacr_ea")
base <- file.path(root, "project", "data", "gwas", "processed", "renal_family")
outdir <- file.path(root, "project", "results", "mr", "renal_family")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fg <- fread(file.path(base, "finngen_aortdis_union_rows.tsv"))
fg <- fg[!duplicated(rsids)]
fg[, Phenotype := "Aortic dissection (FinnGen R13 I9_AORTDIS)"]

all_results <- list()
all_qc <- list()
all_het <- list()
all_pleio <- list()

for (id in ids) {
  cat("processing", id, "\n")
  flush.console()
  x <- fread(file.path(base, paste0(id, "_clumped_instruments.tsv")))
  x[, Phenotype := paste0(trait, " [", id, "]")]
  exp_dat <- format_data(
    as.data.frame(x), type = "exposure", snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele", eaf_col = "eaf",
    pval_col = "pval", samplesize_col = "n", chr_col = "chr", pos_col = "pos",
    phenotype_col = "Phenotype"
  )
  exp_dat$exposure <- id
  y <- fg[rsids %chin% x$SNP]
  out_dat <- format_data(
    as.data.frame(y), type = "outcome", snp_col = "rsids", beta_col = "beta", se_col = "sebeta",
    effect_allele_col = "alt", other_allele_col = "ref", eaf_col = "af_alt", pval_col = "pval",
    phenotype_col = "Phenotype"
  )
  out_dat$outcome <- "Aortic dissection (FinnGen R13 I9_AORTDIS)"
  harm <- harmonise_data(exp_dat, out_dat, action = 2)
  harm <- harm[harm$mr_keep, ]
  harm$F_statistic <- (harm$beta.exposure / harm$se.exposure)^2
  fwrite(harm, file.path(outdir, paste0("harmonised_", id, ".tsv")), sep = "\t")

  nsnp <- nrow(harm)
  methods <- if (nsnp >= 3) {
    c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
  } else if (nsnp == 2) {
    c("mr_ivw")
  } else if (nsnp == 1) {
    c("mr_wald_ratio")
  } else character()
  if (length(methods)) {
    res <- mr(harm, method_list = methods)
    res$exposure_id <- id
    res$OR <- exp(res$b)
    res$OR_lci95 <- exp(res$b - 1.96 * res$se)
    res$OR_uci95 <- exp(res$b + 1.96 * res$se)
    all_results[[id]] <- res
  }
  if (nsnp >= 2) {
    h <- tryCatch(mr_heterogeneity(harm), error = function(e) data.frame())
    if (nrow(h)) { h$exposure_id <- id; all_het[[id]] <- h }
  }
  if (nsnp >= 3) {
    p <- tryCatch(mr_pleiotropy_test(harm), error = function(e) data.frame())
    if (nrow(p)) { p$exposure_id <- id; all_pleio[[id]] <- p }
  }
  all_qc[[id]] <- data.frame(
    exposure_id = id, clumped = nrow(x), outcome_matched = nrow(y), harmonised_kept = nsnp,
    mean_F = if (nsnp) mean(harm$F_statistic) else NA_real_,
    min_F = if (nsnp) min(harm$F_statistic) else NA_real_, stringsAsFactors = FALSE
  )
}

results <- rbindlist(all_results, fill = TRUE)
qc <- rbindlist(all_qc, fill = TRUE)
het <- rbindlist(all_het, fill = TRUE)
pleio <- rbindlist(all_pleio, fill = TRUE)

results[, primary_ivw := method == "Inverse variance weighted"]
results[, q_primary_family := NA_real_]
idx <- which(results$primary_ivw)
results$q_primary_family[idx] <- p.adjust(results$pval[idx], method = "BH")
results[, evidence_grade := fifelse(primary_ivw & q_primary_family < 0.05, "FDR_significant", "not_FDR_significant")]

fwrite(results, file.path(outdir, "renal_family_mr_results.tsv"), sep = "\t")
fwrite(qc, file.path(outdir, "renal_family_qc.tsv"), sep = "\t")
fwrite(het, file.path(outdir, "renal_family_heterogeneity.tsv"), sep = "\t")
fwrite(pleio, file.path(outdir, "renal_family_egger_intercept.tsv"), sep = "\t")
print(qc)
print(results[primary_ivw == TRUE, .(exposure_id, nsnp, b, se, pval, OR, OR_lci95, OR_uci95, q_primary_family, evidence_grade)])
