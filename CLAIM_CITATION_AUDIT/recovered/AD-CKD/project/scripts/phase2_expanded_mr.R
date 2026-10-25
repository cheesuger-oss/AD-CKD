#!/usr/bin/env Rscript

map_chip_mgi_outcome <- function(x) {
  data.frame(
    SNP = x$SNPID,
    beta = x$BETA,
    se = x$SE,
    effect_allele = x$ALT,
    other_allele = x$REF,
    pval = x$P,
    samplesize = x$N,
    stringsAsFactors = FALSE
  )
}

combine_result_tables <- function(...) {
  rows <- list(...)
  rownames_out <- NULL
  out <- do.call(rbind, rows)
  rownames(out) <- rownames_out
  out
}

with_analysis_seed <- function(seed, expression) {
  set.seed(seed)
  force(expression)
}

run_expanded_mr <- function(root) {
  root <- normalizePath(root)
  lib <- file.path(root, "project", "Rlib-4.5")
  renv_lib <- file.path(root, "renv", "library", "macos", "R-4.5", "aarch64-apple-darwin20")
  .libPaths(c(lib, renv_lib, .libPaths()))

  suppressPackageStartupMessages({
    library(data.table)
    library(TwoSampleMR)
    library(ggplot2)
  })

  datadir <- file.path(root, "project", "data", "gwas", "processed")
  outdir <- file.path(root, "project", "results", "mr", "expanded_outcomes")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  egfr <- fread(file.path(datadir, "egfr_clumped_instruments.tsv"))
  egfr[, Phenotype := "eGFRcrea (Stanzick 2021)"]
  exposure <- format_data(as.data.frame(egfr), type = "exposure",
    snp_col = "RSID", beta_col = "Effect", se_col = "StdErr.GC",
    effect_allele_col = "Allele1", other_allele_col = "Allele2",
    eaf_col = "Freq1", pval_col = "P.value.GC", samplesize_col = "n",
    chr_col = "chr", pos_col = "pos", phenotype_col = "Phenotype")
  exposure$exposure <- "eGFRcrea (Stanzick 2021)"

  fg <- fread(file.path(datadir, "finngen_thaortaneur_instrument_rows.tsv"))
  fg <- fg[!duplicated(rsids)]
  fg[, Phenotype := "Thoracic aortic aneurysm (FinnGen R13 I9_THAORTANEUR)"]
  fg_outcome <- format_data(as.data.frame(fg), type = "outcome",
    snp_col = "rsids", beta_col = "beta", se_col = "sebeta",
    effect_allele_col = "alt", other_allele_col = "ref",
    eaf_col = "af_alt", pval_col = "pval", phenotype_col = "Phenotype")
  fg_outcome$outcome <- fg$Phenotype[[1]]

  chip_raw <- fread(file.path(datadir, "chip_mgi_taa_instrument_rows.tsv"))
  chip_raw <- chip_raw[!duplicated(SNPID)]
  chip <- map_chip_mgi_outcome(as.data.frame(chip_raw))
  chip$Phenotype <- "Thoracic aortic aneurysm (CHIP-MGI 2021)"
  chip_outcome <- format_data(chip, type = "outcome",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    pval_col = "pval", samplesize_col = "samplesize", phenotype_col = "Phenotype")
  chip_outcome$outcome <- chip$Phenotype[[1]]

  analyse <- function(outcome, prefix, matches) {
    harm <- harmonise_data(exposure, outcome, action = 2)
    harm <- harm[harm$mr_keep, ]
    harm$F_statistic <- (harm$beta.exposure / harm$se.exposure)^2
    methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
    result <- with_analysis_seed(20260808, mr(harm, method_list = methods))
    result$OR <- exp(result$b)
    result$OR_lci95 <- exp(result$b - 1.96 * result$se)
    result$OR_uci95 <- exp(result$b + 1.96 * result$se)
    heterogeneity <- mr_heterogeneity(harm)
    pleiotropy <- mr_pleiotropy_test(harm)
    single <- mr_singlesnp(harm)

    fwrite(harm, file.path(outdir, paste0(prefix, "_harmonised.csv")))
    fwrite(result, file.path(outdir, paste0(prefix, "_mr_results.csv")))
    fwrite(heterogeneity, file.path(outdir, paste0(prefix, "_heterogeneity.csv")))
    fwrite(pleiotropy, file.path(outdir, paste0(prefix, "_egger_intercept.csv")))
    ggsave(file.path(outdir, paste0(prefix, "_scatter.pdf")),
      mr_scatter_plot(result, harm)[[1]], width = 8, height = 6)

    qc <- data.frame(
      outcome = unique(harm$outcome), source_matches = matches,
      harmonised_kept = nrow(harm), mean_F = mean(harm$F_statistic),
      min_F = min(harm$F_statistic), stringsAsFactors = FALSE
    )
    list(result = result, qc = qc, heterogeneity = heterogeneity, pleiotropy = pleiotropy)
  }

  fg_res <- analyse(fg_outcome, "finngen_thaortaneur", nrow(fg))
  chip_res <- analyse(chip_outcome, "chip_mgi_taa", nrow(chip_raw))
  combined_results <- combine_result_tables(fg_res$result, chip_res$result)
  fwrite(combined_results,
    file.path(outdir, "expanded_outcome_mr_results.csv"))
  fwrite(rbind(fg_res$qc, chip_res$qc), file.path(outdir, "expanded_outcome_qc.csv"))
  print(rbind(fg_res$qc, chip_res$qc))
  print(combined_results)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  run_expanded_mr(if (length(args)) args[[1]] else ".")
}
