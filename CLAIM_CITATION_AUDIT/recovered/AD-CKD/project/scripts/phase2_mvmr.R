#!/usr/bin/env Rscript

.phase2_root_guess <- if (dir.exists(file.path("..", "..", "project", "Rlib-4.5"))) {
  normalizePath(file.path("..", ".."))
} else {
  normalizePath(".")
}
.libPaths(c(
  file.path(.phase2_root_guess, "project", "Rlib-4.5"),
  file.path(.phase2_root_guess, "renv", "library", "macos", "R-4.5", "aarch64-apple-darwin20"),
  .libPaths()
))

analyse_mvmr <- function(beta_x, se_x, beta_y, se_y, rsids) {
  if (!requireNamespace("MVMR", quietly = TRUE)) stop("MVMR package is required")
  formatted <- MVMR::format_mvmr(
    BXGs = beta_x, BYG = beta_y, seBXGs = se_x, seBYG = se_y, RSID = rsids
  )
  estimate_matrix <- suppressWarnings(MVMR::ivw_mvmr(formatted, gencov = 0))
  estimates <- data.frame(
    estimate = estimate_matrix[, "Estimate"],
    se = estimate_matrix[, "Std. Error"],
    pval = estimate_matrix[, "Pr(>|t|)"],
    row.names = c("eGFRcrea", "SBP_IRNT")
  )
  estimates$ci_low <- estimates$estimate - 1.96 * estimates$se
  estimates$ci_high <- estimates$estimate + 1.96 * estimates$se
  conditional_f <- suppressWarnings(MVMR::strength_mvmr(formatted, gencov = 0))
  colnames(conditional_f) <- c("eGFRcrea", "SBP_IRNT")
  q_stat <- suppressWarnings(MVMR::pleiotropy_mvmr(formatted, gencov = 0))
  list(estimates = estimates, conditional_f = conditional_f, q_stat = q_stat, formatted = formatted)
}


combine_exposure_data <- function(...) {
  as.data.frame(data.table::rbindlist(list(...), use.names = TRUE, fill = TRUE))
}


extract_mvdat_rsids <- function(mvdat) {
  rownames(mvdat$exposure_beta)
}


run_phase2_mvmr <- function(root) {
  lib <- file.path(root, "project", "Rlib-4.5")
  renv_lib <- file.path(root, "renv", "library", "macos", "R-4.5", "aarch64-apple-darwin20")
  .libPaths(c(lib, renv_lib, .libPaths()))
  suppressPackageStartupMessages({
    library(data.table)
    library(TwoSampleMR)
    library(MVMR)
  })
  datadir <- file.path(root, "project", "data", "gwas", "processed")
  outdir <- file.path(root, "project", "results", "mr", "mvmr")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  selected <- scan(file.path(datadir, "mvmr_joint_clumped_rsids.txt"), what = "character", quiet = TRUE)
  joint <- fread(file.path(datadir, "mvmr_joint_clump_input.tsv"))
  joint <- joint[SNP %chin% selected, .(SNP, CHR = as.character(CHR), BP = as.character(BP))]

  egfr <- fread(file.path(datadir, "mvmr_egfr_rows.tsv"))
  egfr <- egfr[RSID %chin% selected]
  setorder(egfr, P.value.GC)
  egfr <- egfr[!duplicated(RSID)]
  egfr[, Phenotype := "eGFRcrea"]

  sbp <- fread(file.path(datadir, "mvmr_sbp_rows.tsv"))
  sbp[, `:=`(chr = as.character(chr), pos = as.character(pos))]
  sbp <- merge(sbp, joint, by.x = c("chr", "pos"), by.y = c("CHR", "BP"), all = FALSE)
  setorder(sbp, -neglog10_pval_EUR)
  sbp <- sbp[!duplicated(SNP)]
  sbp[, `:=`(pval = 10^(-neglog10_pval_EUR), Phenotype = "SBP_IRNT")]

  outcome <- fread(file.path(datadir, "mvmr_aortdis_rows.tsv"))
  outcome[, SNP := vapply(strsplit(rsids, ",", fixed = TRUE), function(ids) {
    hit <- intersect(ids, selected)
    if (length(hit)) hit[[1]] else NA_character_
  }, character(1))]
  outcome <- outcome[!is.na(SNP)]
  setorder(outcome, pval)
  outcome <- outcome[!duplicated(SNP)]
  outcome[, Phenotype := "Aortic dissection (FinnGen R13 I9_AORTDIS)"]

  egfr_dat <- format_data(as.data.frame(egfr), type = "exposure",
    snp_col = "RSID", beta_col = "Effect", se_col = "StdErr.GC",
    effect_allele_col = "Allele1", other_allele_col = "Allele2",
    eaf_col = "Freq1", pval_col = "P.value.GC", samplesize_col = "n",
    chr_col = "chr", pos_col = "pos", phenotype_col = "Phenotype")
  sbp_dat <- format_data(as.data.frame(sbp), type = "exposure",
    snp_col = "SNP", beta_col = "beta_EUR", se_col = "se_EUR",
    effect_allele_col = "alt", other_allele_col = "ref",
    eaf_col = "af_EUR", pval_col = "pval", chr_col = "chr", pos_col = "pos",
    phenotype_col = "Phenotype")
  outcome_dat <- format_data(as.data.frame(outcome), type = "outcome",
    snp_col = "SNP", beta_col = "beta", se_col = "sebeta",
    effect_allele_col = "alt", other_allele_col = "ref",
    eaf_col = "af_alt", pval_col = "pval", phenotype_col = "Phenotype")

  exposure_dat <- combine_exposure_data(egfr_dat, sbp_dat)
  mvdat <- mv_harmonise_data(exposure_dat, outcome_dat, harmonise_strictness = 2)
  result <- analyse_mvmr(
    mvdat$exposure_beta, mvdat$exposure_se, mvdat$outcome_beta,
    mvdat$outcome_se, extract_mvdat_rsids(mvdat)
  )

  fwrite(as.data.table(result$formatted), file.path(outdir, "harmonised_mvmr.tsv"), sep = "\t")
  fwrite(as.data.table(result$estimates, keep.rownames = "exposure"),
    file.path(outdir, "mvmr_direct_effects.tsv"), sep = "\t")
  fwrite(as.data.table(result$conditional_f, keep.rownames = "metric"),
    file.path(outdir, "conditional_f.tsv"), sep = "\t")
  fwrite(data.table(Q = result$q_stat$Qstat, pval = result$q_stat$Qpval),
    file.path(outdir, "pleiotropy_q.tsv"), sep = "\t")
  qc <- data.table(
    joint_clumped = length(selected), egfr_rows = nrow(egfr), sbp_rows = nrow(sbp),
    outcome_rows = nrow(outcome), harmonised_complete = nrow(mvdat$exposure_beta),
    gencov_assumption = "0 (approximation; exposure-association covariance unavailable)"
  )
  fwrite(qc, file.path(outdir, "mvmr_qc_summary.tsv"), sep = "\t")
  print(qc)
  print(result$estimates)
  invisible(result)
}


if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
  run_phase2_mvmr(root)
}
