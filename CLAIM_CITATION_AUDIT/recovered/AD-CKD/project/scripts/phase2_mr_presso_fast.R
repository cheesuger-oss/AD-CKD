#!/usr/bin/env Rscript

presso_observed_rss <- function(beta_x, beta_y, se_y) {
  weights <- 1 / se_y^2
  sum_xy <- sum(weights * beta_x * beta_y)
  sum_xx <- sum(weights * beta_x^2)
  slopes <- (sum_xy - weights * beta_x * beta_y) /
    (sum_xx - weights * beta_x^2)
  residuals <- beta_y - slopes * beta_x
  list(rss = sum(weights * residuals^2), slopes = slopes, residuals = residuals)
}


mr_presso_fast_univariable <- function(
  beta_x, beta_y, se_x, se_y, distributions = 1000L, seed = 20260808L
) {
  stopifnot(length(beta_x) == length(beta_y), length(beta_x) == length(se_x),
    length(beta_x) == length(se_y), distributions > length(beta_x))
  observed <- presso_observed_rss(beta_x, beta_y, se_y)
  predicted_y <- observed$slopes * beta_x
  simulated_rss <- numeric(distributions)
  simulated_residual_sq <- matrix(NA_real_, nrow = distributions, ncol = length(beta_x))
  set.seed(seed)
  for (iteration in seq_len(distributions)) {
    random_x <- rnorm(length(beta_x), beta_x, se_x)
    random_y <- rnorm(length(beta_y), predicted_y, se_y)
    random <- presso_observed_rss(random_x, random_y, se_y)
    simulated_rss[[iteration]] <- random$rss
    simulated_residual_sq[iteration, ] <- random$residuals^2
  }
  outlier_raw <- colMeans(sweep(simulated_residual_sq, 2, observed$residuals^2, `>`))
  outlier_bonferroni <- pmin(outlier_raw * length(beta_x), 1)
  weights <- 1 / se_y^2
  fit <- lm(beta_y ~ -1 + beta_x, weights = weights)
  coefficient <- summary(fit)$coefficients[1, ]
  list(
    method = "MR-PRESSO univariable algorithm (vectorised exact leave-one-out RSS)",
    distributions = as.integer(distributions), seed = as.integer(seed),
    global_rss = observed$rss,
    global_p = mean(simulated_rss > observed$rss),
    outlier_p_raw = outlier_raw,
    outlier_p_bonferroni = outlier_bonferroni,
    outlier_indices = which(outlier_bonferroni <= 0.05),
    distortion_test = NA_real_,
    estimate = unname(coefficient[["Estimate"]]),
    se = unname(coefficient[["Std. Error"]]),
    pval = unname(coefficient[["Pr(>|t|)"]])
  )
}


run_phase2_presso <- function(root) {
  harm <- read.csv(file.path(root, "project", "data", "gwas", "processed", "harmonised_egfr_aortdis.csv"))
  result <- mr_presso_fast_univariable(
    harm$beta.exposure, harm$beta.outcome, harm$se.exposure, harm$se.outcome,
    distributions = 1000L, seed = 20260808L
  )
  outdir <- file.path(root, "project", "results", "mr")
  saveRDS(result, file.path(outdir, "mr_presso_1000.rds"))
  summary <- data.frame(
    distributions = result$distributions, seed = result$seed,
    global_rss = result$global_rss,
    global_p = result$global_p,
    outlier_count = length(result$outlier_indices),
    min_outlier_p_bonferroni = min(result$outlier_p_bonferroni),
    estimate = result$estimate, se = result$se, pval = result$pval,
    stringsAsFactors = FALSE
  )
  write.table(summary, file.path(outdir, "mr_presso_1000_summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  capture.output(str(result), file = file.path(outdir, "mr_presso_1000.txt"))
  print(summary)
  invisible(result)
}


if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
  run_phase2_presso(root)
}
