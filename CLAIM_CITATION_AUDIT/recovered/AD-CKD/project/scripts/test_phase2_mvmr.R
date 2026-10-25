#!/usr/bin/env Rscript

source("phase2_mvmr.R")

set.seed(42)
n <- 30
beta_x <- cbind(rnorm(n, 0.05, 0.02), rnorm(n, 0.03, 0.02))
se_x <- matrix(0.01, nrow = n, ncol = 2)
beta_y <- as.numeric(beta_x %*% c(-0.2, 0.1) + rnorm(n, 0, 0.01))
se_y <- rep(0.02, n)

result <- analyse_mvmr(beta_x, se_x, beta_y, se_y, paste0("rs", seq_len(n)))
stopifnot(identical(rownames(result$estimates), c("eGFRcrea", "SBP_IRNT")))
stopifnot(all(c("estimate", "se", "pval", "ci_low", "ci_high") %in% colnames(result$estimates)))
stopifnot(identical(colnames(result$conditional_f), c("eGFRcrea", "SBP_IRNT")))
stopifnot(nrow(result$conditional_f) == 1L)
combined <- combine_exposure_data(data.frame(a = 1, b = 2), data.frame(a = 3, c = 4))
stopifnot(nrow(combined) == 2L)
stopifnot(identical(names(combined), c("a", "b", "c")))
mock_mvdat <- list(exposure_beta = matrix(1:4, nrow = 2, dimnames = list(c("rsA", "rsB"), NULL)))
stopifnot(identical(extract_mvdat_rsids(mock_mvdat), c("rsA", "rsB")))
cat("phase2_mvmr regression test passed\n")
