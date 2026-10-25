#!/usr/bin/env Rscript

source("phase2_mr_presso_fast.R")

x <- c(0.10, 0.20, 0.15, 0.08, 0.12)
y <- c(0.03, 0.06, 0.02, 0.01, 0.04)
sx <- rep(0.01, length(x))
sy <- rep(0.02, length(x))
w <- 1 / sy^2
brute_slopes <- vapply(seq_along(x), function(i) {
  sum(w[-i] * x[-i] * y[-i]) / sum(w[-i] * x[-i]^2)
}, numeric(1))
brute_rss <- sum(w * (y - brute_slopes * x)^2)

fast_rss <- presso_observed_rss(x, y, sy)
stopifnot(isTRUE(all.equal(fast_rss$rss, brute_rss, tolerance = 1e-12)))

a <- mr_presso_fast_univariable(x, y, sx, sy, distributions = 100, seed = 42)
b <- mr_presso_fast_univariable(x, y, sx, sy, distributions = 100, seed = 42)
stopifnot(identical(a, b))
stopifnot(a$distributions == 100L)
stopifnot(length(a$outlier_p_bonferroni) == length(x))
cat("fast MR-PRESSO regression test passed\n")
