source("phase2_expanded_mr.R")

chip <- data.frame(
  SNPID = "rs1", REF = "A", ALT = "G", BETA = 0.2, SE = 0.1, P = 0.05, N = 19646
)
mapped <- map_chip_mgi_outcome(chip)

stopifnot(mapped$SNP[[1]] == "rs1")
stopifnot(mapped$effect_allele[[1]] == "G")
stopifnot(mapped$other_allele[[1]] == "A")
stopifnot(!("eaf" %in% names(mapped)))

combined <- combine_result_tables(
  data.frame(method = "IVW", b = 0.1),
  data.frame(method = "Egger", b = 0.2)
)
stopifnot(nrow(combined) == 2)
stopifnot(identical(combined$method, c("IVW", "Egger")))

draw_a <- with_analysis_seed(20260808, runif(3))
draw_b <- with_analysis_seed(20260808, runif(3))
stopifnot(identical(draw_a, draw_b))

cat("expanded MR mapping tests passed\n")
