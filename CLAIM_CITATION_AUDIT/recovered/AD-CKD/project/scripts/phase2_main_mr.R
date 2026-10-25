#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
lib <- file.path(root, "project", "Rlib-4.5")
renv_lib <- file.path(root, "renv", "library", "macos", "R-4.5", "aarch64-apple-darwin20")
.libPaths(c(lib, renv_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
  library(MRPRESSO)
  library(ggplot2)
})

downloads <- "inputs"
egfr_file <- file.path(downloads, "metal_eGFR_meta_ea1.TBL.map.annot.gc.gz")
out_file <- file.path(downloads, "finngen_R13_I9_AORTDIS.gz")
bfile <- file.path(root, "project", "data", "gwas", "ld", "EUR")
plink <- file.path(root, "project", "tools", "plink", "plink")
outdir <- file.path(root, "project", "results", "mr")
datadir <- file.path(root, "project", "data", "gwas", "processed")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)

egfr_subset <- file.path(datadir, "egfr_gws_candidates_v2.tsv")
egfr <- fread(egfr_subset,
  select = c("Allele1", "Allele2", "n", "Freq1", "Effect", "P.value.GC", "StdErr.GC", "chr", "pos", "RSID"))
egfr <- egfr[is.finite(P.value.GC) & P.value.GC < 5e-8 & grepl("^rs", RSID)]
setorder(egfr, P.value.GC)
egfr <- egfr[!duplicated(RSID)]
fwrite(egfr, file.path(datadir, "egfr_gws_candidates.tsv"), sep = "\t")

clumped_ids <- readLines(file.path(datadir, "egfr_clumped_rsids.txt"))
egfr_c <- egfr[match(clumped_ids, RSID), nomatch = 0]
egfr_c[, Phenotype := "eGFRcrea (Stanzick 2021)"]
fwrite(egfr_c, file.path(datadir, "egfr_clumped_instruments.tsv"), sep = "\t")

exp_dat <- format_data(as.data.frame(egfr_c), type = "exposure",
  snp_col = "RSID", beta_col = "Effect", se_col = "StdErr.GC",
  effect_allele_col = "Allele1", other_allele_col = "Allele2",
  eaf_col = "Freq1", pval_col = "P.value.GC", samplesize_col = "n",
  chr_col = "chr", pos_col = "pos", phenotype_col = "Phenotype")
exp_dat$exposure <- "eGFRcrea (Stanzick 2021)"

fg <- fread(file.path(datadir, "finngen_aortdis_instrument_rows_v2.tsv"),
  select = c("rsids", "ref", "alt", "pval", "beta", "sebeta", "af_alt"))
fg <- fg[rsids %chin% egfr_c$RSID]
fg <- fg[!duplicated(rsids)]
fg[, Phenotype := "Aortic dissection (FinnGen R13 I9_AORTDIS)"]
fwrite(fg, file.path(datadir, "finngen_aortdis_instrument_rows.tsv"), sep = "\t")

out_dat <- format_data(as.data.frame(fg), type = "outcome",
  snp_col = "rsids", beta_col = "beta", se_col = "sebeta",
  effect_allele_col = "alt", other_allele_col = "ref",
  eaf_col = "af_alt", pval_col = "pval", phenotype_col = "Phenotype")
out_dat$outcome <- "Aortic dissection (FinnGen R13 I9_AORTDIS)"

harm <- harmonise_data(exp_dat, out_dat, action = 2)
harm <- harm[harm$mr_keep, ]
harm$F_statistic <- (harm$beta.exposure / harm$se.exposure)^2
write.csv(harm, file.path(datadir, "harmonised_egfr_aortdis.csv"), row.names = FALSE)

methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
res <- mr(harm, method_list = methods)
res$OR <- exp(res$b)
res$OR_lci95 <- exp(res$b - 1.96 * res$se)
res$OR_uci95 <- exp(res$b + 1.96 * res$se)
write.csv(res, file.path(outdir, "main_mr_results.csv"), row.names = FALSE)
write.csv(mr_heterogeneity(harm), file.path(outdir, "heterogeneity.csv"), row.names = FALSE)
write.csv(mr_pleiotropy_test(harm), file.path(outdir, "egger_intercept.csv"), row.names = FALSE)
write.csv(mr_singlesnp(harm), file.path(outdir, "single_snp.csv"), row.names = FALSE)
write.csv(mr_leaveoneout(harm), file.path(outdir, "leave_one_out.csv"), row.names = FALSE)

ggsave(file.path(outdir, "scatter.pdf"), mr_scatter_plot(res, harm)[[1]], width = 8, height = 6)
ggsave(file.path(outdir, "forest.pdf"), mr_forest_plot(mr_singlesnp(harm))[[1]], width = 8, height = 9)
ggsave(file.path(outdir, "funnel.pdf"), mr_funnel_plot(mr_singlesnp(harm))[[1]], width = 8, height = 6)

set.seed(20260808)
presso <- tryCatch(mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE, data = harm, NbDistribution = 1000, SignifThreshold = 0.05),
  error = function(e) list(error = conditionMessage(e)))
capture.output(str(presso), file = file.path(outdir, "mr_presso.txt"))

summary <- data.frame(
  gws_candidates = nrow(egfr), clumped_instruments = nrow(egfr_c),
  outcome_matches = nrow(fg), harmonised_kept = nrow(harm),
  mean_F = mean(harm$F_statistic), min_F = min(harm$F_statistic),
  stringsAsFactors = FALSE)
write.csv(summary, file.path(outdir, "main_mr_qc_summary.csv"), row.names = FALSE)
print(summary)
print(res)
