#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
.libPaths(c(
  file.path(root, "project", "Rlib-4.5"),
  file.path(root, "renv", "library", "macos", "R-4.5", "aarch64-apple-darwin20"),
  .libPaths()
))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(scales) })

mr_dir <- file.path(root, "project", "results", "mr")
gen_dir <- file.path(root, "project", "results", "shared_genetics")
fig_dir <- file.path(root, "project", "results", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_classic(base_size = 8.2) + theme(
  plot.title = element_text(face = "bold", size = 9.2),
  plot.subtitle = element_text(size = 7.2),
  axis.title = element_text(face = "bold"),
  plot.margin = margin(5, 7, 5, 5)
)

# Panel A: all nine prespecified forward-MR estimates.
main <- read.csv(file.path(mr_dir, "main_mr_results.csv"), check.names = FALSE)
main <- main[main$method == "Inverse variance weighted", ][1, ]
renal <- read.delim(file.path(mr_dir, "renal_family", "renal_family_mr_results.tsv"), check.names = FALSE)
renal <- renal[renal$primary_ivw, ]
family_q <- unique(renal$q_primary_family)
stopifnot(length(family_q) == 1L)

forest <- rbind(
  data.frame(exposure_id = "egfrcrea", nsnp = main$nsnp, OR = main$OR,
             low = main$OR_lci95, high = main$OR_uci95, p = main$pval, q = family_q),
  data.frame(exposure_id = renal$exposure_id, nsnp = renal$nsnp, OR = renal$OR,
             low = renal$OR_lci95, high = renal$OR_uci95, p = renal$pval,
             q = renal$q_primary_family)
)
labels <- c(
  egfrcrea = "eGFRcrea", egfrcys = "eGFRcys", serum_cystatin_c = "Serum cystatin C",
  bun_ea = "Blood urea nitrogen", uacr_ea = "UACR", microalbuminuria = "Microalbuminuria",
  ckd_ea = "Chronic kidney disease", rapid3 = "Rapid3 decline", ckdi25 = "CKDi25 decline"
)
order_ids <- c("egfrcrea", "egfrcys", "serum_cystatin_c", "bun_ea", "uacr_ea",
               "microalbuminuria", "ckd_ea", "rapid3", "ckdi25")
forest$label <- sprintf("%s  (n=%d)", labels[forest$exposure_id], forest$nsnp)
forest$label <- factor(forest$label, levels = rev(sprintf("%s  (n=%d)", labels[order_ids], forest$nsnp[match(order_ids, forest$exposure_id)])))

p_a <- ggplot(forest, aes(OR, label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = low, xmax = high), width = 0.16, linewidth = 0.48) +
  geom_point(shape = 21, size = 2.4, fill = "#0072B2") +
  scale_x_log10(limits = c(0.025, 12), breaks = c(0.03, 0.1, 0.3, 1, 3, 10)) +
  labs(
    title = "No renal-trait estimate\npassed correction",
    subtitle = sprintf("Nine primary IVW estimates\nAll BH q=%.3f", family_q),
    x = "Odds ratio for aortic dissection (log scale)", y = NULL
  ) + theme_pub

# Panel B: MVMR conditional effects on their native beta scale.
mvmr <- read.delim(file.path(mr_dir, "mvmr", "mvmr_direct_effects.tsv"), check.names = FALSE)
cf <- read.delim(file.path(mr_dir, "mvmr", "conditional_f.tsv"), check.names = FALSE)
cf_values <- as.numeric(cf[1, c("eGFRcrea", "SBP_IRNT")])
mvmr$label <- c("eGFRcrea", "Systolic blood pressure (IRNT)")
mvmr$conditional_f <- cf_values[match(mvmr$exposure, c("eGFRcrea", "SBP_IRNT"))]
mvmr$label_n <- sprintf("%s  (conditional F=%.1f)", mvmr$label, mvmr$conditional_f)
mvmr$label_n <- factor(mvmr$label_n, levels = rev(mvmr$label_n))
p_b <- ggplot(mvmr, aes(estimate, label_n)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.12, linewidth = 0.5) +
  geom_point(shape = 21, size = 2.5, fill = "#D55E00") +
  labs(
    title = "Conditional multivariable\nMR estimates",
    subtitle = "Mutually adjusted effects\nNot a mediation test",
    x = "Conditional beta (95% CI)", y = NULL
  ) + theme_pub

# Panel C: compact status matrix for complementary analyses.
ldsc <- read.delim(file.path(gen_dir, "ldsc", "univariate_h2_qc.tsv"), check.names = FALSE)
aad <- ldsc[ldsc$trait == "AAD", ][1, ]
reverse <- read.delim(file.path(gen_dir, "coloc", "reverse_mr_results.tsv"), check.names = FALSE)
estimable <- reverse[reverse$status == "ESTIMABLE", ]
under <- reverse[reverse$status != "ESTIMABLE", ]
coloc <- read.delim(file.path(gen_dir, "coloc", "coloc_primary_summary.tsv"), check.names = FALSE)
coloc_priors <- read.delim(file.path(gen_dir, "coloc", "coloc_abf_all_priors.tsv"), check.names = FALSE)
chr14 <- coloc[coloc$window_id %in% c("egfrcys_W042", "cystatin_c_W097"), ]
chr20 <- coloc[coloc$window_id %in% c("egfrcys_W087", "cystatin_c_W188"), ]

status <- data.frame(
  analysis = factor(c("LD-score regression", "Colocalization", "Reverse MR"),
                    levels = rev(c("LD-score regression", "Colocalization", "Reverse MR"))),
  status = c("LDSC GATE", "PRIOR SENSITIVE", "REVERSE MR LIMITED"),
  detail = c(
    sprintf("AAD h2 Z=%.2f < 4; genetic correlation was not estimated", aad$z),
    sprintf("chr14 PP.H4 %.3f / %.3f was prior-sensitive; chr20 PP.H3 %.3f / %.3f favoured distinct variants",
            chr14$PP.H4[match("egfrcys", chr14$trait)], chr14$PP.H4[match("cystatin_c", chr14$trait)],
            chr20$PP.H3[match("egfrcys", chr20$trait)], chr20$PP.H3[match("cystatin_c", chr20$trait)]),
    sprintf("%d outcomes each used 2 instruments; none passed BH; %d outcome was underidentified",
            nrow(estimable), nrow(under))
  ),
  fill = c("#E69F00", "#CC79A7", "#56B4E9")
)
status$card <- c(
  sprintf("LDSC GATE\nAAD h2 Z=%.2f < 4\nGenetic correlation was not estimated", aad$z),
  sprintf("PRIOR SENSITIVE\nchr14 PP.H4 %.3f / %.3f\nvaried with the prior\nchr20 PP.H3 %.3f / %.3f\nfavoured distinct variants",
          chr14$PP.H4[match("egfrcys", chr14$trait)], chr14$PP.H4[match("cystatin_c", chr14$trait)],
          chr20$PP.H3[match("egfrcys", chr20$trait)], chr20$PP.H3[match("cystatin_c", chr20$trait)]),
  sprintf("REVERSE MR LIMITED\n%d outcomes each used 2 instruments\nNone passed BH correction\n%d outcome was underidentified",
          nrow(estimable), nrow(under))
)
p_c <- ggplot(status, aes(y = analysis)) +
  geom_tile(aes(x = 0.5, fill = analysis), width = 0.94, height = 0.76,
            colour = "black", linewidth = 0.25) +
  geom_text(aes(x = 0.5, label = card), size = 1.95, lineheight = 0.9) +
  scale_fill_manual(values = setNames(status$fill, status$analysis), guide = "none") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "Complementary analyses\nconstrained interpretation",
    subtitle = "Distinct evidential states\nare not interchangeable",
    x = NULL, y = NULL
  ) +
  theme_void(base_size = 8.2) +
  theme(plot.title = element_text(face = "bold", size = 9.2),
        plot.subtitle = element_text(size = 7.2), plot.margin = margin(5, 7, 5, 5))

right_column <- p_b / p_c + plot_layout(heights = c(0.82, 1.18))
figure3 <- (p_a | right_column) +
  plot_layout(widths = c(1.08, 1.22)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 10))
  )

base <- file.path(fig_dir, "Figure3_genetic_evidence_boundary")
ggsave(paste0(base, ".pdf"), figure3, width = 183, height = 155, units = "mm",
       device = "pdf", bg = "white", useDingbats = FALSE)
ggsave(paste0(base, ".png"), figure3, width = 183, height = 155, units = "mm",
       dpi = 600, bg = "white")

write.table(forest, paste0(base, "_source_data_A_forward_MR.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(mvmr, paste0(base, "_source_data_B_MVMR.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
prior14 <- coloc_priors[coloc_priors$window_id %in% c("egfrcys_W042", "cystatin_c_W097"), ]
empty_source <- function(n) data.frame(
  analysis = rep(NA_character_, n), trait = rep(NA_character_, n), window_id = rep(NA_character_, n),
  p12 = rep(NA_real_, n), estimate = rep(NA_real_, n), se = rep(NA_real_, n),
  ci_low = rep(NA_real_, n), ci_high = rep(NA_real_, n), z = rep(NA_real_, n),
  pp_h3 = rep(NA_real_, n), pp_h4 = rep(NA_real_, n), n_instruments = rep(NA_real_, n),
  p = rep(NA_real_, n), q_bh = rep(NA_real_, n), heterogeneity_q = rep(NA_real_, n),
  heterogeneity_p = rep(NA_real_, n), status = rep(NA_character_, n), stringsAsFactors = FALSE
)
ldsc_source <- empty_source(1); ldsc_source$analysis <- "LD-score regression"
ldsc_source$trait <- "AAD"; ldsc_source$estimate <- aad$h2_observed; ldsc_source$se <- aad$se
ldsc_source$ci_low <- aad$ci_low; ldsc_source$ci_high <- aad$ci_high; ldsc_source$z <- aad$z
ldsc_source$status <- aad$qc_status
prior_source <- empty_source(nrow(prior14)); prior_source$analysis <- "Colocalization chr14 prior sensitivity"
prior_source$trait <- prior14$trait; prior_source$window_id <- prior14$window_id
prior_source$p12 <- prior14$p12; prior_source$pp_h3 <- prior14$PP.H3; prior_source$pp_h4 <- prior14$PP.H4
prior_source$status <- prior14$status
chr20_source <- empty_source(nrow(chr20)); chr20_source$analysis <- "Colocalization chr20 primary prior"
chr20_source$trait <- chr20$trait; chr20_source$window_id <- chr20$window_id
chr20_source$pp_h3 <- chr20$PP.H3; chr20_source$pp_h4 <- chr20$PP.H4
chr20_source$status <- chr20$classification
reverse_source <- empty_source(nrow(reverse)); reverse_source$analysis <- "Reverse MR"
reverse_source$trait <- reverse$trait; reverse_source$estimate <- reverse$beta; reverse_source$se <- reverse$se
reverse_source$ci_low <- reverse$ci_low; reverse_source$ci_high <- reverse$ci_high
reverse_source$n_instruments <- reverse$n_instruments; reverse_source$p <- reverse$p
reverse_source$q_bh <- reverse$q_bh; reverse_source$heterogeneity_q <- reverse$Q
reverse_source$heterogeneity_p <- reverse$Q_p; reverse_source$status <- reverse$status
status_source <- rbind(ldsc_source, prior_source, chr20_source, reverse_source)
write.table(status_source,
            paste0(base, "_source_data_C_status.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
source_index <- data.frame(
  panel = c("A", "B", "C"),
  file = c("Figure3_genetic_evidence_boundary_source_data_A_forward_MR.tsv",
           "Figure3_genetic_evidence_boundary_source_data_B_MVMR.tsv",
           "Figure3_genetic_evidence_boundary_source_data_C_status.tsv"),
  claim_ids = c("C4", "C5", "C12-C15")
)
write.table(source_index, paste0(base, "_source_data.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("FIGURE3_COMPLETE\n")
