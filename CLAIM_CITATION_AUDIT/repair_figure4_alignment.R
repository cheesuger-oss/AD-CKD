#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
.libPaths(c(
  file.path(root, "project", "Rlib-4.5"),
  file.path(root, "renv", "library", "macos", "R-4.5", "aarch64-apple-darwin20"),
  .libPaths()
))
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

analysis_dir <- file.path(root, "CLAIM_CITATION_AUDIT/recovered/AD-CKD/project/results/figures")
figure_dir <- file.path(root, "CLAIM_CITATION_AUDIT/aad_revision_candidate/figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_classic(base_size = 8.5) +
  theme(
    plot.title = element_text(face = "bold", size = 8.8),
    plot.subtitle = element_text(size = 7.5),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.margin = margin(5, 7, 5, 5)
  )

palette <- c(GSE146638 = "#009E73", GSE153434 = "#D55E00", GSE38750 = "#0072B2")
dataset_labels <- c(
  GSE146638 = "Uraemic calcification model",
  GSE153434 = "Human aortic dissection",
  GSE38750 = "Human uremic vascular"
)

# Panel A: standardized immune-programme effects with approximate 95% CIs.
immune <- read.delim(file.path(analysis_dir, "Figure4_immune_classifier_robustness_source_data_A_immune.tsv"), check.names = FALSE)
immune$dataset_label <- factor(dataset_labels[immune$dataset], levels = unname(dataset_labels))
immune$program_name <- factor(immune$program_name, levels = rev(unique(immune$program_name)))
pd <- position_dodge(width = 0.55)
p_a <- ggplot(immune, aes(standardized_effect, program_name, colour = dataset_label, group = dataset_label)) +
  geom_vline(xintercept = 0, colour = "grey65", linewidth = 0.35) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high, linetype = dataset_label), position = pd, width = 0, linewidth = 0.42) +
  geom_point(aes(shape = q_value < 0.05), position = pd, size = 2) +
  scale_colour_manual(values = setNames(palette[names(dataset_labels)], unname(dataset_labels))) +
  scale_linetype_manual(values = setNames(c('solid','dashed','dotted'),unname(dataset_labels)),name='Dataset') +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(
    title = "Immune-programme differences varied across datasets",
    subtitle = "Approximate Wald 95% CIs for uncorrected standardized mean differences",
    x = "Standardized difference", y = NULL, colour = "Dataset"
  ) + theme_pub + guides(colour = guide_legend(nrow = 2),linetype=guide_legend(nrow=2))

# Both layers must assign exactly the same y coordinate to each dataset/programme.
built_a <- ggplot_build(p_a)
stopifnot(identical(built_a$data[[2]]$y, built_a$data[[3]]$y))

# Panel B: observed internal discrimination against the complete random-gene null.
frozen_random <- read.delim(file.path(analysis_dir, "Figure4_immune_classifier_robustness_source_data_B_random_auc.tsv"), check.names = FALSE)
random_auc <- frozen_random[frozen_random$record_type == "size_matched_random_gene", ]
observed_auc <- frozen_random$auc[frozen_random$record_type == "observed_renal_derived"]
perfect_random_n <- sum(random_auc$auc == 1)
p_b <- ggplot(random_auc, aes(auc)) +
  geom_histogram(binwidth = 0.025, boundary = 0, fill = "#999999", colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = observed_auc, colour = "#D55E00", linewidth = 0.8) +
  annotate("text", x = observed_auc, y = Inf, label = sprintf("Observed OOF AUC = %.2f", observed_auc),
           colour = "#D55E00", hjust = 1.03, vjust = 1.5, size = 2.7) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(
    title = "Random gene sets frequently matched\nperfect internal separation",
    subtitle = "100 size-matched controls; fixed preprocessing",
    x = "Out-of-fold ROC AUC", y = "Random gene sets (n)"
  ) + theme_pub
p_b <- p_b + annotate(
  "text", x = 0.98, y = 84,
  label = sprintf("%d/100 controls had AUC=1.00\n95th percentile=1.00", perfect_random_n),
  hjust = 1, size = 2.6, colour = "#444444"
)

# Panel C: locked external validation, retaining the non-evaluable cohort visibly.
external <- read.delim(file.path(analysis_dir, "Figure4_immune_classifier_robustness_source_data_C_external.tsv"), check.names = FALSE,
                       na.strings = c("", "NA"))
external$dataset <- factor(external$dataset, levels = rev(c("GSE52093", "GSE98770", "GSE147026", "GSE190635")))
evaluated <- external[external$status == "EVALUATED", ]
not_ready <- external[external$status != "EVALUATED", ]
p_c <- ggplot(external, aes(auc, dataset)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60") +
  geom_errorbar(data = evaluated, aes(xmin = auc_ci_low, xmax = auc_ci_high), width = 0.16, linewidth = 0.5) +
  geom_point(data = evaluated, shape = 21, size = 2.6, fill = "#CC79A7") +
  geom_text(data = not_ready, aes(x = 0.02, label = "NOT EVALUABLE: no expression measurements"),
            hjust = 0, colour = "#666666", size = 2.6) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_discrete(limits = rev(c("GSE52093", "GSE98770", "GSE147026", "GSE190635"))) +
  labs(
    title = "The locked signature did not\ngeneralize consistently",
    subtitle = "Fixed weighted-rank scores; no cohort refitting",
    x = "ROC AUC (bootstrap 95% CI)", y = NULL
  ) + theme_pub

figure4 <- p_a / (p_b | p_c) +
  plot_layout(heights = c(1.2, 1), widths = c(0.9, 1.1)) +
  plot_annotation(
    tag_levels = "A",
    caption = paste(strwrap(paste(
      "Filled immune-effect points passed BH q<0.05 within the corresponding six-programme dataset family; open points did not.",
      "Immune scores represent expression programmes, not immune-cell proportions.",
      "Expression preprocessing used the full development cohort; feature restriction, scaling and tuning were fold-internal.",
      "The classifier analysis concerns retrospective tissue discrimination and does not evaluate future aortic-dissection risk."
    ), width = 150), collapse = "\n"),
    theme = theme(
      plot.tag = element_text(face = "bold", size = 10),
      plot.caption = element_text(size = 6.8, hjust = 0)
    )
  )

base <- file.path(figure_dir, "Figure4_immune_classifier_robustness")
embedded_pdf <- function(filename, width, height, bg, ...) {
  stopifnot(capabilities("aqua"))
  grDevices::quartz(type="pdf", file=filename, width=width, height=height, bg=bg)
}
ggsave(paste0(base, ".pdf"), figure4, width = 183, height = 190, units = "mm",
       device = embedded_pdf, bg = "white")
ggsave(paste0(base, ".png"), figure4, width = 183, height = 190, units = "mm",
       dpi = 600, bg = "white")

write.table(immune, paste0(base, "_source_data_A_immune.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
random_source <- rbind(
  data.frame(record_type = "observed_renal_derived", iteration = NA_integer_, auc = observed_auc),
  data.frame(record_type = "size_matched_random_gene", iteration = random_auc$iteration, auc = random_auc$auc)
)
write.table(random_source, paste0(base, "_source_data_B_random_auc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(external, paste0(base, "_source_data_C_external.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
source_index <- data.frame(
  panel = c("A", "B", "C"),
  file = c(
    "Figure4_immune_classifier_robustness_source_data_A_immune.tsv",
    "Figure4_immune_classifier_robustness_source_data_B_random_auc.tsv",
    "Figure4_immune_classifier_robustness_source_data_C_external.tsv"
  ),
  contents = c(
    "immune programme effects, approximate confidence intervals and BH results",
    "observed and 100 size-matched random-gene nested OOF AUC values",
    "all four frozen external cohorts, status, sample sizes, coverage, AUC and bootstrap intervals"
  ),
  claim_ids = c("C16", "C17", "C18-C19, C21")
)
write.table(source_index, paste0(base, "_source_data.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("FIGURE4_COMPLETE\n")
