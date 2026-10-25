#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
.libPaths(c(file.path(root, "project", "Rlib-4.5"),
            file.path(root, "renv", "library", "macos", "R-4.5", "aarch64-apple-darwin20"),
            .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(limma); library(ggplot2); library(patchwork) })

geo <- file.path(root, "project", "results", "geo")
out <- file.path(root, "project", "results", "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20260807)

theme_qc <- theme_classic(base_size = 7.5) + theme(
  plot.title = element_text(face = "bold", size = 8.5), plot.subtitle = element_text(size = 6.8),
  axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"),
  plot.margin = margin(4, 5, 4, 4)
)
read_meta <- function(path) { x <- fread(path); x$expression_column <- as.character(x$expression_column); x }
collapse_sum <- function(mat, ids) rowsum(mat, group = ids, reorder = FALSE, na.rm = TRUE)

m1 <- read_meta(file.path(geo, "GSE146638_sample_metadata.tsv"))
x1 <- fread(file.path(root, "project", "data", "geo", "processed", "GSE146638_expression.tsv.gz"))
ids1 <- x1[[1]]; mat1 <- as.matrix(x1[, -1]); storage.mode(mat1) <- "double"; rownames(mat1) <- ids1
mat1 <- collapse_sum(mat1, ids1); m1 <- m1[m1$include_primary == "yes", ]
log1 <- log2(mat1[rowMeans(mat1[, m1$expression_column, drop = FALSE]) > 1, m1$expression_column, drop = FALSE] + 1)
m1$primary_group <- factor(m1$primary_group, levels = c("Control", "Uremic"))

m2 <- read_meta(file.path(geo, "GSE325608_sample_metadata.tsv"))
x2 <- fread(file.path(root, "project", "data", "geo", "GSE325608_raw_counts.txt.gz"))
ids2 <- x2[[1]]; mat2 <- as.matrix(x2[, -1]); storage.mode(mat2) <- "double"; rownames(mat2) <- ids2
m2 <- m2[m2$include_primary == "yes", ]; mat2 <- mat2[, m2$expression_column, drop = FALSE]
mat2 <- mat2[rowSums(mat2 >= 10) >= 5, , drop = FALSE]
m2$sex <- factor(m2$sex, levels = c("Male", "Female")); m2$primary_group <- factor(m2$primary_group, levels = c("Control", "AAD"))
design2 <- model.matrix(~ sex + primary_group, data = m2); v2 <- voom(mat2, design2, plot = TRUE, save.plot = TRUE)

pca_data <- function(expr, meta, dataset) {
  vars <- apply(expr, 1, var); keep <- head(order(vars, decreasing = TRUE), min(1000, length(vars)))
  pc <- prcomp(t(expr[keep, , drop = FALSE]), center = TRUE, scale. = FALSE)
  pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  dat <- cbind(as.data.frame(meta, check.names = FALSE), as.data.frame(pc$x[, 1:2, drop = FALSE]))
  dat$dataset <- dataset
  attr(dat, "pct") <- pct[1:2]; dat
}
pc1 <- pca_data(log1, m1, "GSE146638"); pc2 <- pca_data(v2$E, m2, "GSE325608")
plot_pca <- function(dat, title) {
  pct <- attr(dat, "pct")
  ggplot(dat, aes(PC1, PC2, colour = primary_group, shape = sex)) + geom_point(size = 2.2, alpha = 0.9) +
    labs(title = title, subtitle = "PCA of the 1,000 most variable analysed genes",
         x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2]), colour = "Group", shape = "Sex") +
    scale_colour_manual(values = c("Control" = "grey40", "Uremic" = "#0072B2", "AAD" = "#D55E00")) +
    scale_shape_discrete(labels = function(x) ifelse(x == "not_reported", "Not reported", x)) + theme_qc
}
p_a <- plot_pca(pc1, "Renal-dysfunction model PCA")
p_b <- plot_pca(pc2, "Experimental AD model PCA")

corr_long <- function(expr, meta, dataset) {
  ord <- order(meta$primary_group, if ("sex" %in% names(meta)) meta$sex else "", meta$expression_column)
  meta <- meta[ord, ]
  z <- cor(expr, method = "spearman"); d <- as.data.frame(as.table(z), stringsAsFactors = FALSE)
  names(d) <- c("sample_x", "sample_y", "rho"); d$dataset <- dataset
  map <- setNames(as.character(meta$primary_group), meta$expression_column)
  compact <- if (dataset == "GSE325608") sub("^read_count_", "", meta$expression_column) else meta$expression_column
  display <- setNames(compact, meta$expression_column)
  group_levels <- levels(meta$primary_group)
  d$group_x <- factor(map[d$sample_x], levels = group_levels)
  d$group_y <- factor(map[d$sample_y], levels = rev(group_levels))
  d$display_x <- factor(display[d$sample_x], levels = compact)
  d$display_y <- factor(display[d$sample_y], levels = rev(compact))
  d
}
cor1 <- corr_long(log1, m1, "GSE146638"); cor2 <- corr_long(v2$E, m2, "GSE325608")
plot_corr <- function(dat, title) {
  ggplot(dat, aes(display_x, display_y, fill = rho)) + geom_tile() +
    scale_fill_viridis_c(option = "C", limits = c(min(dat$rho), 1), name = "Spearman rho") +
    facet_grid(group_y ~ group_x, scales = "free", space = "free") +
    labs(title = title, subtitle = "All retained primary-analysis samples", x = NULL, y = NULL) +
    theme_minimal(base_size = 6.7) + theme(
      plot.title = element_text(face = "bold", size = 8.5), plot.subtitle = element_text(size = 6.8),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.0),
      axis.text.y = element_text(size = 5.0), strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold", size = 5.8), panel.spacing = grid::unit(0.7, "mm"),
      panel.grid = element_blank(), plot.margin = margin(4, 5, 4, 4))
}
p_c <- plot_corr(cor1, "Renal-model sample correlation")
p_d <- plot_corr(cor2, "AD-model sample correlation")

voom <- data.frame(mean_log2_expression = v2$voom.xy$x, sqrt_residual_sd = v2$voom.xy$y)
voom_line <- data.frame(mean_log2_expression = v2$voom.line$x, fitted_sqrt_residual_sd = v2$voom.line$y)
p_e <- ggplot(voom, aes(mean_log2_expression, sqrt_residual_sd)) + geom_point(size = 0.45, alpha = 0.45, colour = "grey35") +
  geom_line(data = voom_line, aes(mean_log2_expression, fitted_sqrt_residual_sd), colour = "#D55E00", linewidth = 0.7) +
  labs(title = "Voom mean–variance diagnostic for the experimental AD model",
       subtitle = "The fitted trend was used to derive observation-level precision weights",
       x = "Mean log2 expression", y = "Square-root residual standard deviation") + theme_qc

fig <- (p_a | p_b) / (p_c | p_d) / p_e + plot_layout(heights = c(0.85, 1.35, 0.75)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 9.5)))
base <- file.path(out, "Supplementary_Figure_S1_animal_QC")
ggsave(paste0(base, ".pdf"), fig, width = 183, height = 230, units = "mm", device = "pdf", bg = "white", useDingbats = FALSE)
ggsave(paste0(base, ".png"), fig, width = 183, height = 230, units = "mm", dpi = 600, bg = "white")
write.table(pc1, paste0(base, "_source_data_A_PCA.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pc2, paste0(base, "_source_data_B_PCA.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cor1, paste0(base, "_source_data_C_correlation.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cor2, paste0(base, "_source_data_D_correlation.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(voom, paste0(base, "_source_data_E_voom_points.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(voom_line, paste0(base, "_source_data_E_voom_line.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
index <- data.frame(panel = LETTERS[1:5], file = c(
  "Supplementary_Figure_S1_animal_QC_source_data_A_PCA.tsv",
  "Supplementary_Figure_S1_animal_QC_source_data_B_PCA.tsv",
  "Supplementary_Figure_S1_animal_QC_source_data_C_correlation.tsv",
  "Supplementary_Figure_S1_animal_QC_source_data_D_correlation.tsv",
  "Supplementary_Figure_S1_animal_QC_source_data_E_voom_points.tsv; Supplementary_Figure_S1_animal_QC_source_data_E_voom_line.tsv"
), role = c("QC", "QC", "QC", "QC", "QC"))
write.table(index, paste0(base, "_source_data.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("SUPPLEMENTARY_FIGURE_S1_COMPLETE\n")
