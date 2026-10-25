#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
.libPaths(c(file.path(root, "project", "Rlib-4.5"), file.path(root, "renv", "library", "macos", "R-4.5", "aarch64-apple-darwin20"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(limma); library(ggplot2); library(patchwork) })
geo <- file.path(root, "project", "results", "geo"); out <- file.path(root, "project", "results", "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
theme_pub <- theme_classic(base_size = 7.6) + theme(plot.title = element_text(face = "bold", size = 8.7), plot.subtitle = element_text(size = 6.8), axis.title = element_text(face = "bold"), plot.margin = margin(4,5,4,4))
zscore <- function(x) { z <- t(scale(t(x))); z[!is.finite(z)] <- 0; pmax(pmin(z, 2.5), -2.5) }

# Reconstruct the primary normalized matrices.
m1 <- fread(file.path(geo, "GSE146638_sample_metadata.tsv")); x1 <- fread(file.path(root, "project", "data", "geo", "processed", "GSE146638_expression.tsv.gz"))
ids1 <- x1[[1]]; mat1 <- rowsum(as.matrix(x1[, -1]), group = ids1, reorder = FALSE, na.rm = TRUE); m1 <- m1[include_primary == "yes"]
m1$primary_group <- factor(m1$primary_group, levels = c("Control", "Uremic"))
log1 <- log2(mat1[, m1$expression_column, drop = FALSE] + 1); d1 <- fread(file.path(geo, "GSE146638_DEG.csv")); top1 <- head(d1[order(-abs(t))]$gene_symbol, 30); h1 <- zscore(log1[top1, , drop = FALSE])
m2all <- fread(file.path(geo, "GSE325608_sample_metadata.tsv")); x2 <- fread(file.path(root, "project", "data", "geo", "GSE325608_raw_counts.txt.gz"))
ids2 <- x2[[1]]; mat2 <- as.matrix(x2[, -1]); storage.mode(mat2) <- "double"; rownames(mat2) <- ids2
m2 <- m2all[include_primary == "yes"]; keep2 <- rowSums(mat2[, m2$expression_column, drop = FALSE] >= 10) >= 5
m2$sex <- factor(m2$sex, levels = c("Male", "Female")); m2$primary_group <- factor(m2$primary_group, levels = c("Control", "AAD"))
v2 <- voom(mat2[keep2, m2$expression_column, drop = FALSE], model.matrix(~ sex + primary_group, data = m2), plot = FALSE)
d2 <- fread(file.path(geo, "GSE325608_DEG.csv")); top2 <- head(d2[order(-abs(t))]$gene_symbol, 30); h2 <- zscore(v2$E[top2, , drop = FALSE])

heat_data <- function(mat, meta, dataset) {
  z <- as.data.frame(as.table(mat), stringsAsFactors = FALSE); names(z) <- c("gene", "sample", "row_z")
  z$dataset <- dataset
  z$group <- factor(setNames(as.character(meta$primary_group), meta$expression_column)[z$sample], levels = levels(meta$primary_group))
  z$gene <- factor(z$gene, levels = rev(rownames(mat))); z$sample <- factor(z$sample, levels = meta$expression_column); z
}
heat_plot <- function(dat, title) {
  ggplot(dat, aes(sample, gene, fill = row_z)) + geom_tile() +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish, name = "Row z-score") +
    facet_grid(~group, scales = "free_x", space = "free_x") +
    labs(title = title, subtitle = "Top 30 DEGs by absolute moderated t statistic", x = "Samples", y = NULL) + theme_pub +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "bold"), panel.spacing.x = grid::unit(1, "mm"))
}
hd1 <- heat_data(h1, m1, "GSE146638"); hd2 <- heat_data(h2, m2, "GSE325608")
p_a <- heat_plot(hd1, "Renal-dysfunction top-DEG heatmap")
p_b <- heat_plot(hd2, "Experimental AD top-DEG heatmap")

# Exposure-matched AAD versus non-AD sensitivity result.
sens <- read.csv(file.path(geo, "GSE325608_AAD_vs_NonAD_sensitivity.csv"), check.names = FALSE)
sens$status <- "Not DEG"; sens$status[sens$adj.P.Val < 0.05 & sens$logFC > 0.5] <- "Up"; sens$status[sens$adj.P.Val < 0.05 & sens$logFC < -0.5] <- "Down"
sens$neglog10q <- -log10(pmax(sens$adj.P.Val, .Machine$double.xmin)); sens$status <- factor(sens$status, levels = c("Down", "Not DEG", "Up"))
p_c <- ggplot(sens, aes(logFC, neglog10q, colour = status)) + geom_point(size = 0.48, alpha = 0.6) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey55") + geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey55") +
  scale_colour_manual(values = c("Down" = "#0072B2", "Not DEG" = "grey78", "Up" = "#D55E00")) +
  labs(title = "Exposure-matched sensitivity contrast", subtitle = "BAPN + AD versus BAPN without AD; 16 versus 12; sex-adjusted; 1,963 DEGs",
       x = expression(log[2]("fold change")), y = expression(-log[10]("BH q")), colour = NULL) + theme_pub + theme(legend.position = "bottom")

primary <- read.csv(file.path(geo, "GSE325608_all_genes.csv"), check.names = FALSE)
compare <- merge(primary[, c("gene_symbol", "logFC")], sens[, c("gene_symbol", "logFC", "adj.P.Val")], by = "gene_symbol", suffixes = c("_primary", "_sensitivity"))
compare$sensitivity_deg <- compare$adj.P.Val < 0.05 & abs(compare$logFC_sensitivity) > 0.5
rho <- cor(compare$logFC_primary, compare$logFC_sensitivity, method = "spearman", use = "complete.obs")
p_d <- ggplot(compare, aes(logFC_primary, logFC_sensitivity, colour = sensitivity_deg)) + geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(size = 0.5, alpha = 0.45) + scale_colour_manual(values = c(`FALSE` = "grey72", `TRUE` = "#CC79A7"), labels = c(`FALSE` = "Not sensitivity DEG", `TRUE` = "Sensitivity DEG")) +
  labs(title = "Primary and sensitivity effects were related", subtitle = sprintf("Spearman rho = %.3f; effects were not identical", rho),
       x = "Primary AD versus untreated log2 fold change", y = "AAD versus BAPN-nonAD log2 fold change", colour = NULL) + theme_pub + theme(legend.position = "bottom")

fig <- (p_a | p_b) / (p_c | p_d) + plot_layout(heights = c(1.05, 0.95)) + plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 9.5)))
base <- file.path(out, "Supplementary_Figure_S2_animal_expression_detail")
ggsave(paste0(base, ".pdf"), fig, width = 183, height = 190, units = "mm", device = "pdf", bg = "white", useDingbats = FALSE)
ggsave(paste0(base, ".png"), fig, width = 183, height = 190, units = "mm", dpi = 600, bg = "white")
write.table(hd1, paste0(base, "_source_data_A_heatmap.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(hd2, paste0(base, "_source_data_B_heatmap.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(sens[, c("gene_symbol", "logFC", "adj.P.Val", "status")], paste0(base, "_source_data_C_sensitivity_volcano.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(compare, paste0(base, "_source_data_D_effect_comparison.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
index <- data.frame(panel = LETTERS[1:4], file = paste0(basename(base), c("_source_data_A_heatmap.tsv", "_source_data_B_heatmap.tsv", "_source_data_C_sensitivity_volcano.tsv", "_source_data_D_effect_comparison.tsv")), claim_ids = c("C1 context", "C1 context", "C1 context", "C1 context"))
write.table(index, paste0(base, "_source_data.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("SUPPLEMENTARY_FIGURE_S2_COMPLETE\n")
