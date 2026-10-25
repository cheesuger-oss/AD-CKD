#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
script_path <- sub("^--file=", "", script_arg)
root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
setwd(root)
.libPaths(c("renv/library/macos/R-4.5/aarch64-apple-darwin20", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  library(pheatmap)
})

out <- "project/results/geo"
fig <- file.path(out, "figures")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
set.seed(20260807)

read_meta <- function(path) {
  x <- fread(path)
  x$expression_column <- as.character(x$expression_column)
  x
}

collapse_sum <- function(mat, ids) {
  stopifnot(nrow(mat) == length(ids))
  rowsum(mat, group = ids, reorder = FALSE, na.rm = TRUE)
}

plot_pca <- function(log_expr, meta, title, filename) {
  vars <- apply(log_expr, 1, var)
  keep <- head(order(vars, decreasing = TRUE), min(1000, length(vars)))
  pc <- prcomp(t(log_expr[keep, , drop = FALSE]), center = TRUE, scale. = FALSE)
  pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  dat <- cbind(meta, as.data.frame(pc$x[, 1:2, drop = FALSE]))
  p <- ggplot(dat, aes(PC1, PC2, color = primary_group, shape = sex)) +
    geom_point(size = 3, alpha = 0.9) +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2]), color = "Group", shape = "Sex") +
    theme_classic(base_size = 12)
  ggsave(filename, p, width = 7, height = 5, dpi = 300)
  fwrite(dat, sub("\\.png$", ".csv", filename))
  invisible(dat)
}

plot_volcano <- function(tab, title, filename) {
  tab$status <- "Not significant"
  tab$status[tab$adj.P.Val < 0.05 & tab$logFC > 0.5] <- "Up"
  tab$status[tab$adj.P.Val < 0.05 & tab$logFC < -0.5] <- "Down"
  tab$neglog10 <- -log10(pmax(tab$adj.P.Val, .Machine$double.xmin))
  p <- ggplot(tab, aes(logFC, neglog10, color = status)) +
    geom_point(alpha = 0.55, size = 1) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50") +
    scale_color_manual(values = c(Down = "#2878B5", `Not significant` = "grey75", Up = "#C82423")) +
    labs(title = title, x = "log2 fold change", y = "-log10(FDR)", color = NULL) +
    theme_classic(base_size = 12)
  ggsave(filename, p, width = 7, height = 5, dpi = 300)
}

sample_qc <- function(raw, log_expr, meta, dataset) {
  lib <- colSums(raw)
  detected <- colSums(raw > 0)
  med <- apply(log_expr, 2, median)
  q <- data.frame(dataset = dataset, expression_column = colnames(raw), library_sum = lib,
                  detected_genes = detected, median_log_expression = med)
  q <- merge(q, meta, by = c("dataset", "expression_column"), all.x = TRUE, sort = FALSE)
  q
}

# GSE146638: processed continuous expression, primary Control vs Uremic.
m1 <- read_meta(file.path(out, "GSE146638_sample_metadata.tsv"))
x1 <- fread("project/data/geo/processed/GSE146638_expression.tsv.gz")
ids1 <- x1[[1]]
mat1 <- as.matrix(x1[, -1])
storage.mode(mat1) <- "double"
rownames(mat1) <- ids1
mat1 <- collapse_sum(mat1, ids1)
sel1 <- m1$include_primary == "yes"
m1p <- m1[sel1, ]
mat1p <- mat1[, m1p$expression_column, drop = FALSE]
keep1 <- rowMeans(mat1p) > 1
log1 <- log2(mat1p[keep1, , drop = FALSE] + 1)
m1p$primary_group <- factor(m1p$primary_group, levels = c("Control", "Uremic"))
design1 <- model.matrix(~ primary_group, data = m1p)
fit1 <- eBayes(lmFit(log1, design1), trend = TRUE)
r1 <- topTable(fit1, coef = "primary_groupUremic", number = Inf, sort.by = "P")
r1$gene_symbol <- rownames(r1)
r1 <- r1[, c("gene_symbol", setdiff(colnames(r1), "gene_symbol"))]
fwrite(r1, file.path(out, "GSE146638_all_genes.csv"))
d1 <- r1[r1$adj.P.Val < 0.05 & abs(r1$logFC) > 0.5, ]
fwrite(d1, file.path(out, "GSE146638_DEG.csv"))
plot_pca(log1, m1p, "GSE146638: Control vs Uremic", file.path(fig, "GSE146638_PCA.png"))
plot_volcano(r1, "GSE146638: Uremic vs Control", file.path(fig, "GSE146638_volcano.png"))
pheatmap(cor(log1, method = "spearman"), annotation_col = data.frame(Group = m1p$primary_group, row.names = m1p$expression_column),
         filename = file.path(fig, "GSE146638_sample_correlation.png"), width = 8, height = 7)

# GSE325608: RSEM estimated counts, main AAD vs Control adjusted for sex.
m2 <- read_meta(file.path(out, "GSE325608_sample_metadata.tsv"))
x2 <- fread("project/data/geo/GSE325608_raw_counts.txt.gz")
ids2 <- x2[[1]]
mat2 <- as.matrix(x2[, -1])
storage.mode(mat2) <- "double"
rownames(mat2) <- ids2
sel2 <- m2$include_primary == "yes"
m2p <- m2[sel2, ]
mat2p <- mat2[, m2p$expression_column, drop = FALSE]
keep2 <- rowSums(mat2p >= 10) >= 5
mat2f <- mat2p[keep2, , drop = FALSE]
m2p$sex <- factor(m2p$sex, levels = c("Male", "Female"))
m2p$primary_group <- factor(m2p$primary_group, levels = c("Control", "AAD"))
design2 <- model.matrix(~ sex + primary_group, data = m2p)
v2 <- voom(mat2f, design2, plot = TRUE, save.plot = TRUE)
png(file.path(fig, "GSE325608_voom_mean_variance.png"), width = 1800, height = 1400, res = 250)
plot(v2$voom.xy, xlab = "Mean log2 expression", ylab = "Sqrt residual SD", main = "GSE325608 voom mean-variance trend", pch = 16, cex = 0.5)
lines(v2$voom.line, col = "red", lwd = 2)
dev.off()
fit2 <- eBayes(lmFit(v2, design2), trend = FALSE)
r2 <- topTable(fit2, coef = "primary_groupAAD", number = Inf, sort.by = "P")
r2$gene_symbol <- rownames(r2)
r2 <- r2[, c("gene_symbol", setdiff(colnames(r2), "gene_symbol"))]
fwrite(r2, file.path(out, "GSE325608_all_genes.csv"))
d2 <- r2[r2$adj.P.Val < 0.05 & abs(r2$logFC) > 0.5, ]
fwrite(d2, file.path(out, "GSE325608_DEG.csv"))
plot_pca(v2$E, m2p, "GSE325608: Control vs AAD", file.path(fig, "GSE325608_PCA.png"))
plot_volcano(r2, "GSE325608: AAD vs Control (sex-adjusted)", file.path(fig, "GSE325608_volcano.png"))
pheatmap(cor(v2$E, method = "spearman"), annotation_col = data.frame(Group = m2p$primary_group, Sex = m2p$sex, row.names = m2p$expression_column),
         filename = file.path(fig, "GSE325608_sample_correlation.png"), width = 11, height = 10)

# Sensitivity: among BAPN-exposed survivors, AAD vs NonAD adjusted for sex.
s2 <- m2$pathology %in% c("NonAD", "AAD")
m2s <- m2[s2, ]
mat2s <- mat2[, m2s$expression_column, drop = FALSE]
keep2s <- rowSums(mat2s >= 10) >= 5
m2s$sex <- factor(m2s$sex, levels = c("Male", "Female"))
m2s$pathology <- factor(m2s$pathology, levels = c("NonAD", "AAD"))
design2s <- model.matrix(~ sex + pathology, data = m2s)
v2s <- voom(mat2s[keep2s, , drop = FALSE], design2s, plot = FALSE)
fit2s <- eBayes(lmFit(v2s, design2s))
r2s <- topTable(fit2s, coef = "pathologyAAD", number = Inf, sort.by = "P")
r2s$gene_symbol <- rownames(r2s)
r2s <- r2s[, c("gene_symbol", setdiff(colnames(r2s), "gene_symbol"))]
fwrite(r2s, file.path(out, "GSE325608_AAD_vs_NonAD_sensitivity.csv"))

qc <- rbindlist(list(
  sample_qc(mat1p, log1, m1p, "GSE146638"),
  sample_qc(mat2p, v2$E, m2p, "GSE325608")
), fill = TRUE, use.names = TRUE)
fwrite(qc, file.path(out, "phase1_sample_qc_metrics.csv"))

summary <- data.frame(
  dataset = c("GSE146638", "GSE325608"),
  method = c("limma on log2(processed expression + 1)", "limma-voom on RSEM estimated counts"),
  design = c("~ primary_group", "~ sex + primary_group"),
  control_n = c(sum(m1p$primary_group == "Control"), sum(m2p$primary_group == "Control")),
  disease_n = c(sum(m1p$primary_group == "Uremic"), sum(m2p$primary_group == "AAD")),
  tested_genes = c(nrow(r1), nrow(r2)),
  deg_up = c(sum(d1$logFC > 0), sum(d2$logFC > 0)),
  deg_down = c(sum(d1$logFC < 0), sum(d2$logFC < 0)),
  stringsAsFactors = FALSE
)
fwrite(summary, file.path(out, "phase1_differential_summary.csv"))

writeLines(capture.output(sessionInfo()), file.path(out, "phase1_R_sessionInfo.txt"))
print(summary)
