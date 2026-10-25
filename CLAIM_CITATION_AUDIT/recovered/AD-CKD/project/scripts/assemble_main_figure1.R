#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
.libPaths(c(file.path(root, "project", "Rlib-4.5"),
            file.path(root, "renv", "library", "macos", "R-4.5", "aarch64-apple-darwin20"),
            .libPaths()))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
geo <- file.path(root, "project", "results", "geo")
fig_dir <- file.path(root, "project", "results", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_classic(base_size = 8.2) + theme(
  plot.title = element_text(face = "bold", size = 9.2), plot.subtitle = element_text(size = 7.1),
  axis.title = element_text(face = "bold"), plot.margin = margin(5, 7, 5, 5)
)

# A: frozen animal contrasts and the exposure-matched sensitivity comparison.
design <- data.frame(
  row = factor(c("Renal dysfunction", "Renal dysfunction", "Experimental AD", "Experimental AD", "BAPN sensitivity", "BAPN sensitivity"),
               levels = rev(c("Renal dysfunction", "Experimental AD", "BAPN sensitivity"))),
  x = rep(c(1, 2), 3),
  group = c("Uremic", "Untreated", "BAPN + AD", "Untreated", "BAPN + AD", "BAPN, no AD"),
  n = c(5, 5, 16, 20, 16, 12),
  fill = c("#56B4E9", "#E8E8E8", "#E69F00", "#E8E8E8", "#E69F00", "#F0E442")
)
design$label <- sprintf("%s\nn=%d", design$group, design$n)
p_a <- ggplot(design, aes(x, row)) +
  geom_tile(aes(fill = interaction(row, x)), width = 0.82, height = 0.66, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = label), size = 2.15, lineheight = 0.95) +
  scale_fill_manual(values = setNames(design$fill, interaction(design$row, design$x)), guide = "none") +
  scale_x_continuous(limits = c(0.5, 2.5), breaks = NULL) +
  labs(title = "Animal aortic transcriptomic contrasts",
       subtitle = "Primary comparisons and an exposure-matched sensitivity analysis", x = NULL, y = NULL) +
  theme_void(base_size = 8.2) +
  theme(plot.title = element_text(face = "bold", size = 9.2), plot.subtitle = element_text(size = 7.1),
        axis.text.y = element_text(size = 7.3), plot.margin = margin(5, 7, 5, 5))

volcano_plot <- function(tab, title, subtitle) {
  tab$status <- "Not DEG"
  tab$status[tab$adj.P.Val < 0.05 & tab$logFC > 0.5] <- "Up"
  tab$status[tab$adj.P.Val < 0.05 & tab$logFC < -0.5] <- "Down"
  tab$neglog10q <- -log10(pmax(tab$adj.P.Val, .Machine$double.xmin))
  tab$status <- factor(tab$status, levels = c("Down", "Not DEG", "Up"))
  p <- ggplot(tab, aes(logFC, neglog10q, colour = status)) +
    geom_point(size = 0.55, alpha = 0.62) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey55") +
    scale_colour_manual(values = c("Down" = "#0072B2", "Not DEG" = "grey78", "Up" = "#D55E00")) +
    labs(title = title, subtitle = subtitle, colour = NULL,
         x = expression(log[2]("fold change")), y = expression(-log[10]("BH q"))) +
    theme_pub + theme(legend.position = "bottom")
  list(plot = p, data = tab)
}
rat <- read.csv(file.path(geo, "GSE146638_all_genes.csv"), check.names = FALSE)
mouse <- read.csv(file.path(geo, "GSE325608_all_genes.csv"), check.names = FALSE)
v_b <- volcano_plot(rat, "Renal-dysfunction contrast", "Uremic versus untreated; 5 versus 5")
v_c <- volcano_plot(mouse, "Experimental AD contrast", "BAPN + AD versus untreated; 16 versus 20; sex-adjusted")
p_b <- v_b$plot; p_c <- v_c$plot

# D: direction-aware human 1:1 orthologue overlap.
rat_map <- read.csv(file.path(geo, "ortholog_rat_human_1to1.csv"), check.names = FALSE)
mouse_map <- read.csv(file.path(geo, "ortholog_mouse_human_1to1.csv"), check.names = FALSE)
rat_deg_tab <- read.csv(file.path(geo, "GSE146638_DEG.csv"), check.names = FALSE)
mouse_deg_tab <- read.csv(file.path(geo, "GSE325608_DEG.csv"), check.names = FALSE)
rat_deg <- merge(subset(rat_deg_tab, adj.P.Val < 0.05 & abs(logFC) > 0.5), rat_map,
                 by.x = "gene_symbol", by.y = "source_symbol")
mouse_deg <- merge(subset(mouse_deg_tab, adj.P.Val < 0.05 & abs(logFC) > 0.5), mouse_map,
                   by.x = "gene_symbol", by.y = "source_symbol")
collapse_lfc <- function(x) x[which.max(abs(x))]
rat_h <- aggregate(logFC ~ human_symbol, rat_deg, collapse_lfc)
mouse_h <- aggregate(logFC ~ human_symbol, mouse_deg, collapse_lfc)
overlap <- merge(rat_h, mouse_h, by = "human_symbol", suffixes = c("_renal", "_ad"))
overlap$direction <- ifelse(overlap$logFC_renal > 0 & overlap$logFC_ad > 0, "Up in both",
                     ifelse(overlap$logFC_renal < 0 & overlap$logFC_ad < 0, "Down in both",
                     ifelse(overlap$logFC_renal > 0, "Up renal / down AD", "Down renal / up AD")))
direction_order <- c("Up in both", "Down in both", "Up renal / down AD", "Down renal / up AD")
direction <- as.data.frame(table(factor(overlap$direction, levels = direction_order)))
names(direction) <- c("direction", "genes")
direction$concordance <- ifelse(direction$direction %in% c("Up in both", "Down in both"), "Concordant", "Discordant")
sens <- read.csv(file.path(geo, "GSE325608_AAD_vs_NonAD_sensitivity.csv"), check.names = FALSE)
sens_n <- sum(sens$adj.P.Val < 0.05 & abs(sens$logFC) > 0.5)
p_d <- ggplot(direction, aes(genes, factor(direction, levels = rev(direction_order)), fill = concordance)) +
  geom_col(width = 0.68) + geom_text(aes(label = genes), hjust = -0.18, size = 2.25) +
  scale_fill_manual(values = c("Concordant" = "#009E73", "Discordant" = "#CC79A7")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = "Direction of the 549-gene overlap",
       subtitle = sprintf("478 concordant and 71 discordant\nAAD versus BAPN-nonAD sensitivity: %s DEGs", format(sens_n, big.mark = ",")),
       x = "Human 1:1 orthologue genes", y = NULL, fill = NULL) + theme_pub +
  theme(legend.position = "bottom")

fig <- (p_a | p_d) / (p_b | p_c) +
  plot_layout(widths = c(0.9, 1.1), heights = c(0.72, 1.28)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 10)))
base <- file.path(fig_dir, "Figure1_animal_transcriptomic_discovery")
ggsave(paste0(base, ".pdf"), fig, width = 183, height = 160, units = "mm", device = "pdf", bg = "white", useDingbats = FALSE)
ggsave(paste0(base, ".png"), fig, width = 183, height = 160, units = "mm", dpi = 600, bg = "white")
write.table(design[, c("row", "group", "n")], paste0(base, "_source_data_A_design.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(v_b$data[, c("gene_symbol", "logFC", "adj.P.Val", "status")], paste0(base, "_source_data_B_renal_volcano.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(v_c$data[, c("gene_symbol", "logFC", "adj.P.Val", "status")], paste0(base, "_source_data_C_AD_volcano.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
direction$sensitivity_comparison <- "AAD versus BAPN-exposed non-AD"
direction$sensitivity_deg_count <- sens_n
write.table(direction, paste0(base, "_source_data_D_direction_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(overlap[, c("human_symbol", "logFC_renal", "logFC_ad", "direction")], paste0(base, "_source_data_D_overlap.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
index <- data.frame(panel = LETTERS[1:4], file = c(
  paste0(basename(base), "_source_data_A_design.tsv"),
  paste0(basename(base), "_source_data_D_direction_summary.tsv; ", basename(base), "_source_data_D_overlap.tsv"),
  paste0(basename(base), "_source_data_B_renal_volcano.tsv"),
  paste0(basename(base), "_source_data_C_AD_volcano.tsv")
), claim_ids = c("C1", "C1", "C1", "C1"))
write.table(index, paste0(base, "_source_data.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("FIGURE1_COMPLETE\n")
