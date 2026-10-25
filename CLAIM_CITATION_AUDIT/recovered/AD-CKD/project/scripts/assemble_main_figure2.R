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

ids <- c("GO:0030198", "GO:0001503", "GO:0006939", "GO:0048659", "GO:0048660", "GO:0048661")
labels <- c(
  "GO:0030198" = "Extracellular matrix\norganization",
  "GO:0001503" = "Ossification",
  "GO:0006939" = "Smooth muscle\ncontraction",
  "GO:0048659" = "Smooth muscle cell\nproliferation",
  "GO:0048660" = "Regulation of smooth\nmuscle proliferation",
  "GO:0048661" = "Positive regulation of\nsmooth muscle proliferation"
)
theme_pub <- theme_classic(base_size = 8.2) + theme(
  plot.title = element_text(face = "bold", size = 9.2), plot.subtitle = element_text(size = 7.1),
  axis.title = element_text(face = "bold"), legend.title = element_blank(),
  plot.margin = margin(5, 7, 5, 5)
)

# A: thresholded over-representation analysis used for the shared-pathway screen.
conv <- read.csv(file.path(geo, "convergence_pathways.csv"), check.names = FALSE)
conv <- conv[conv$ID %in% ids & conv$focus_match, ]
ora <- rbind(
  data.frame(pathway_id = conv$ID, model = "Renal dysfunction", q = conv$FDR_UREMIC),
  data.frame(pathway_id = conv$ID, model = "Aortic dissection", q = conv$FDR_AD)
)
ora$pathway <- factor(labels[ora$pathway_id], levels = rev(labels[ids]))
p_a <- ggplot(ora, aes(-log10(q), pathway, colour = model, shape = model)) +
  geom_vline(xintercept = -log10(0.01), linetype = "dashed", colour = "grey55") +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c("Renal dysfunction" = "#0072B2", "Aortic dissection" = "#D55E00")) +
  scale_shape_manual(values = c("Renal dysfunction" = 16, "Aortic dissection" = 17)) +
  labs(title = "Selected programmes overlapped in the animal models",
       subtitle = "DEG over-representation; dashed line, BH q=0.01",
       x = expression(-log[10]("BH q")), y = NULL) + theme_pub +
  theme(legend.position = "bottom")

# B: all-ranked GSEA is a separate sensitivity layer and need not mirror ORA.
gsea <- read.csv(file.path(geo, "GSEA_all_results.csv"), check.names = FALSE)
gsea <- gsea[gsea$ID %in% ids & gsea$dataset %in% c("Uremic", "AAD"), ]
gsea$model <- factor(ifelse(gsea$dataset == "Uremic", "Renal dysfunction", "Aortic dissection"),
                     levels = c("Renal dysfunction", "Aortic dissection"))
gsea$pathway_id <- gsea$ID
gsea$pathway <- factor(labels[gsea$pathway_id], levels = rev(labels[ids]))
gsea$passed <- gsea$p.adjust < 0.05
p_b <- ggplot(gsea, aes(NES, pathway, colour = model, shape = passed)) +
  geom_vline(xintercept = 0, colour = "grey60") + geom_point(size = 2.4) +
  scale_colour_manual(values = c("Renal dysfunction" = "#0072B2", "Aortic dissection" = "#D55E00")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "BH q<0.05", `FALSE` = "BH q>=0.05")) +
  labs(title = "All-ranked evidence\nwas more selective",
       subtitle = "Filled points passed BH correction",
       x = "Normalized enrichment score", y = NULL) + theme_pub +
  guides(colour = "none", shape = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.3), legend.key.width = grid::unit(3.5, "mm"))

# C: prespecified six-pathway test in human vascular tissue.
human <- read.csv(file.path(geo, "human_validation", "human_validation_six_pathway_GSEA.csv"), check.names = FALSE)
human$pathway_id <- human$GO_ID
human$pathway <- factor(labels[human$pathway_id], levels = rev(labels[ids]))
human$dataset_label <- factor(
  c("GSE153434" = "Human AD aorta",
    "GSE38750" = "Uremic vascular tissue\n(cross-site)",
    "GSE38751" = "Uremic internal mammary artery\n(same-site)")[human$dataset],
  levels = c("Human AD aorta", "Uremic vascular tissue\n(cross-site)",
             "Uremic internal mammary artery\n(same-site)"))
human$passed <- human$BH_q < 0.05
p_c <- ggplot(human, aes(NES, pathway, colour = dataset_label, shape = passed)) +
  geom_vline(xintercept = 0, colour = "grey60") + geom_point(size = 2.4) +
  scale_colour_manual(values = c("Human AD aorta" = "#D55E00",
                                 "Uremic vascular tissue\n(cross-site)" = "#0072B2",
                                 "Uremic internal mammary artery\n(same-site)" = "#009E73")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "BH q<0.05", `FALSE` = "BH q>=0.05")) +
  labs(title = "Human support varied across\narterial datasets",
       subtitle = "Six prespecified programmes tested within each dataset",
       x = "Normalized enrichment score", y = NULL) + theme_pub +
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.3), legend.key.width = grid::unit(3.5, "mm"))

# D: design constraints that qualify interpretation of both renal vascular datasets.
design <- data.frame(
  dataset = factor(c("GSE153434", "GSE153434", "GSE38750", "GSE38750", "GSE38751", "GSE38751"),
                   levels = c("GSE38751", "GSE38750", "GSE153434")),
  group = c("AAD", "Control", "Uremia", "Control", "Uremia", "Control"),
  tissue = c("Ascending aorta", "Ascending aorta", "Iliac artery", "Renal artery",
             "Internal mammary artery", "Internal mammary artery"),
  n = c(10, 10, 15, 19, 5, 5), x = rep(c(1, 2), 3)
)
design$label <- sprintf("%s\n%s\nn=%d", design$group, design$tissue, design$n)
p_d <- ggplot(design, aes(x, dataset)) +
  geom_tile(aes(fill = dataset), width = 0.84, height = 0.62, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = label), size = 2.05, lineheight = 0.93) +
  scale_fill_manual(values = c("GSE153434" = "#E8E8E8", "GSE38750" = "#F0E442",
                               "GSE38751" = "#CCEDDF"), guide = "none") +
  scale_x_continuous(limits = c(0.45, 2.55), breaks = NULL) +
  labs(title = "Sampling design limits attribution",
       subtitle = "GSE38751 matches arterial site but is small and\ncomes from the same study as GSE38750", x = NULL, y = NULL) +
  theme_void(base_size = 8.2) +
  theme(plot.title = element_text(face = "bold", size = 9.2), plot.subtitle = element_text(size = 7.1),
        axis.text.y = element_text(size = 7.2), plot.margin = margin(5, 7, 5, 5))

fig <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(widths = c(1.05, 0.95), heights = c(1, 1)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 10)))

base <- file.path(fig_dir, "Figure2_pathway_convergence_human_support")
ggsave(paste0(base, ".pdf"), fig, width = 183, height = 175, units = "mm", device = "pdf", bg = "white", useDingbats = FALSE)
ggsave(paste0(base, ".png"), fig, width = 183, height = 175, units = "mm", dpi = 600, bg = "white")
ora$pathway_name <- gsub("\n", " ", as.character(ora$pathway))
write.table(ora[, c("pathway_id", "pathway_name", "model", "q")],
            paste0(base, "_source_data_A_ORA.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(gsea[, c("pathway_id", "Description", "model", "NES", "pvalue", "p.adjust", "passed")],
            paste0(base, "_source_data_B_animal_GSEA.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
human$pathway_name <- gsub("\n", " ", as.character(human$pathway))
write.table(human[, c("dataset", "pathway_id", "pathway_name", "set_size", "ES", "NES", "P_value", "BH_q", "passed")],
            paste0(base, "_source_data_C_human_GSEA.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(design[, c("dataset", "group", "tissue", "n")],
            paste0(base, "_source_data_D_design.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
index <- data.frame(panel = LETTERS[1:4], file = paste0(basename(base), c("_source_data_A_ORA.tsv", "_source_data_B_animal_GSEA.tsv", "_source_data_C_human_GSEA.tsv", "_source_data_D_design.tsv")), claim_ids = c("C2", "C3", "C10-C11", "C10-C11"))
write.table(index, paste0(base, "_source_data.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("FIGURE2_COMPLETE\n")
