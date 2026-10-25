#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
.libPaths(c(file.path(root, "project", "Rlib-4.5"),
            file.path(root, "renv", "library", "macos", "R-4.5", "aarch64-apple-darwin20"),
            .libPaths()))
suppressPackageStartupMessages({ library(ggplot2); library(scales); library(patchwork) })

network_dir <- file.path(root, "project", "results", "network")
vp_dir <- file.path(root, "project", "results", "virtual_perturbation")
mr_dir <- file.path(root, "project", "results", "mr")
fig_dir <- file.path(root, "project", "results", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pal <- c(other="#999999", AD_only="#0072B2", URE_only="#D55E00", intersection="#009E73")
theme_pub <- theme_classic(base_size=9) + theme(
  plot.title=element_text(face="bold", size=10), axis.title=element_text(face="bold"),
  legend.title=element_text(face="bold"), plot.margin=margin(5,7,5,5),
  strip.background=element_blank(), strip.text=element_text(face="bold")
)
wrap_text <- function(x, width=42) vapply(x, function(z) paste(strwrap(z, width=width), collapse="\n"), character(1))

save_pair <- function(plot, stem, width_mm=183, height_mm=150) {
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width=width_mm, height=height_mm,
         units="mm", device="pdf", bg="white", useDingbats=FALSE)
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width=width_mm, height=height_mm,
         units="mm", dpi=600, bg="white")
}

# Figure 3: mapping, degree topology, and enrichment.
coverage <- data.frame(
  status=factor(c("Connected", "Mapped, isolated", "Unmapped"), levels=c("Connected", "Mapped, isolated", "Unmapped")),
  count=c(511,67,1), colour=c("#0072B2","#E69F00","#999999")
)
p3a <- ggplot(coverage, aes(status, count, fill=status)) +
  geom_col(width=.68, colour="black", linewidth=.25) + geom_text(aes(label=count), vjust=-.35, size=3) +
  scale_fill_manual(values=setNames(coverage$colour, coverage$status), guide="none") +
  scale_y_continuous(limits=c(0,550), expand=expansion(mult=c(0,.02))) +
  labs(title="STRING coverage", subtitle="Frozen input: n=579", x=NULL, y="Input genes (n)") + theme_pub +
  theme(axis.text.x=element_text(angle=25,hjust=1))

hub <- read.csv(file.path(network_dir, "hub_genes.csv"), check.names=FALSE)
audit <- read.delim(file.path(network_dir, "gene_source_audit.tsv"), check.names=FALSE)
names(audit) <- tolower(names(audit))
gene_col <- if ("gene" %in% names(audit)) "gene" else names(audit)[1]
hub$module <- "other"
hub$module[hub$gene %in% c("FBN1","ACTA2","COL3A1","COL5A1","EFEMP2","EMILIN1","ENPP1","LOX","MFAP5","MYLK","PLOD3","PRKG1","SMAD2","SMAD3","TGFB2","TGFB3","TGFBR1","TGFBR2","THBS2","THSD4")] <- "AD_only"
hub$module[hub$gene %in% c("ELN","COL1A1","MMP1","MMP2","MMP9","TIMP1","BGN","MGP","RUNX2","SPP1","BMP2","IBSP","NOX4","NFKB1","NFE2L2","ACE","AGTR1","AGT","TGFB1","ENG","BMP4")] <- "URE_only"
hub$module[hub$gene %in% c("COL3A1","LOX","SMAD3","TGFBR1","TGFBR2")] <- "intersection"
top_degree <- hub[order(-hub$degree, hub$gene),][1:15,]
top_degree$gene <- factor(top_degree$gene, levels=rev(top_degree$gene))
p3b <- ggplot(top_degree, aes(gene, degree, fill=module)) +
  geom_col(width=.72, colour="black", linewidth=.15) + coord_flip() +
  scale_fill_manual(values=pal, breaks=c("intersection","AD_only","URE_only","other"),
                    labels=c("AD–URE overlap","AD-labelled","Uremic mechanism","Other")) +
  scale_y_continuous(expand=expansion(mult=c(0,.04))) +
  labs(title="Degree-ranked nodes", x=NULL, y="Degree", fill="Seed annotation") + theme_pub +
  theme(legend.position="bottom", legend.text=element_text(size=7), legend.title=element_text(size=7))

enr <- read.csv(file.path(network_dir, "enrichment_table.csv"), check.names=FALSE)
enr <- enr[order(enr$p.adjust, -enr$Count),][1:12,]
enr$Description <- wrap_text(enr$Description, 42)
enr$Description <- factor(enr$Description, levels=rev(enr$Description))
enr$neglogq <- -log10(pmax(enr$p.adjust, .Machine$double.xmin))
p3c <- ggplot(enr, aes(FoldEnrichment, Description, size=Count, colour=neglogq)) +
  geom_point(alpha=.9) + scale_colour_viridis_c(option="C", direction=1) +
  scale_size_continuous(range=c(2.5,7)) +
  labs(title="GO Biological Process enrichment", x="Fold enrichment", y=NULL,
       colour=expression(-log[10](BH~adjusted~P)), size="Gene count") + theme_pub +
  theme(legend.position="right")

fig3 <- (p3a | p3b) / p3c + plot_layout(heights=c(.9,1.15), widths=c(.72,1.28)) +
  plot_annotation(tag_levels="A", caption=paste(strwrap("STRING v12.0; human functional network; required score >=0.4. Topology and enrichment are descriptive, not causal evidence.", 115), collapse="\n"),
                  theme=theme(plot.tag=element_text(face="bold",size=11), plot.caption=element_text(size=7,hjust=0)))
save_pair(fig3, "Figure3_network_enrichment", 183, 170)

# Figure 4: frozen-network robustness.
bet <- read.csv(file.path(vp_dir, "betweenness.csv"), check.names=FALSE)
top_bet <- bet[order(-bet$betweenness_normalized, bet$gene),][1:15,]
top_bet$gene <- factor(top_bet$gene, levels=rev(top_bet$gene))
p4a <- ggplot(top_bet, aes(gene, betweenness_normalized, fill=module)) +
  geom_col(width=.72, colour="black", linewidth=.15) + coord_flip() +
  scale_fill_manual(values=pal, breaks=c("intersection","AD_only","URE_only","other"),
                    labels=c("AD-URE overlap","AD-labelled","Uremic mechanism","Other")) +
  labs(title="Betweenness-ranked nodes", x=NULL,
      y="Normalized betweenness", fill="Seed annotation") + theme_pub + theme(legend.position="bottom")

ko <- read.csv(file.path(vp_dir, "ko_fragmentation.csv"), check.names=FALSE)
ko$label <- c("Baseline","COL3A1","LOX","SMAD3","TGFBR1","TGFBR2","All five")
ko$label <- factor(ko$label, levels=ko$label)
p4b <- ggplot(ko, aes(label, largest_fraction)) + geom_col(width=.68, fill="#56B4E9", colour="black", linewidth=.25) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2), labels=percent_format(accuracy=1), expand=c(0,0)) +
  labs(title="Largest component after removal", subtitle="Six components in every condition", x=NULL, y="Largest-component fraction") + theme_pub +
  theme(axis.text.x=element_text(angle=35,hjust=1))

pairs <- read.csv(file.path(vp_dir, "ure_ad_pair_distances.csv"), check.names=FALSE)
dist_tab <- as.data.frame(table(pairs$distance)); names(dist_tab) <- c("distance","count")
dist_tab$distance <- as.numeric(as.character(dist_tab$distance))
p4c <- ggplot(dist_tab, aes(factor(distance), count)) + geom_col(width=.68, fill="#E69F00", colour="black", linewidth=.25) +
  geom_text(aes(label=count), vjust=-.3, size=3) + scale_y_continuous(limits=c(0,220), expand=c(0,0)) +
  labs(title="URE-only to AD-only distances", x="Unweighted path distance", y="Gene pairs (n=315)") + theme_pub

cov <- read.csv(file.path(vp_dir, "path_coverage.csv"), check.names=FALSE)
cov$gene <- factor(cov$gene, levels=cov$gene[order(cov$percentage)])
p4d <- ggplot(cov, aes(gene, percentage)) + geom_col(width=.68, fill="#009E73", colour="black", linewidth=.25) + coord_flip() +
  geom_text(aes(label=sprintf("%.1f%%",percentage)), hjust=-.15, size=3) +
  scale_y_continuous(limits=c(0,38), expand=c(0,0)) +
  labs(title="Shortest paths via overlap gene", x=NULL, y="URE-AD pairs (%)") + theme_pub

fig4 <- (p4a | p4b) / (p4c | p4d) +
  plot_annotation(tag_levels="A", caption=paste(strwrap("In silico deletion and shortest paths use an undirected STRING association graph. They do not model biological knockout, directionality, causality, or treatment efficacy.", 120), collapse="\n"),
                  theme=theme(plot.tag=element_text(face="bold",size=11),plot.caption=element_text(size=7,hjust=0)))
save_pair(fig4, "Figure4_virtual_perturbation", 183, 160)

# Supplementary forest plot: complete renal phenotype family, all non-FDR-significant.
ext <- read.delim(file.path(mr_dir, "renal_family", "renal_family_mr_results.tsv"), check.names=FALSE)
ext <- ext[ext$method == "Inverse variance weighted",]
labels <- c(serum_cystatin_c="Serum cystatin C", egfrcys="eGFRcys", bun_ea="BUN", ckd_ea="CKD",
            microalbuminuria="Microalbuminuria", rapid3="Rapid3", ckdi25="CKDi25", uacr_ea="UACR")
forest <- data.frame(phenotype=unname(labels[ext$exposure_id]), OR=ext$OR, low=ext$OR_lci95,
                     high=ext$OR_uci95, p=ext$pval, q=0.685111596064231, stringsAsFactors=FALSE)
main <- read.csv(file.path(mr_dir, "main_mr_results.csv"), check.names=FALSE)
ivw <- main[main$method == "Inverse variance weighted",][1,]
forest <- rbind(data.frame(phenotype="eGFRcrea", OR=ivw$OR, low=ivw$OR_lci95, high=ivw$OR_uci95,
                           p=ivw$pval, q=0.685111596064231), forest)
forest$phenotype <- factor(forest$phenotype, levels=rev(forest$phenotype))
pS <- ggplot(forest, aes(OR, phenotype)) + geom_vline(xintercept=1, linetype="dashed", colour="#666666") +
  geom_errorbar(aes(xmin=low,xmax=high), width=.18, linewidth=.5, orientation="y") + geom_point(size=2.4, shape=21, fill="#0072B2") +
  scale_x_log10(breaks=c(.05,.1,.25,.5,1,2,5,10), limits=c(.03,10)) +
  labs(title="Supplementary Figure: renal phenotype family MR",
       subtitle="All nine primary estimates: BH-FDR q=0.685\nNone was statistically significant",
       x="Odds ratio for FinnGen aortic dissection (log scale)", y=NULL,
       caption="Points are primary IVW estimates; bars are 95% confidence intervals.\nNon-significance does not prove absence of a clinically relevant effect.") +
  theme_pub + theme(plot.subtitle=element_text(size=9), plot.caption=element_text(hjust=0))
save_pair(pS, "Supplementary_Figure_renal_MR_forest", 150, 115)

writeLines(capture.output(sessionInfo()), file.path(fig_dir, "figure_sessionInfo.txt"))
cat("FIGURES_COMPLETE", fig_dir, "\n")
