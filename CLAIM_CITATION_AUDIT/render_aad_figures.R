#!/usr/bin/env Rscript
# Label-only revision of original R figures using archived frozen Source Data.
audit <- normalizePath("CLAIM_CITATION_AUDIT")
# Use the existing framework R 4.5 library when running R 4.5; do not load
# compiled packages from the separate Homebrew R 4.6 installation.
if (getRversion() >= "4.6.0") .libPaths(c(file.path(audit,"Rlib"), .libPaths()))
src <- file.path(audit,"recovered/AD-CKD/project")
out <- file.path(audit,"aad_revision_candidate/figures")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
frozen <- file.path(src,"results/figures")
read_frozen <- function(name) read.delim(file.path(frozen,name),check.names=FALSE,stringsAsFactors=FALSE)
embedded_pdf <- function(filename,width,height,bg,...) {
  stopifnot(capabilities("aqua"))
  grDevices::quartz(type="pdf",file=filename,width=width,height=height,bg=bg)
}
rewrite <- function(text, old, new) {
  stopifnot(grepl(old,text,fixed=TRUE)); sub(old,new,text,fixed=TRUE)
}
run_original <- function(name, prepare) {
  if ('f2-only' %in% commandArgs(trailingOnly=TRUE) && name!='assemble_main_figure2.R') return(invisible(NULL))
  if ('s1-only' %in% commandArgs(trailingOnly=TRUE) && name!='assemble_supplementary_figure_s1.R') return(invisible(NULL))
  if ('s2-only' %in% commandArgs(trailingOnly=TRUE) && name!='assemble_supplementary_figure_s2.R') return(invisible(NULL))
  text <- paste(readLines(file.path(src,"scripts",name)),collapse="\n")
  start <- regexpr("theme_pub <-|ids <-|theme_qc <-",text)[1]
  text <- substr(text,start,nchar(text))
  text <- prepare(text)
  text <- gsub('device = "pdf"','device = embedded_pdf',text,fixed=TRUE)
  # Display labels only; frozen group identifiers and numerical inputs stay intact.
  text <- gsub('Renal-dysfunction contrast','Uraemic calcification contrast',text,fixed=TRUE)
  text <- gsub('Renal-dysfunction model PCA','Uraemic calcification model PCA',text,fixed=TRUE)
  text <- gsub('Renal-dysfunction top-DEG heatmap','Uraemic calcification top-DEG heatmap',text,fixed=TRUE)
  text <- gsub('"Not reported"','"Not encoded"',text,fixed=TRUE)
  # Display terminology only; identifiers such as logFC_ad remain stable.
  text <- gsub("Experimental AD","Experimental AAD",text,fixed=TRUE)
  text <- gsub("BAPN + AD","BAPN + AAD",text,fixed=TRUE)
  text <- gsub("BAPN, no AD","BAPN NonAD",text,fixed=TRUE)
  text <- gsub("down AD","down AAD",text,fixed=TRUE)
  text <- gsub("up AD","up AAD",text,fixed=TRUE)
  e <- new.env(parent=globalenv()); e$fig_dir <- out; e$out <- out
  eval(parse(text=text),envir=e)
}
library(ggplot2);library(patchwork)
run_original("assemble_main_figure1.R",function(text) {
  text <- gsub('"Renal dysfunction"','"Uraemic calcification"',text,fixed=TRUE)
  text <- rewrite(text,'labs(title = "Animal aortic transcriptomic contrasts"',
    'scale_y_discrete(labels = function(x) sub("Uraemic calcification", "Uraemic\\ncalcification", x, fixed=TRUE)) + labs(title = "Animal aortic transcriptomic contrasts"')
  text <- gsub('; 16 versus 20; sex-adjusted','; 16 versus 20\\nSex-adjusted',text,fixed=TRUE)
  text <- rewrite(text,'rat <- read.csv(file.path(geo, "GSE146638_all_genes.csv"), check.names = FALSE)',
    'rat <- read_frozen("Figure1_animal_transcriptomic_discovery_source_data_B_renal_volcano.tsv")')
  text <- rewrite(text,'mouse <- read.csv(file.path(geo, "GSE325608_all_genes.csv"), check.names = FALSE)',
    'mouse <- read_frozen("Figure1_animal_transcriptomic_discovery_source_data_C_AD_volcano.tsv")')
  a <- regexpr('rat_map <-',text,fixed=TRUE)[1]
  b <- regexpr('p_d <- ggplot',text,fixed=TRUE)[1]
  text <- paste0(substr(text,1,a-1),
    'overlap <- read_frozen("Figure1_animal_transcriptomic_discovery_source_data_D_overlap.tsv")\n',
    'direction <- read_frozen("Figure1_animal_transcriptomic_discovery_source_data_D_direction_summary.tsv")\n',
    'direction$direction <- gsub(" AD", " AAD", direction$direction, fixed=TRUE)\n',
    'overlap$direction <- gsub(" AD", " AAD", overlap$direction, fixed=TRUE)\n',
    'direction_order <- c("Up in both", "Down in both", "Up renal / down AAD", "Down renal / up AAD")\n',
    'sens_n <- unique(direction$sensitivity_deg_count); stopifnot(length(sens_n)==1,sens_n==1963)\n',
    substr(text,b,nchar(text)))
  text
})
run_original("assemble_main_figure2.R",function(text) {
  a <- regexpr('conv <-',text,fixed=TRUE)[1]; b <- regexpr('ora$pathway <-',text,fixed=TRUE)[1]
  text <- paste0(substr(text,1,a-1),'ora <- read_frozen("Figure2_pathway_convergence_human_support_source_data_A_ORA.tsv")\n',
    'ora$model[ora$model!="Renal dysfunction"] <- "AAD lesions"\n',substr(text,b,nchar(text)))
  a <- regexpr('gsea <-',text,fixed=TRUE)[1]; b <- regexpr('gsea$pathway <-',text,fixed=TRUE)[1]
  text <- paste0(substr(text,1,a-1),'gsea <- read_frozen("Figure2_pathway_convergence_human_support_source_data_B_animal_GSEA.tsv")\n',
    'gsea$model[gsea$model!="Renal dysfunction"] <- "AAD lesions"\n',substr(text,b,nchar(text)))
  text <- rewrite(text,'human <- read.csv(file.path(geo, "human_validation", "human_validation_six_pathway_GSEA.csv"), check.names = FALSE)\nhuman$pathway_id <- human$GO_ID',
    'human <- read_frozen("Figure2_pathway_convergence_human_support_source_data_C_human_GSEA.tsv")')
  text <- gsub('"Aortic dissection"','"AAD lesions"',text,fixed=TRUE)
  text <- gsub('scale_colour_manual(values = c("Renal dysfunction"',
    'scale_colour_manual(labels = c("Renal dysfunction"="Uraemic calcification", "AAD lesions"="AAD lesions"), values = c("Renal dysfunction"',text,fixed=TRUE)
  text <- gsub('scale_shape_manual(values = c("Renal dysfunction"',
    'scale_shape_manual(labels = c("Renal dysfunction"="Uraemic calcification", "AAD lesions"="AAD lesions"), values = c("Renal dysfunction"',text,fixed=TRUE)
  # Human cases are dissection, not the composite animal AAD phenotype.
  text <- rewrite(text,'group = c("AAD", "Control"','group = c("AD", "Control"')
  text <- rewrite(text,'design$label <- sprintf("%s\\n%s\\nn=%d", design$group, design$tissue, design$n)',
    'design$label <- sprintf("%s\\n%s\\nn=%d", design$group, sub("Internal mammary artery", "Internal mammary\\nartery", design$tissue, fixed=TRUE), design$n)')
  # Redundant encoding: outline shape identifies dataset; fill identifies BH status.
  text <- rewrite(text,'aes(NES, pathway, colour = model, shape = passed)',
    'aes(NES, pathway, colour = model, shape = model, fill = passed)')
  text <- rewrite(text,'aes(NES, pathway, colour = dataset_label, shape = passed)',
    'aes(NES, pathway, colour = dataset_label, shape = dataset_label, fill = passed)')
  old_scale <- 'scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "BH q<0.05", `FALSE` = "BH q>=0.05"))'
  fill_scale <- 'scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white"), labels = c(`TRUE` = "BH q<0.05", `FALSE` = "BH q>=0.05"))'
  text <- rewrite(text,old_scale,paste0('scale_shape_manual(values = c("Renal dysfunction" = 21, "AAD lesions" = 24), labels = c("Renal dysfunction"="Uraemic calcification", "AAD lesions"="AAD lesions")) + ',fill_scale))
  text <- rewrite(text,old_scale,paste0('scale_shape_manual(values = c(21, 24, 22)) + ',fill_scale))
  text <- rewrite(text,'guides(colour = "none", shape = guide_legend(nrow = 1))',
    'guides(colour = "none", shape = "none", fill = guide_legend(nrow = 1, override.aes = list(shape = 21, colour = "black")))')
  text <- rewrite(text,'guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))',
    'guides(colour = "none", shape = guide_legend(nrow = 1, override.aes = list(fill = "white", colour = "black")), fill = guide_legend(nrow = 1, override.aes = list(shape = 21, colour = "black")))')
  text
})
cat("Label-only Figure 1 and Figure 2 candidates generated from frozen values.\n")
run_original("assemble_supplementary_figure_s1.R",function(text) {
  text <- rewrite(text,'colour = primary_group, shape = sex','fill = primary_group, shape = sex')
  text <- rewrite(text,'colour = "Group", shape = "Sex"','fill = "Group", shape = "Sex"')
  text <- rewrite(text,'scale_colour_manual(values = c("Control" = "grey40", "Uremic" = "#0072B2", "AAD" = "#D55E00"))',
    'scale_fill_manual(values = c("Control" = "white", "Uremic" = "#0072B2", "AAD" = "#D55E00")) + guides(fill = guide_legend(override.aes = list(shape = 21)))')
  text <- rewrite(text,'scale_shape_discrete(labels = function(x) ifelse(x == "not_reported", "Not reported", x))',
    'scale_shape_manual(values = c(Male = 21, Female = 24, not_reported = 21), labels = function(x) ifelse(x == "not_reported", "Not reported", x))')
  a <- regexpr('read_meta <-',text,fixed=TRUE)[1]; b <- regexpr('plot_pca <-',text,fixed=TRUE)[1]
  text <- paste0(substr(text,1,a-1),
    'pc1 <- read_frozen("Supplementary_Figure_S1_animal_QC_source_data_A_PCA.tsv")\n',
    'pc2 <- read_frozen("Supplementary_Figure_S1_animal_QC_source_data_B_PCA.tsv")\n',
    'pc1$primary_group <- factor(pc1$primary_group,levels=c("Control","Uremic"))\n',
    'pc2$primary_group <- factor(pc2$primary_group,levels=c("Control","AAD"))\n',
    'pc2$sex <- factor(pc2$sex,levels=c("Male","Female"))\n',
    '# Variance percentages preserved from the frozen published panel, not recomputed.\n',
    'attr(pc1,"pct") <- c(59.6,13.6); attr(pc2,"pct") <- c(59.0,5.1)\n',
    substr(text,b,nchar(text)))
  a <- regexpr('corr_long <-',text,fixed=TRUE)[1]; b <- regexpr('plot_corr <-',text,fixed=TRUE)[1]
  text <- paste0(substr(text,1,a-1),
    'cor1 <- read_frozen("Supplementary_Figure_S1_animal_QC_source_data_C_correlation.tsv")\n',
    'cor2 <- read_frozen("Supplementary_Figure_S1_animal_QC_source_data_D_correlation.tsv")\n',
    'restore_corr <- function(d,pc) {\n',
    '  meta <- pc[order(pc$primary_group,pc$sex,pc$expression_column),]\n',
    '  compact <- if (unique(d$dataset)=="GSE325608") sub("^read_count_","",meta$expression_column) else meta$expression_column\n',
    '  d$group_x <- factor(d$group_x,levels=levels(pc$primary_group)); d$group_y <- factor(d$group_y,levels=rev(levels(pc$primary_group)))\n',
    '  d$display_x <- factor(d$display_x,levels=compact); d$display_y <- factor(d$display_y,levels=rev(compact));d}\n',
    'cor1 <- restore_corr(cor1,pc1);cor2 <- restore_corr(cor2,pc2)\n',
    substr(text,b,nchar(text)))
  a <- regexpr('voom <- data.frame',text,fixed=TRUE)[1]; b <- regexpr('p_e <-',text,fixed=TRUE)[1]
  text <- paste0(substr(text,1,a-1),
    'voom <- read_frozen("Supplementary_Figure_S1_animal_QC_source_data_E_voom_points.tsv")\n',
    'voom_line <- read_frozen("Supplementary_Figure_S1_animal_QC_source_data_E_voom_line.tsv")\n',substr(text,b,nchar(text)))
  text <- gsub('AD-model sample','AAD-model sample',text,fixed=TRUE)
  gsub('experimental AD model','experimental AAD model',text,fixed=TRUE)
})
run_original("assemble_supplementary_figure_s2.R",function(text) {
  text <- rewrite(text,'`FALSE` = "grey72", `TRUE` = "#CC79A7"','`FALSE` = "grey80", `TRUE` = "#6A3D9A"')
  text <- gsub('scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,', 'scale_fill_viridis_c(option = "E", direction = -1, breaks = c(-2.5, 0, 2.5),',text,fixed=TRUE)
  text <- gsub('Top 30 DEGs by absolute moderated t statistic','Top 30 DEGs by absolute moderated t statistic\\nLight = negative; dark = positive; limits at +/-2.5',text,fixed=TRUE)
  a <- regexpr('zscore <-',text,fixed=TRUE)[1]; b <- regexpr('heat_plot <-',text,fixed=TRUE)[1]
  text <- paste0(substr(text,1,a-1),substr(text,b,nchar(text)))
  text <- rewrite(text,'hd1 <- heat_data(h1, m1, "GSE146638"); hd2 <- heat_data(h2, m2, "GSE325608")',paste0(
    'restore_heat <- function(d,groups) {d$gene <- factor(d$gene,levels=rev(unique(d$gene))); d$sample <- factor(d$sample,levels=unique(d$sample));d$group <- factor(d$group,levels=groups);d}\n',
    'hd1 <- restore_heat(read_frozen("Supplementary_Figure_S2_animal_expression_detail_source_data_A_heatmap.tsv"),c("Control","Uremic"))\n',
    'hd2 <- restore_heat(read_frozen("Supplementary_Figure_S2_animal_expression_detail_source_data_B_heatmap.tsv"),c("Control","AAD"))'))
  text <- rewrite(text,'sens <- read.csv(file.path(geo, "GSE325608_AAD_vs_NonAD_sensitivity.csv"), check.names = FALSE)',
    'sens <- read_frozen("Supplementary_Figure_S2_animal_expression_detail_source_data_C_sensitivity_volcano.tsv")')
  a <- regexpr('primary <-',text,fixed=TRUE)[1]; b <- regexpr('p_d <- ggplot',text,fixed=TRUE)[1]
  text <- paste0(substr(text,1,a-1),
    'compare <- read_frozen("Supplementary_Figure_S2_animal_expression_detail_source_data_D_effect_comparison.tsv")\n',
    'rho <- cor(compare$logFC_primary,compare$logFC_sensitivity,method="spearman",use="complete.obs")\n',substr(text,b,nchar(text)))
  text <- gsub('BAPN without AD','BAPN NonAD',text,fixed=TRUE)
  text <- gsub('Primary AD versus','Primary AAD versus',text,fixed=TRUE)
  text <- rewrite(text,'; 16 versus 12; sex-adjusted; 1,963 DEGs','; 16 versus 12\\nSex-adjusted; 1,963 DEGs')
  text
})
