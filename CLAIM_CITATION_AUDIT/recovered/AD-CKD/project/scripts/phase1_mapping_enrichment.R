#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, timeout = 180)
script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
script_path <- sub("^--file=", "", script_arg)
root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
setwd(root)
.libPaths(c("renv/library/macos/R-4.5/aarch64-apple-darwin20", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(biomaRt)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(pheatmap)
})

out <- "project/results/geo"
cache <- "project/data/geo/processed"
fig <- file.path(out, "figures")
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)

fetch_orthologs <- function(genes, dataset, prefix, cache_file) {
  if (file.exists(cache_file)) {
    message("Using cached ortholog mapping: ", cache_file)
    return(fread(cache_file))
  }
  mart <- useEnsembl(biomart = "genes", dataset = dataset)
  attrs <- c(
    "external_gene_name",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_type",
    "hsapiens_homolog_orthology_confidence"
  )
  chunks <- split(unique(genes), ceiling(seq_along(unique(genes)) / 500))
  ans <- rbindlist(lapply(seq_along(chunks), function(i) {
    message(prefix, " BioMart chunk ", i, "/", length(chunks))
    getBM(attributes = attrs, filters = "external_gene_name", values = chunks[[i]], mart = mart)
  }), fill = TRUE)
  setnames(ans, c("source_symbol", "human_symbol", "orthology_type", "orthology_confidence"))
  ans <- unique(ans[nzchar(source_symbol) & nzchar(human_symbol)])
  fwrite(ans, cache_file)
  ans
}

keep_one2one <- function(map) {
  unique(map[orthology_type == "ortholog_one2one" & orthology_confidence == 1,
             .(source_symbol, human_symbol, orthology_type, orthology_confidence)])
}

map_symbols <- function(symbols, map) {
  unique(map[source_symbol %in% symbols, human_symbol])
}

run_go <- function(genes, universe, label) {
  if (length(genes) < 10) return(data.table())
  x <- enrichGO(
    gene = genes, universe = universe, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
    minGSSize = 10, maxGSSize = 500, readable = FALSE
  )
  ans <- as.data.table(as.data.frame(x))
  if (nrow(ans)) ans[, `:=`(database = "GO_BP", dataset = label)]
  ans
}

to_entrez <- function(symbols) {
  x <- suppressMessages(bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  unique(x$ENTREZID)
}

run_kegg <- function(genes, universe, label, term2gene = NULL, term2name = NULL) {
  g <- to_entrez(genes)
  u <- to_entrez(universe)
  if (length(g) < 10) return(data.table())
  if (!is.null(term2gene) && !is.null(term2name)) {
    x <- enricher(gene = g, universe = u, TERM2GENE = term2gene, TERM2NAME = term2name,
                  pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1,
                  minGSSize = 10, maxGSSize = 500)
  } else {
    x <- enrichKEGG(
      gene = g, universe = u, organism = "hsa", keyType = "kegg",
      pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1,
      minGSSize = 10, maxGSSize = 500, use_internal_data = FALSE
    )
  }
  ans <- as.data.table(as.data.frame(x))
  if (nrow(ans)) ans[, `:=`(database = "KEGG", dataset = label)]
  ans
}

r1 <- fread(file.path(out, "GSE146638_all_genes.csv"))
r2 <- fread(file.path(out, "GSE325608_all_genes.csv"))
d1 <- r1[adj.P.Val < 0.05 & abs(logFC) > 0.5]
d2 <- r2[adj.P.Val < 0.05 & abs(logFC) > 0.5]

rat_all <- fetch_orthologs(r1$gene_symbol, "rnorvegicus_gene_ensembl", "rat", file.path(cache, "rat_human_biomart_all.csv"))
mouse_all <- fetch_orthologs(r2$gene_symbol, "mmusculus_gene_ensembl", "mouse", file.path(cache, "mouse_human_biomart_all.csv"))
fwrite(rat_all, file.path(out, "ortholog_rat_human_all.csv"))
fwrite(mouse_all, file.path(out, "ortholog_mouse_human_all.csv"))
rat_map <- keep_one2one(rat_all)
mouse_map <- keep_one2one(mouse_all)
fwrite(rat_map, file.path(out, "ortholog_rat_human_1to1.csv"))
fwrite(mouse_map, file.path(out, "ortholog_mouse_human_1to1.csv"))

uremic <- map_symbols(d1$gene_symbol, rat_map)
ad <- map_symbols(d2$gene_symbol, mouse_map)
univ_uremic <- map_symbols(r1$gene_symbol, rat_map)
univ_ad <- map_symbols(r2$gene_symbol, mouse_map)

intersection <- sort(intersect(uremic, ad))
fwrite(data.table(gene = intersection), file.path(out, "intersection_genes.csv"))
fwrite(data.table(dataset = "GSE146638", human_gene = sort(uremic)), file.path(out, "GSE146638_human_1to1_DEG.csv"))
fwrite(data.table(dataset = "GSE325608", human_gene = sort(ad)), file.path(out, "GSE325608_human_1to1_DEG.csv"))

go1 <- run_go(uremic, univ_uremic, "UREMIC")
go2 <- run_go(ad, univ_ad, "AD")
fwrite(go1, file.path(out, "GOBP_UREMIC.csv"))
fwrite(go2, file.path(out, "GOBP_AD.csv"))

kegg_status <- "VERIFIED"
kegg_error <- ""
k1 <- k2 <- data.table()
tryCatch({
  download_dir <- path.expand("~/Downloads")
  link_candidates <- c(file.path(cache, "kegg_pathway_gene_link.tsv"),
                       file.path(download_dir, c("kegg_pathway_gene_link.tsv", "KEGG_pathway-gene_link.tsv")))
  name_candidates <- c(file.path(cache, "kegg_pathway_names.tsv"),
                       file.path(download_dir, c("kegg_pathway_names.tsv", "KEGG _pahtway_names.tsv")))
  link_file <- link_candidates[file.exists(link_candidates)][1]
  name_file <- name_candidates[file.exists(name_candidates)][1]
  if (file.exists(link_file) && file.exists(name_file) && length(readLines(link_file, n = 11000)) >= 10000) {
    links <- fread(link_file, header = FALSE, col.names = c("term", "gene"))
    names <- fread(name_file, header = FALSE, col.names = c("term", "name"))
    links[, `:=`(term = sub("^path:", "", term), gene = sub("^hsa:", "", gene))]
    names[, `:=`(term = sub("^path:", "", term), name = sub(" - Homo sapiens \\(human\\)$", "", name))]
    k1 <- run_kegg(uremic, univ_uremic, "UREMIC", links, names)
    k2 <- run_kegg(ad, univ_ad, "AD", links, names)
    kegg_status <- "VERIFIED_LOCAL_KEGG_REST_SNAPSHOT"
  } else {
    k1 <- run_kegg(uremic, univ_uremic, "UREMIC")
    k2 <- run_kegg(ad, univ_ad, "AD")
  }
}, error = function(e) {
  kegg_status <<- "UNAVAILABLE"
  kegg_error <<- conditionMessage(e)
})
fwrite(k1, file.path(out, "KEGG_UREMIC.csv"))
fwrite(k2, file.path(out, "KEGG_AD.csv"))

standardize <- function(a, b, db) {
  if (!nrow(a) || !nrow(b)) return(data.table())
  aa <- a[, .(ID, Description_UREMIC = Description, FDR_UREMIC = p.adjust, GeneRatio_UREMIC = GeneRatio)]
  bb <- b[, .(ID, Description_AD = Description, FDR_AD = p.adjust, GeneRatio_AD = GeneRatio)]
  z <- merge(aa, bb, by = "ID")
  z[, `:=`(database = db, Description = fifelse(nzchar(Description_UREMIC), Description_UREMIC, Description_AD))]
  z
}

shared <- rbindlist(list(standardize(go1, go2, "GO_BP"), standardize(k1, k2, "KEGG")), fill = TRUE)
shared <- shared[FDR_UREMIC < 0.01 & FDR_AD < 0.01]
focus_pattern <- paste(c("extracellular matrix", "collagen", "ossification", "calcification", "oxidative stress",
                         "TGF", "smooth muscle", "AGE-RAGE", "blood vessel", "elastic fiber", "proteoglycan"), collapse = "|")
if (nrow(shared)) {
  shared[, focus_match := grepl(focus_pattern, Description, ignore.case = TRUE)]
  shared[, max_FDR := pmax(FDR_UREMIC, FDR_AD)]
  setorder(shared, -focus_match, max_FDR)
}
fwrite(shared, file.path(out, "convergence_pathways.csv"))

if (nrow(shared)) {
  top <- head(shared, 30)
  hm <- as.matrix(top[, .(-log10(pmax(FDR_UREMIC, 1e-300)), -log10(pmax(FDR_AD, 1e-300)))])
  rownames(hm) <- make.unique(substr(top$Description, 1, 70))
  colnames(hm) <- c("Uremic", "AAD")
  ann <- data.frame(Focus = ifelse(top$focus_match, "Hypothesis-related", "Other"), row.names = rownames(hm))
  pheatmap(hm, cluster_cols = FALSE, annotation_row = ann,
           filename = file.path(fig, "convergence_pathway_heatmap.png"), width = 9, height = max(5, 0.25 * nrow(hm) + 2))
}

scored <- if (nrow(shared)) sum(shared$focus_match) else 0
total <- nrow(shared)
support <- if (scored >= 3 && total >= 5) "SUPPORTED" else if (total >= 5) "PARTIAL_SUPPORT" else "WEAK_CONVERGENCE"
summary <- data.table(
  metric = c("rat_tested_genes", "rat_1to1_universe", "rat_DEG", "rat_human_1to1_DEG",
             "mouse_tested_genes", "mouse_1to1_universe", "mouse_DEG", "mouse_human_1to1_DEG",
             "intersection_genes", "shared_FDR_lt_0.01_pathways", "focus_pathways", "support", "KEGG_status", "KEGG_error"),
  value = c(nrow(r1), length(univ_uremic), nrow(d1), length(uremic), nrow(r2), length(univ_ad), nrow(d2), length(ad),
            length(intersection), total, scored, support, kegg_status, kegg_error)
)
fwrite(summary, file.path(out, "phase1_convergence_summary.tsv"), sep = "\t")
print(summary)
