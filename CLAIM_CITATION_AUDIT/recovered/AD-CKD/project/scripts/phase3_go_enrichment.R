#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
.libPaths(c(file.path(root, "project", "Rlib-4.5"),
            file.path(root, "renv", "library", "macos", "R-4.5", "aarch64-apple-darwin20"),
            .libPaths()))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

outdir <- file.path(root, "project", "results", "network")
genes <- unique(scan(file.path(outdir, "union_genes.txt"), what = character(), quiet = TRUE))
background <- sort(unique(keys(org.Hs.eg.db, keytype = "SYMBOL")))
writeLines(background, file.path(outdir, "go_background_symbols.txt"))
background_entrez <- sort(unique(keys(org.Hs.eg.db, keytype = "ENTREZID")))
writeLines(background_entrez, file.path(outdir, "go_background_entrez.txt"))

mapping <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = FALSE)
write.table(mapping, file.path(outdir, "go_symbol_entrez_mapping.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")
mapped <- unique(mapping$ENTREZID[!is.na(mapping$ENTREZID)])

ego <- enrichGO(
  gene = mapped, universe = background_entrez, keyType = "ENTREZID", OrgDb = org.Hs.eg.db,
  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  readable = TRUE
)
result <- as.data.frame(ego)
write.csv(result, file.path(outdir, "enrichment_table.csv"), row.names = FALSE, na = "")

metadata <- data.frame(
  key = c("ontology", "input_symbol_count", "mapped_entrez_count", "background_definition",
          "background_entrez_count", "p_adjust_method", "reported_term_count", "interpretation"),
  value = c("GO Biological Process", length(genes), length(mapped),
            "all SYMBOL keys in installed org.Hs.eg.db; explicit universe supplied to enrichGO",
            length(background_entrez), "BH", nrow(result),
            "over-representation is associative and does not establish causal pathways")
)
write.table(metadata, file.path(outdir, "enrichment_metadata.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("input", length(genes), "mapped", length(mapped), "background", length(background_entrez),
    "terms", nrow(result), "BH_q_lt_0.05", sum(result$p.adjust < 0.05, na.rm = TRUE), "\n")
