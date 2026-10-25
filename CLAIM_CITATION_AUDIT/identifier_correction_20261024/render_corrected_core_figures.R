#!/usr/bin/env Rscript
# Run from the revision root. Historical source files remain unchanged.
audit <- normalizePath('CLAIM_CITATION_AUDIT')
correction <- file.path(audit, 'identifier_correction_20261024')
stopifnot(file.exists(file.path(correction, 'statistical_recalculation_audit.json')))
txt <- paste(readLines(file.path(audit, 'render_aad_figures.R')), collapse='\n')
replace_once <- function(x, old, new) {
  stopifnot(length(gregexpr(old,x,fixed=TRUE)[[1]])==1L, grepl(old,x,fixed=TRUE))
  sub(old,new,x,fixed=TRUE)
}
txt <- replace_once(txt, 'out <- file.path(audit,"aad_revision_candidate/figures")',
  'out <- file.path(correction,"figures")')
txt <- replace_once(txt,
  'read_frozen <- function(name) read.delim(file.path(frozen,name),check.names=FALSE,stringsAsFactors=FALSE)',
  paste0('read_frozen <- function(name) { p <- file.path(correction,"corrected_figure_source_data",name); ',
  'if (!file.exists(p)) p <- file.path(frozen,name); ',
  'read.delim(p,check.names=FALSE,stringsAsFactors=FALSE) }'))
txt <- replace_once(txt, 'attr(pc1,"pct") <- c(59.6,13.6); attr(pc2,"pct") <- c(59.0,5.1)',
  'attr(pc1,"pct") <- read.delim(file.path(correction,"results/corrected_PCA_variance.tsv"))$percent; attr(pc2,"pct") <- c(59.0,5.1)')
txt <- gsub('Variance percentages preserved from the frozen published panel, not recomputed.',
  'Rat variance percentages recomputed from corrected matrix; mouse unchanged.',txt,fixed=TRUE)
txt <- gsub('Label-only Figure 1 and Figure 2 candidates generated from frozen values.',
  'Identifier-corrected Figure 1 and Figure 2 generated; package release pending.',txt,fixed=TRUE)
# Wrap the existing long title inside its own panel; no wording or data change.
txt <- replace_once(txt, 'text <- prepare(text)', paste0('text <- prepare(text)\n',
 '  text <- gsub("Selected programmes overlapped in the animal models", ',
 '"Selected programmes overlapped\\\\nin the animal models", text, fixed=TRUE)'))
# Explicit rectangle guides avoid embedded bitmap colourbars in vector PDFs.
library(ggplot2)
scale_fill_viridis_c <- function(...) ggplot2::scale_fill_viridis_c(...,guide=guide_colourbar(display='rectangles'))
scale_fill_gradient2 <- function(...) ggplot2::scale_fill_gradient2(...,guide=guide_colourbar(display='rectangles'))
scale_fill_gradient <- function(...) ggplot2::scale_fill_gradient(...,guide=guide_colourbar(display='rectangles'))
eval(parse(text=txt), envir=globalenv())
writeLines(capture.output(sessionInfo()),file.path(correction,'core_figure_sessionInfo.txt'))
