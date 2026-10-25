#!/usr/bin/env Rscript
audit <- normalizePath('CLAIM_CITATION_AUDIT')
correction <- file.path(audit,'identifier_correction_20261024')
replace_once <- function(x,old,new) {
 stopifnot(grepl(old,x,fixed=TRUE));sub(old,new,x,fixed=TRUE)
}
txt <- paste(readLines(file.path(audit,'repair_figure4_alignment.R')),collapse='\n')
txt <- replace_once(txt,'figure_dir <- file.path(root, "CLAIM_CITATION_AUDIT/aad_revision_candidate/figures")',
 'figure_dir <- file.path(correction,"figures")')
txt <- replace_once(txt,
 'file.path(analysis_dir, "Figure4_immune_classifier_robustness_source_data_A_immune.tsv")',
 'file.path(correction, "corrected_figure_source_data/Figure4_immune_classifier_robustness_source_data_A_immune.tsv")')
eval(parse(text=txt),envir=globalenv())
txt <- paste(readLines(file.path(audit,'render_figure5_current_synthesis.R')),collapse='\n')
txt <- replace_once(txt,'out <- "CLAIM_CITATION_AUDIT/aad_revision_candidate/figures"',
 'out <- file.path(correction,"figures")')
txt <- replace_once(txt,'549 shared differentially expressed genes and 155 shared pathways',
 '549 shared differentially expressed genes and 154 shared pathways')
eval(parse(text=txt),envir=globalenv())
