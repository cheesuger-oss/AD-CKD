# Export-only contract: preserve frozen scientific panels and colour mappings.
# Original R plotting workflow retained; only colourbar display becomes vector.
suppressPackageStartupMessages(library(ggplot2))
audit <- normalizePath('CLAIM_CITATION_AUDIT')
vector_out <- file.path(audit,'vector_export_stage')
dir.create(vector_out,showWarnings=FALSE)
scale_fill_viridis_c <- function(...) ggplot2::scale_fill_viridis_c(...,guide=guide_colourbar(display='rectangles'))
scale_colour_viridis_c <- function(...) ggplot2::scale_colour_viridis_c(...,guide=guide_colourbar(display='rectangles'))
ggsave <- function(filename, ...) {
  if(!grepl('\\.pdf$',filename) || grepl('Supplementary_Figure_S7_',filename)) return(invisible(NULL))
  ggplot2::ggsave(file.path(vector_out,basename(filename)), ...)
}
mode <- tail(commandArgs(trailingOnly=TRUE),1)
if(mode %in% c('s1-only','s2-only')) source(file.path(audit,'render_aad_figures.R'))
if(mode=='s4') source(file.path(audit,'repair_s4_frozen_layout.R'))
if(mode=='remaining') {
  txt<-paste(readLines(file.path(audit,'reexport_remaining_pdf_fonts.R')),collapse='\n')
  old<-"scripts<-c('assemble_main_figure3.R','assemble_supplementary_figure_s5.R','assemble_supplementary_figures_s6_s7.R','assemble_supplementary_figure_s8.R')"
  stopifnot(grepl(old,txt,fixed=TRUE))
  txt<-sub(old,"scripts<-c('assemble_supplementary_figures_s6_s7.R','assemble_supplementary_figure_s8.R')",txt,fixed=TRUE)
  eval(parse(text=txt),envir=globalenv())
}
