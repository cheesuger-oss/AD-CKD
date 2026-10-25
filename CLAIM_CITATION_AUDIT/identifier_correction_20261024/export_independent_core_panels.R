# Reuse original corrected R plot objects inside an isolated archive extraction.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>=1)
mode <- if(length(args)>1) args[2] else 'core'
revision_home <- normalizePath(args[1])
stopifnot(requireNamespace('ggplot2',quietly=TRUE),requireNamespace('patchwork',quietly=TRUE),capabilities('aqua'))
panel_home <- file.path(revision_home,'independent_panels');dir.create(panel_home,showWarnings=FALSE)
scratch <- normalizePath(file.path(revision_home,'../..'))
setwd(scratch)
panel_records <- list()
capture_panels <- function(name,e) {
 maps <- list(
  'assemble_main_figure1.R'=c(A='p_a',B='p_d',C='p_b',D='p_c'),
  'assemble_main_figure2.R'=c(A='p_a',B='p_b',C='p_c',D='p_d'),
  'assemble_supplementary_figure_s1.R'=c(A='p_a',B='p_b',C='p_c',D='p_d',E='p_e'),
  'assemble_supplementary_figure_s2.R'=c(A='p_a',B='p_b',C='p_c',D='p_d'),
  'f3'=c(A='p_a',B='p_b',C='p_c'),'f4'=c(A='p_a',B='p_b',C='p_c'),
  's4a'=c(A='p1',B='p2',C='p3',D='p4'),'s4b'=c(A='p1',B='p2',C='p3',D='p4'),
  's5a'=c(A='p1',B='p2',C='p3'),'s5b'=c(A='p4',B='p5'),'s5c'=c(A='p6',B='p7',C='p8',D='p9'),
  's6'=c(A='p1',B='p2',C='p3'),'s7'=c(A='p4',B='p5',C='p6',D='p7'),
  's8'=c(A='p1',B='p2',C='p3',D='p4',E='p5',F='p6'))
 ids <- c('assemble_main_figure1.R'='Figure1','assemble_main_figure2.R'='Figure2',
  'assemble_supplementary_figure_s1.R'='FigureS1','assemble_supplementary_figure_s2.R'='FigureS2',
  f3='Figure3',f4='Figure4',s4a='FigureS4a',s4b='FigureS4b',s5a='FigureS5a',s5b='FigureS5b',s5c='FigureS5c',s6='FigureS6',s7='FigureS7',s8='FigureS8')
 for (letter in names(maps[[name]])) {
  object <- maps[[name]][[letter]];original <- get(object,envir=e)
  p <- original + ggplot2::labs(tag=letter) + ggplot2::theme(plot.tag=ggplot2::element_text(face='bold',size=10),plot.margin=ggplot2::margin(8,8,8,8))
  if(name=='assemble_main_figure2.R' && letter=='B') p <- p + ggplot2::guides(shape=ggplot2::guide_legend(nrow=1,override.aes=list(fill='white',colour='black')))
  # Shared-axis suppression is valid in the composite, but not in separate files.
  if(name=='assemble_main_figure2.R' && letter %in% c('B','C')) p <- p + ggplot2::guides(y=ggplot2::guide_axis()) + ggplot2::theme(axis.text.y=ggplot2::element_text(size=8,colour='black'),axis.text.y.left=ggplot2::element_text(size=8,colour='black'),axis.ticks.y=ggplot2::element_line())
  stopifnot(identical(ggplot2::ggplot_build(original)$data,ggplot2::ggplot_build(p)$data))
  stem <- paste0(ids[[name]],'_',letter)
  height <- if ((ids[[name]]=='FigureS1' && letter %in% c('C','D')) || (ids[[name]]=='FigureS2' && letter %in% c('A','B'))) 145 else 110
  if(name %in% c('s4a','s4b') && letter %in% c('B','D')) height<-145
  if(name=='s8' && letter %in% c('E','F')) height<-145
  pdf <- file.path(panel_home,paste0(stem,'.pdf'))
  grDevices::quartz(type='pdf',file=pdf,width=183/25.4,height=height/25.4,bg='white');print(p);dev.off()
  stopifnot(system2('/opt/homebrew/bin/pdftoppm',c('-scale-to','1700','-singlefile','-png',shQuote(pdf),shQuote(file.path(panel_home,stem))))==0)
  panel_records[[stem]] <<- data.frame(panel=stem,source_script=name,original_plot_object=object,width_mm=183,height_mm=height,built_layer_data_identical=TRUE,visual_review='PENDING')
  cat(stem,'exported\n')
 }
}
if(mode=='core') {
runner <- readLines(file.path(revision_home,'render_corrected_core_figures.R'))
needle <- 'eval(parse(text=txt), envir=globalenv())'
at <- which(runner==needle);stopifnot(length(at)==1)
runner[at] <- paste0('txt <- sub("eval(parse(text=text),envir=e)", "eval(parse(text=text),envir=e);capture_panels(name,e)",txt,fixed=TRUE)\n',needle)
eval(parse(text=paste(runner,collapse='\n')),envir=globalenv())
} else {
 library(ggplot2);library(patchwork)
 audit <- normalizePath('CLAIM_CITATION_AUDIT')
 scale_fill_viridis_c <- function(...) ggplot2::scale_fill_viridis_c(...,guide=guide_colourbar(display='rectangles'))
 scale_colour_viridis_c <- function(...) ggplot2::scale_colour_viridis_c(...,guide=guide_colourbar(display='rectangles'))
 txt<-paste(readLines(file.path(audit,'reexport_remaining_pdf_fonts.R')),collapse='\n')
 txt<-sub("scripts<-c('assemble_main_figure3.R','assemble_supplementary_figure_s5.R','assemble_supplementary_figures_s6_s7.R','assemble_supplementary_figure_s8.R')", "scripts<-c('assemble_main_figure3.R','assemble_supplementary_figure_s5.R','assemble_supplementary_figures_s6_s7.R')",txt,fixed=TRUE)
 old<-'eval(parse(text=txt),envir=new.env(parent=globalenv()))'
 stopifnot(grepl(old,txt,fixed=TRUE))
 txt<-sub(old,paste0("e<-new.env(parent=globalenv());eval(parse(text=txt),envir=e); ",
 "keys<-switch(name,assemble_main_figure3.R='f3',assemble_supplementary_figure_s5.R=c('s5a','s5b','s5c'),assemble_supplementary_figures_s6_s7.R=c('s6','s7'));for(key in keys)capture_panels(key,e)"),txt,fixed=TRUE)
 eval(parse(text=txt),envir=globalenv())
 txt<-paste(readLines(file.path(audit,'repair_figure4_alignment.R')),collapse='\n')
 txt<-sub('root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")','root <- normalizePath(".")',txt,fixed=TRUE)
 txt<-sub('file.path(analysis_dir, "Figure4_immune_classifier_robustness_source_data_A_immune.tsv")','file.path(audit,"identifier_correction_20261024/corrected_figure_source_data/Figure4_immune_classifier_robustness_source_data_A_immune.tsv")',txt,fixed=TRUE)
 e<-new.env(parent=globalenv());eval(parse(text=txt),envir=e);capture_panels('f4',e)
 txt<-paste(readLines(file.path(audit,'repair_s4_frozen_layout.R')),collapse='\n')
 txt<-sub('root <- normalizePath(commandArgs(trailingOnly=TRUE)[1])','root <- normalizePath(".")',txt,fixed=TRUE)
 txt<-sub("for(panel in c('pca','cor','volcano','heat'))", "capture_panels(paste0('s4',suffix),environment());for(panel in c('pca','cor','volcano','heat'))",txt,fixed=TRUE)
 eval(parse(text=txt),envir=new.env(parent=globalenv()))
 txt<-paste(readLines(file.path(revision_home,'render_corrected_s8.R')),collapse='\n')
 txt<-sub("eval(parse(text=txt),envir=globalenv())","e<-new.env(parent=globalenv());eval(parse(text=txt),envir=e);capture_panels('s8',e)",txt,fixed=TRUE)
 # S8 supplies its vector guide explicitly, so do not add a second guide argument.
 scale_fill_viridis_c <- ggplot2::scale_fill_viridis_c
 eval(parse(text=txt),envir=globalenv())
}
write.table(do.call(rbind,panel_records),file.path(panel_home,paste0(mode,'_panel_export.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(panel_home,paste0(mode,'_export_sessionInfo.txt')))
writeLines(scratch,file.path(panel_home,paste0(mode,'_isolated_export_location.txt')))
