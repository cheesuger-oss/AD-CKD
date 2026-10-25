#!/usr/bin/env Rscript
audit <- normalizePath('CLAIM_CITATION_AUDIT')
correction <- file.path(audit,'identifier_correction_20261024')
root <- normalizePath(file.path(audit,'font_export_inputs/AD-CKD'))
txt <- paste(readLines(file.path(audit,'recovered/AD-CKD/project/scripts/assemble_supplementary_figure_s8.R')),collapse='\n')
txt <- sub('args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))normalizePath(args[[1]])else normalizePath(".")','',txt,fixed=TRUE)
txt <- sub('out<-file.path(root,"project","results","figures")','out<-file.path(correction,"figures")',txt,fixed=TRUE)
txt <- sub('file.path(src,"immune_program_scores.tsv")','file.path(correction,"corrected_figure_source_data/Supplementary_Figure_S8_source_data_A_immune_scores.tsv")',txt,fixed=TRUE)
txt <- gsub('scale_fill_gradient2(low="#0072B2",mid="white",high="#D55E00",',
 'scale_fill_viridis_c(option="E",direction=-1,breaks=c(-2.5,0,2.5),guide=guide_colourbar(display="rectangles"),',txt,fixed=TRUE)
txt <- gsub('Within-dataset z scores; not immune-cell proportions',
 'Within-dataset z scores; not immune-cell proportions\\nLight = negative; dark = positive; limits at +/-2.5',txt,fixed=TRUE)
txt <- gsub('height=210','height=240',txt,fixed=TRUE)
txt <- gsub('heights=c(1.0,.9,1.25)','heights=c(1.0,.9,1.65)',txt,fixed=TRUE)
txt <- gsub('Nested OOF','Model OOF',txt,fixed=TRUE)
txt <- gsub('device="pdf"','device=embedded_pdf',txt,fixed=TRUE)
embedded_pdf <- function(filename,width,height,bg,...) {
 stopifnot(capabilities('aqua'))
 grDevices::quartz(type='pdf',file=filename,width=width,height=height,bg=bg)
}
eval(parse(text=txt),envir=globalenv())
