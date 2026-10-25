# Re-render archived plot-only generators into an isolated staging directory.
# Inputs are existing analysis outputs; no fit, permutation or cohort selection is run.
audit<-normalizePath('CLAIM_CITATION_AUDIT')
staging_root<-normalizePath(file.path(audit,'font_export_inputs/AD-CKD'))
dir.create(file.path(staging_root,'project/results/figures'),recursive=TRUE,showWarnings=FALSE)
embedded_pdf<-function(filename,width,height,bg,...) {
 stopifnot(capabilities('aqua'))
 grDevices::quartz(type='pdf',file=filename,width=width,height=height,bg=bg)
}
scripts<-c('assemble_main_figure3.R','assemble_supplementary_figure_s5.R','assemble_supplementary_figures_s6_s7.R','assemble_supplementary_figure_s8.R')
if('s8-only' %in% commandArgs(trailingOnly=TRUE)) scripts<-tail(scripts,1)
if('s5-only' %in% commandArgs(trailingOnly=TRUE)) scripts<-'assemble_supplementary_figure_s5.R'
for(name in scripts) {
 lines<-readLines(file.path(audit,'recovered/AD-CKD/project/scripts',name))
 lines[grepl('^root <-|^args<-',lines)]<-'root <- staging_root'
 txt<-paste(lines,collapse='\n')
 txt<-gsub('device = "pdf"','device = embedded_pdf',txt,fixed=TRUE)
 txt<-gsub('device="pdf"','device=embedded_pdf',txt,fixed=TRUE)
 if(name=='assemble_supplementary_figure_s5.R') {
  txt<-gsub('colour=outcome_short','colour=outcome_short,shape=outcome_short',txt,fixed=TRUE)
  txt<-gsub('colour=qc_status','colour=qc_status,shape=qc_status',txt,fixed=TRUE)
  txt<-gsub('colour=trait_label','colour=trait_label,shape=trait_label',txt,fixed=TRUE)
  txt<-gsub('colour=hypothesis','colour=hypothesis,shape=hypothesis',txt,fixed=TRUE)
  txt<-gsub('colour=NULL','colour=NULL,shape=NULL',txt,fixed=TRUE)
  txt<-gsub('colour="QC status"','colour="QC status",shape="QC status"',txt,fixed=TRUE)
  txt<-gsub('yintercept=ivw,colour="#D55E00",linewidth=.5','yintercept=ivw,colour="#D55E00",linewidth=.5,linetype="dashed"',txt,fixed=TRUE)
  txt<-gsub('orange line is the full IVW estimate','dashed line is the full IVW estimate',txt,fixed=TRUE)
 }
 if(name=='assemble_supplementary_figure_s8.R') {
  txt<-gsub('scale_fill_gradient2(low="#0072B2",mid="white",high="#D55E00",','scale_fill_viridis_c(option="E",direction=-1,breaks=c(-2.5,0,2.5),',txt,fixed=TRUE)
  txt<-gsub('Within-dataset z scores; not immune-cell proportions','Within-dataset z scores; not immune-cell proportions\\nLight = negative; dark = positive; limits at +/-2.5',txt,fixed=TRUE)
  txt<-gsub('height=210','height=240',txt,fixed=TRUE)
  txt<-gsub('heights=c(1.0,.9,1.25)','heights=c(1.0,.9,1.65)',txt,fixed=TRUE)
  # The model steps, not raw-expression preprocessing, were nested.
  txt<-gsub('Nested OOF','Model OOF',txt,fixed=TRUE)
 }
 eval(parse(text=txt),envir=new.env(parent=globalenv()))
}
