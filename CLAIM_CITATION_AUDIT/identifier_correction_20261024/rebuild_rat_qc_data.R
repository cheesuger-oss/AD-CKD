options(stringsAsFactors=FALSE)
library(data.table)
root<-normalizePath('CLAIM_CITATION_AUDIT/identifier_correction_20261024')
scripts<-file.path(root,'../recovered/AD-CKD/project/scripts')
load_function<-function(file,name){
 statements<-as.list(parse(file.path(scripts,file)))
 found<-Filter(function(x)is.call(x)&&as.character(x[[1]])=='<-'&&identical(x[[2]],as.name(name)),statements)
 stopifnot(length(found)==1);eval(found[[1]],envir=.GlobalEnv)
}
load_function('assemble_supplementary_figure_s1.R','pca_data')
load_function('assemble_supplementary_figure_s1.R','corr_long')
load_function('assemble_supplementary_figure_s2.R','zscore')
load_function('assemble_supplementary_figure_s2.R','heat_data')
for(label in c('baseline','corrected')){
 m<-readRDS(file.path(root,'results',paste0(label,'_expression_model.rds')))
 pc<-pca_data(m$expr,m$metadata,'GSE146638');cr<-corr_long(m$expr,m$metadata,'GSE146638')
 deg<-fread(file.path(root,'results',paste0(label,'_DEG.csv')))
 top<-head(deg[order(-abs(t))]$gene_symbol,30)
 heat<-heat_data(zscore(m$expr[top,,drop=FALSE]),m$metadata,'GSE146638')
 fwrite(pc,file.path(root,'results',paste0(label,'_PCA.tsv')),sep='\t')
 fwrite(data.table(component=c('PC1','PC2'),percent=attr(pc,'pct')),file.path(root,'results',paste0(label,'_PCA_variance.tsv')),sep='\t')
 fwrite(cr,file.path(root,'results',paste0(label,'_correlation.tsv')),sep='\t')
 fwrite(heat,file.path(root,'results',paste0(label,'_heatmap.tsv')),sep='\t')
 scores<-fread(file.path(root,'results',paste0(label,'_immune_scores.tsv')))
 scores[,z:=as.numeric(scale(score)),by=.(dataset,program_id)];scores[,sample_label:=paste(dataset,sample,sep=' | ')]
 fwrite(scores,file.path(root,'results',paste0(label,'_immune_scores_display.tsv')),sep='\t')
}
