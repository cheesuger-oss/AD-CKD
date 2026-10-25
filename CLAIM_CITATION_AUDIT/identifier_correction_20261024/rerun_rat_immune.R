options(stringsAsFactors=FALSE)
root<-normalizePath('CLAIM_CITATION_AUDIT/identifier_correction_20261024')
suppressPackageStartupMessages({library(limma);library(edgeR);library(data.table)})
original<-parse(file.path(root,'../recovered/AD-CKD/project/scripts/phase4_immune_ml_analysis.R'))
for(name in c('collapse_symbol','rank_scores')){
 definition<-Filter(function(x)is.call(x)&&as.character(x[[1]])=='<-'&&identical(x[[2]],as.name(name)),as.list(original))
 stopifnot(length(definition)==1);eval(definition[[1]])
}
membership<-read.delim(file.path(root,'frozen/results/phase4/immune_program_membership.tsv'),check.names=FALSE)
rat_meta<-read.delim(file.path(root,'frozen/results/geo/GSE146638_sample_metadata.tsv'),check.names=FALSE)
rat_meta<-rat_meta[rat_meta$include_primary=='yes',]
for(label in c('baseline','corrected')){
 path<-if(label=='baseline')file.path(root,'../identifier_source_audit/GSE146638_expression.tsv.gz') else file.path(root,'GSE146638_expression_corrected.tsv.gz')
 map_path<-if(label=='baseline')file.path(root,'../supplementary_source_inputs/AD-CKD/project/results/geo/ortholog_rat_human_1to1.csv') else file.path(root,'ortholog_rat_human_1to1_corrected.csv')
 rat<-read.delim(gzfile(path),check.names=FALSE)
 rat_mat<-as.matrix(rat[match(unique(rat$gene_symbol),rat$gene_symbol),rat_meta$expression_column]);storage.mode(rat_mat)<-'numeric';rownames(rat_mat)<-unique(rat$gene_symbol)
 rat_map<-read.csv(map_path);mi<-match(rownames(rat_mat),rat_map$source_symbol);rat_h<-collapse_symbol(rat_mat,rat_map$human_symbol[mi])
 y_rat<-DGEList(rat_h);y_rat<-calcNormFactors(y_rat);rat_expr<-cpm(y_rat,log=TRUE,prior.count=2)
 group<-factor(ifelse(rat_meta$primary_group=='Uremic','Case','Control'),levels=c('Control','Case'))
 result<-rank_scores(rat_expr,membership,'GSE146638',group)
 fwrite(result$contrast,file.path(root,'results',paste0(label,'_immune_contrasts.tsv')),sep='\t')
 fwrite(result$long,file.path(root,'results',paste0(label,'_immune_scores.tsv')),sep='\t')
 print(result$contrast[,c('program_id','observed_gene_n','standardized_effect','q_value')])
}
writeLines(capture.output(sessionInfo()),file.path(root,'results/immune_sessionInfo.txt'))
