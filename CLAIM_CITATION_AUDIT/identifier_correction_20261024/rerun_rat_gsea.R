options(stringsAsFactors=FALSE)
root<-normalizePath('CLAIM_CITATION_AUDIT/identifier_correction_20261024')
.libPaths(c(file.path(root,'Rlib'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(clusterProfiler);library(org.Hs.eg.db)})
# Use the original rank builder verbatim from the frozen generator, without executing its plotting or mouse analyses.
original<-parse(file.path(root,'../recovered/AD-CKD/project/scripts/phase1_gsea.R'))
rank_definition<-Filter(function(x)is.call(x)&&as.character(x[[1]])=='<-'&&identical(x[[2]],as.name('make_rank')),as.list(original))
stopifnot(length(rank_definition)==1);eval(rank_definition[[1]])
links<-fread(file.path(root,'frozen/data/geo/processed/kegg_pathway_gene_link.tsv'),header=FALSE,col.names=c('term','gene'))
names_k<-fread(file.path(root,'frozen/data/geo/processed/kegg_pathway_names.tsv'),header=FALSE,col.names=c('term','name'))
links[,`:=`(term=sub('^path:','',term),gene=sub('^hsa:','',gene))]
names_k[,`:=`(term=sub('^path:','',term),name=sub(' - Homo sapiens \\(human\\)$','',name))]
for(label in c('baseline','corrected')){
 set.seed(20260808)
 map_file<-if(label=='baseline')file.path(root,'../supplementary_source_inputs/AD-CKD/project/results/geo/ortholog_rat_human_1to1.csv') else file.path(root,'ortholog_rat_human_1to1_corrected.csv')
 gl<-make_rank(file.path(root,'results',paste0(label,'_all_genes.csv')),map_file)
 fwrite(data.table(ENTREZID=names(gl),t=as.numeric(gl)),file.path(root,'results',paste0(label,'_GSEA_rank.csv')))
 message(label,' GO GSEA')
 go<-gseGO(geneList=gl,OrgDb=org.Hs.eg.db,keyType='ENTREZID',ont='BP',by='fgsea',eps=0,minGSSize=10,maxGSSize=500,pvalueCutoff=1,pAdjustMethod='BH',verbose=FALSE)
 message(label,' KEGG GSEA')
 kg<-GSEA(geneList=gl,TERM2GENE=links,TERM2NAME=names_k,by='fgsea',eps=0,minGSSize=10,maxGSSize=500,pvalueCutoff=1,pAdjustMethod='BH',verbose=FALSE)
 for(db in c('GO_BP','KEGG')){
  tab<-as.data.table(as.data.frame(if(db=='GO_BP')go else kg));tab[,`:=`(database=db,dataset='Uremic')]
  fwrite(tab,file.path(root,'results',paste0(label,'_GSEA_',db,'.csv')))
 }
}
writeLines(capture.output(sessionInfo()),file.path(root,'results/gsea_sessionInfo.txt'))
