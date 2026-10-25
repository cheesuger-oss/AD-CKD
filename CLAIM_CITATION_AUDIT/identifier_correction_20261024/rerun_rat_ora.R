options(stringsAsFactors=FALSE)
root<-normalizePath('CLAIM_CITATION_AUDIT/identifier_correction_20261024')
.libPaths(c(file.path(root,'Rlib'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(clusterProfiler);library(org.Hs.eg.db)})
oldmap<-fread(file.path(root,'../supplementary_source_inputs/AD-CKD/project/results/geo/ortholog_rat_human_1to1.csv'))
newmap<-fread(file.path(root,'ortholog_rat_human_1to1_corrected.csv'))
links<-fread(file.path(root,'frozen/data/geo/processed/kegg_pathway_gene_link.tsv'),header=FALSE,col.names=c('term','gene'))
names<-fread(file.path(root,'frozen/data/geo/processed/kegg_pathway_names.tsv'),header=FALSE,col.names=c('term','name'))
links[,`:=`(term=sub('^path:','',term),gene=sub('^hsa:','',gene))]
names[,`:=`(term=sub('^path:','',term),name=sub(' - Homo sapiens \\(human\\)$','',name))]
to_entrez<-function(s)unique(suppressMessages(bitr(s,fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Hs.eg.db))$ENTREZID)
summaries<-list()
for(label in c('baseline','corrected')){
 r<-fread(file.path(root,'results',paste0(label,'_all_genes.csv')))
 m<-if(label=='baseline')oldmap else newmap
 d<-r[adj.P.Val<.05 & abs(logFC)>.5]
 genes<-unique(m[source_symbol%in%d$gene_symbol,human_symbol]);universe<-unique(m[source_symbol%in%r$gene_symbol,human_symbol])
 go<-enrichGO(gene=genes,universe=universe,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont='BP',pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,minGSSize=10,maxGSSize=500,readable=FALSE)
 kg<-enricher(gene=to_entrez(genes),universe=to_entrez(universe),TERM2GENE=links,TERM2NAME=names,pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,minGSSize=10,maxGSSize=500)
 for(db in c('GO_BP','KEGG')){
  tab<-as.data.table(as.data.frame(if(db=='GO_BP')go else kg));tab[,`:=`(database=db,dataset='UREMIC')]
  fwrite(tab,file.path(root,'results',paste0(label,'_',db,'_UREMIC.csv')))
 }
 summaries[[label]]<-data.table(version=label,DEG_human=length(genes),universe_human=length(universe),GO_tested=nrow(as.data.frame(go)),KEGG_tested=nrow(as.data.frame(kg)))
}
fwrite(rbindlist(summaries),file.path(root,'results/ora_change_summary.csv'));print(rbindlist(summaries))
writeLines(capture.output(sessionInfo()),file.path(root,'results/ora_sessionInfo.txt'))
