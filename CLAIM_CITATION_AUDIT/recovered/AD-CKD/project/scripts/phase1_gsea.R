#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
args<-commandArgs(trailingOnly=TRUE); root<-if(length(args))normalizePath(args[[1]]) else normalizePath(".")
.libPaths(c(file.path(root,"project","Rlib-4.5"),file.path(root,"renv","library","macos","R-4.5","aarch64-apple-darwin20"),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(clusterProfiler);library(org.Hs.eg.db);library(enrichplot);library(ggplot2);library(patchwork)})
set.seed(20260808)
geo<-file.path(root,"project","results","geo"); fig<-file.path(geo,"figures","extended");dir.create(fig,recursive=TRUE,showWarnings=FALSE)

make_rank<-function(result_file,map_file){
 r<-fread(result_file);m<-fread(map_file);z<-merge(r,m,by.x="gene_symbol",by.y="source_symbol");z[,abs_t:=abs(t)]
 setorder(z,-abs_t);z<-z[!duplicated(human_symbol)]
 conv<-suppressMessages(bitr(z$human_symbol,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db))
 z<-merge(z,conv,by.x="human_symbol",by.y="SYMBOL");setorder(z,-abs_t);z<-z[!duplicated(ENTREZID)]
 ranks<-z$t;names(ranks)<-z$ENTREZID;sort(ranks,decreasing=TRUE)
}
ranks<-list(
 Uremic=make_rank(file.path(geo,"GSE146638_all_genes.csv"),file.path(geo,"ortholog_rat_human_1to1.csv")),
 AAD=make_rank(file.path(geo,"GSE325608_all_genes.csv"),file.path(geo,"ortholog_mouse_human_1to1.csv"))
)

link_candidates<-c(file.path(root,"project","data","geo","processed","kegg_pathway_gene_link.tsv"),file.path(root,"Downloads","KEGG_pathway-gene_link.tsv"),"inputs/KEGG_pathway-gene_link.tsv")
name_candidates<-c(file.path(root,"project","data","geo","processed","kegg_pathway_names.tsv"),file.path(root,"Downloads","KEGG _pahtway_names.tsv"),"inputs/KEGG _pahtway_names.tsv")
link_file<-link_candidates[file.exists(link_candidates)][1];name_file<-name_candidates[file.exists(name_candidates)][1]
if(is.na(link_file)||is.na(name_file))stop("Frozen KEGG TERM2GENE/TERM2NAME files unavailable")
links<-fread(link_file,header=FALSE,col.names=c("term","gene"));names_k<-fread(name_file,header=FALSE,col.names=c("term","name"))
links[,`:=`(term=sub("^path:","",term),gene=sub("^hsa:","",gene))];names_k[,`:=`(term=sub("^path:","",term),name=sub(" - Homo sapiens \\(human\\)$","",name))]

results<-list(); plots<-list()
for(ds in names(ranks)){
 gl<-ranks[[ds]]
 go<-gseGO(geneList=gl,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont="BP",by="fgsea",eps=0,minGSSize=10,maxGSSize=500,pvalueCutoff=1,pAdjustMethod="BH",verbose=FALSE)
 kg<-GSEA(geneList=gl,TERM2GENE=links,TERM2NAME=names_k,by="fgsea",eps=0,minGSSize=10,maxGSSize=500,pvalueCutoff=1,pAdjustMethod="BH",verbose=FALSE)
 go_df<-as.data.table(as.data.frame(go));kg_df<-as.data.table(as.data.frame(kg));if(nrow(go_df))go_df[,`:=`(dataset=ds,database="GO_BP")];if(nrow(kg_df))kg_df[,`:=`(dataset=ds,database="KEGG")]
 fwrite(go_df,file.path(geo,paste0("GSEA_GOBP_",toupper(ds),".csv")));fwrite(kg_df,file.path(geo,paste0("GSEA_KEGG_",toupper(ds),".csv")))
 results[[paste(ds,"GO")]]<-go_df;results[[paste(ds,"KEGG")]]<-kg_df
 for(obj_name in c("GO","KEGG")){
   obj<-if(obj_name=="GO")go else kg;df<-as.data.frame(obj)
   if(nrow(df)){
     pick<-df$ID[order(df$p.adjust,-abs(df$NES))][1]
     chosen<-df[match(pick,df$ID),]
     fwrite(as.data.table(chosen),file.path(fig,paste0("Supplementary_Figure_S3_",toupper(ds),"_",obj_name,"_selected_source_data.tsv")),sep="\t")
     page_label<-c("AAD_GO"="S3c","AAD_KEGG"="S3d","Uremic_GO"="S3e","Uremic_KEGG"="S3f")[[paste(ds,obj_name,sep="_")]]
     boundary<-if(ds=="Uremic"&&obj_name=="KEGG")"\nKEGG label only; not evidence of infection" else ""
     plot_title<-paste0(page_label," | ",ds," ",obj_name,"\n",chosen$Description,"; NES=",sprintf("%.2f",chosen$NES),"; BH q=",format(chosen$p.adjust,digits=2,scientific=TRUE),boundary)
     p<-gseaplot2(obj,geneSetID=pick,title=plot_title,base_size=7.5,subplots=1:3,pvalue_table=FALSE,color="#0072B2")
     ggsave(file.path(fig,paste0("GSEA_",toupper(ds),"_",obj_name,".pdf")),p,width=183,height=125,units="mm",device="pdf",useDingbats=FALSE,bg="white",limitsize=FALSE)
     ggsave(file.path(fig,paste0("GSEA_",toupper(ds),"_",obj_name,".png")),p,width=183,height=125,units="mm",dpi=600,bg="white",limitsize=FALSE)
   }
 }
}
all<-rbindlist(results,fill=TRUE);fwrite(all,file.path(geo,"GSEA_all_results.csv"))
summary<-all[,.(tested_terms=.N,BH_lt_0.05=sum(p.adjust<.05,na.rm=TRUE),positive_NES=sum(NES>0&p.adjust<.05,na.rm=TRUE),negative_NES=sum(NES<0&p.adjust<.05,na.rm=TRUE)),by=.(dataset,database)]
fwrite(summary,file.path(geo,"GSEA_summary.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(geo,"GSEA_sessionInfo.txt"));print(summary)
