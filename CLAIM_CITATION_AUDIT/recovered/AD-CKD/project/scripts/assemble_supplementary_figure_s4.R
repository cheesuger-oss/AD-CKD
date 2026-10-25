#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
.libPaths(c(file.path(root,"project","Rlib-4.5"), file.path(root,"renv","library","macos","R-4.5","aarch64-apple-darwin20"), "/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library", .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma); library(AnnotationDbi)
  library(org.Hs.eg.db); library(hgu133a2.db); library(ggplot2); library(patchwork)
})
set.seed(20260810)
hv <- file.path(root,"project","results","geo","human_validation")
out <- file.path(root,"project","results","figures"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
theme_pub <- theme_classic(base_size=7.2)+theme(plot.title=element_text(face="bold",size=8.5),plot.subtitle=element_text(size=6.4),axis.title=element_text(face="bold"),plot.margin=margin(4,5,4,4))
zscore <- function(x){z<-t(scale(t(x)));z[!is.finite(z)]<-0;pmax(pmin(z,2.5),-2.5)}
collapse_by_id <- function(mat,ids,rule=c("sum","maxvar")){rule<-match.arg(rule);keep<-!is.na(ids)&nzchar(ids);mat<-mat[keep,,drop=FALSE];ids<-as.character(ids[keep]);if(rule=="sum")return(rowsum(mat,ids,reorder=FALSE));v<-apply(mat,1,var);o<-order(ids,-v,rownames(mat));f<-!duplicated(ids[o]);ans<-mat[o[f],,drop=FALSE];rownames(ans)<-ids[o[f]];ans}
read_series_matrix <- function(path){x<-readLines(gzfile(path),warn=FALSE);a<-grep("^!series_matrix_table_begin",x)+1L;b<-grep("^!series_matrix_table_end",x)-1L;con<-textConnection(x[a:b]);on.exit(close(con));read.delim(con,check.names=FALSE,quote="\"",row.names=1)}

# Reconstruct the frozen analysed matrices without changing the fitted contrasts.
downloads <- file.path(path.expand("~"),"Downloads")
aad <- read.delim(gzfile(file.path(downloads,"GSE153434_all.counts.txt.gz")),check.names=FALSE)
sample_cols <- c(paste0("Normal",c(1,10,2:9)),paste0("Patient",c(10,1:9)))
cnt <- as.matrix(aad[,sample_cols]); storage.mode(cnt)<-"numeric"
ent <- AnnotationDbi::mapIds(org.Hs.eg.db,keys=aad$Symbol,keytype="SYMBOL",column="ENTREZID",multiVals="first")
cnt <- collapse_by_id(cnt,ent,"sum"); group1<-factor(c(rep("Control",10),rep("AAD",10)),levels=c("Control","AAD"))
y<-DGEList(cnt);keep<-filterByExpr(y,group=group1);y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y);v<-voom(y,model.matrix(~group1),plot=FALSE)
expr1<-v$E; de1<-fread(file.path(hv,"GSE153434_all_genes.csv")); meta1<-data.frame(sample=colnames(expr1),group=group1,display=c(paste0("C",sprintf("%02d",1:10)),paste0("AAD",sprintf("%02d",1:10))))

mic<-read_series_matrix(file.path(downloads,"GSE38750_series_matrix.txt.gz")); pm<-AnnotationDbi::select(hgu133a2.db,keys=rownames(mic),keytype="PROBEID",columns=c("ENTREZID","SYMBOL"));pm<-pm[!duplicated(pm$PROBEID),]
expr2<-collapse_by_id(as.matrix(mic),pm$ENTREZID[match(rownames(mic),pm$PROBEID)],"maxvar");group2<-factor(c(rep("Control",19),rep("Uremic",15)),levels=c("Control","Uremic"));de2<-fread(file.path(hv,"GSE38750_all_genes.csv"));meta2<-data.frame(sample=colnames(expr2),group=group2,display=c(paste0("C",sprintf("%02d",1:19)),paste0("U",sprintf("%02d",1:15))))

make_page <- function(expr,de,meta,dataset,confounded=FALSE){
  vars<-apply(expr,1,var);use<-names(sort(vars,decreasing=TRUE))[seq_len(min(1000,length(vars)))];pc<-prcomp(t(expr[use,,drop=FALSE]));ve<-100*pc$sdev^2/sum(pc$sdev^2)
  pd<-data.frame(dataset=dataset,sample=rownames(pc$x),display=meta$display,group=meta$group,PC1=pc$x[,1],PC2=pc$x[,2],pca_input_n_genes=length(use),variance_PC1_pct=ve[1],variance_PC2_pct=ve[2])
  pca_sub<-if(confounded)"1,000 most variable genes; uremic iliac versus control renal artery" else "1,000 most variable analysed genes; no PCA-based exclusion"
  p1<-ggplot(pd,aes(PC1,PC2,colour=group,shape=group))+geom_point(size=2)+scale_colour_manual(values=c(Control="#0072B2",AAD="#D55E00",Uremic="#CC79A7"))+labs(tag="A",title="PCA",subtitle=pca_sub,x=sprintf("PC1 (%.1f%%)",ve[1]),y=sprintf("PC2 (%.1f%%)",ve[2]),colour=NULL,shape=NULL)+theme_pub+theme(legend.position="bottom")
  co<-cor(expr,method="pearson");cd<-as.data.frame(as.table(co));names(cd)<-c("sample_x","sample_y","correlation");cd$display_x<-meta$display[match(cd$sample_x,meta$sample)];cd$display_y<-meta$display[match(cd$sample_y,meta$sample)];cd$group_x<-meta$group[match(cd$sample_x,meta$sample)];cd$group_y<-meta$group[match(cd$sample_y,meta$sample)];cd$display_x<-factor(cd$display_x,levels=meta$display);cd$display_y<-factor(cd$display_y,levels=rev(meta$display))
  p2<-ggplot(cd,aes(display_x,display_y,fill=correlation))+geom_tile()+scale_fill_viridis_c(option="C",limits=c(min(cd$correlation),1),name="Pearson r")+labs(tag="B",title="Sample correlation",subtitle="Metadata order; display only",x=NULL,y=NULL)+theme_pub+theme(axis.text.x=element_text(angle=90,hjust=1,size=5.1),axis.text.y=element_text(size=5.1),axis.ticks=element_blank())
  de$status<-"Not DEG";de$status[de$FDR<.05&de$logFC>.5]<-"Up";de$status[de$FDR<.05&de$logFC < -0.5]<-"Down";de$status<-factor(de$status,levels=c("Down","Not DEG","Up"));de$neglogq<--log10(pmax(de$FDR,.Machine$double.xmin))
  p3<-ggplot(de,aes(logFC,neglogq,colour=status))+geom_point(size=.45,alpha=.55)+geom_vline(xintercept=c(-.5,.5),linetype="dashed",colour="grey55")+geom_hline(yintercept=-log10(.05),linetype="dashed",colour="grey55")+scale_colour_manual(values=c(Down="#0072B2","Not DEG"="grey78",Up="#D55E00"))+labs(tag="C",title="Differential-expression overview",subtitle="All tested genes; BH q<0.05 and abs(log2FC)>0.5",x="log2 fold change",y=expression(-log[10](BH~q)),colour=NULL)+theme_pub+theme(legend.position="bottom")
  top<-as.character(head(de[order(-abs(t))]$ENTREZID,40));hm<-zscore(expr[intersect(top,rownames(expr)),,drop=FALSE]);hd<-as.data.frame(as.table(hm));names(hd)<-c("gene_id","sample","row_z");sym<-setNames(de$Symbol,as.character(de$ENTREZID));gid<-as.character(hd$gene_id);rowid<-rownames(hm);display_gene<-ifelse(is.na(sym[gid])|sym[gid]=="",gid,sym[gid]);row_gene<-ifelse(is.na(sym[rowid])|sym[rowid]=="",rowid,sym[rowid]);hd$gene<-factor(display_gene,levels=rev(unique(row_gene)));hd$display<-factor(meta$display[match(hd$sample,meta$sample)],levels=meta$display);hd$group<-meta$group[match(hd$sample,meta$sample)]
  sub4<-"Top 40 by absolute moderated t; selected from this contrast"
  p4<-ggplot(hd,aes(display,gene,fill=row_z))+geom_tile()+scale_fill_gradient2(low="#0072B2",mid="white",high="#D55E00",midpoint=0,limits=c(-2.5,2.5),oob=scales::squish,name="Row z-score")+facet_grid(~group,scales="free_x",space="free_x")+labs(tag="D",title="Selected-gene heatmap",subtitle=sub4,x="Samples",y=NULL)+theme_pub+theme(axis.text.x=element_blank(),axis.text.y=element_text(size=5.0),axis.ticks.x=element_blank(),strip.background=element_blank(),strip.text=element_text(face="bold"))
  page_sub<-if(confounded)"Uremic iliac artery versus control renal artery; disease and arterial bed are completely confounded" else "Ascending-aorta tissue: 10 controls and 10 aortic-dissection cases"
  page<-((p1|p2)/(p3|p4))+plot_annotation(title=dataset,subtitle=page_sub,theme=theme(plot.title=element_text(face="bold",size=10),plot.subtitle=element_text(size=7.2,face=if(confounded)"bold" else "plain",colour=if(confounded)"#A33A2B" else "black",margin=margin(b=4)))) & theme(plot.tag=element_text(face="bold",size=9),plot.tag.position=c(0,1))
  list(plot=page,pca=pd,cor=cd,volcano=de[,.(ENTREZID,Symbol,logFC,t,P.Value,FDR,status)],heat=hd)
}
r1<-make_page(expr1,de1,meta1,"GSE153434: human ascending aorta",FALSE);r2<-make_page(expr2,de2,meta2,"GSE38750: human vascular tissue",TRUE)
save_page<-function(obj,stem){ggsave(file.path(out,paste0(stem,".pdf")),obj$plot,width=183,height=170,units="mm",device="pdf",useDingbats=FALSE,bg="white");ggsave(file.path(out,paste0(stem,".png")),obj$plot,width=183,height=170,units="mm",dpi=600,bg="white")}
save_page(r1,"Supplementary_Figure_S4a_human_AAD_expression_QC");save_page(r2,"Supplementary_Figure_S4b_human_uremic_expression_QC")
write_obj<-function(obj,prefix){for(nm in c("pca","cor","volcano","heat"))fwrite(obj[[nm]],file.path(out,paste0(prefix,"_source_data_",nm,".tsv")),sep="\t")}
write_obj(r1,"Supplementary_Figure_S4a");write_obj(r2,"Supplementary_Figure_S4b")
idx<-data.frame(page=c("S4a","S4b"),visual_file=c("Supplementary_Figure_S4a_human_AAD_expression_QC.pdf","Supplementary_Figure_S4b_human_uremic_expression_QC.pdf"),source_prefix=c("Supplementary_Figure_S4a_source_data_","Supplementary_Figure_S4b_source_data_"),claim_ids="C10-C11 context",boundary=c("Same arterial site; descriptive tissue expression","Disease status completely confounded with arterial site"))
fwrite(idx,file.path(out,"Supplementary_Figure_S4_source_data.tsv"),sep="\t");writeLines(capture.output(sessionInfo()),file.path(out,"Supplementary_Figure_S4_sessionInfo.txt"));cat("SUPPLEMENTARY_FIGURE_S4_COMPLETE\n")
