#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
.libPaths(c(file.path(root,"project","Rlib-4.5"), file.path(root,"renv","library","macos","R-4.5","aarch64-apple-darwin20"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(limma); library(ggplot2); library(patchwork) })

geo <- file.path(root,"project","results","geo")
fig <- file.path(geo,"figures","extended"); dir.create(fig,recursive=TRUE,showWarnings=FALSE)
pal <- c(Control="#0072B2",Uremic="#D55E00",AAD="#CC79A7")
theme_pub <- theme_classic(base_size=9)+theme(plot.title=element_text(face="bold",size=10),axis.title=element_text(face="bold"),plot.margin=margin(5,7,5,5))
wrap_text<-function(x,width=42)vapply(x,function(z)paste(strwrap(z,width=width),collapse="\n"),character(1))
save_pair <- function(p,stem,w=183,h=130){
 ggsave(file.path(fig,paste0(stem,".pdf")),p,width=w,height=h,units="mm",device="pdf",useDingbats=FALSE,bg="white")
 ggsave(file.path(fig,paste0(stem,".png")),p,width=w,height=h,units="mm",dpi=600,bg="white")
}
zscore <- function(x){ z<-t(scale(t(x))); z[!is.finite(z)]<-0; pmax(pmin(z,2.5),-2.5) }

# Exact two-set Venn counts (diagram is not area-proportional).
u <- fread(file.path(geo,"GSE146638_human_1to1_DEG.csv"))$human_gene
a <- fread(file.path(geo,"GSE325608_human_1to1_DEG.csv"))$human_gene
venn_counts <- data.frame(uremic_only=length(setdiff(u,a)),overlap=length(intersect(u,a)),aad_only=length(setdiff(a,u)),uremic_total=length(u),aad_total=length(a))
fwrite(venn_counts,file.path(geo,"venn_counts.tsv"),sep="\t")
theta<-seq(0,2*pi,length.out=500)
circles<-rbind(data.frame(x=-.75+1.35*cos(theta),y=1.35*sin(theta),set="Uremic"),data.frame(x=.75+1.35*cos(theta),y=1.35*sin(theta),set="AAD"))
pvenn<-ggplot(circles,aes(x,y,colour=set))+geom_path(linewidth=1.3)+
 annotate("text",x=-1.15,y=0,label=venn_counts$uremic_only,size=6)+annotate("text",x=0,y=0,label=venn_counts$overlap,size=6,fontface="bold")+
 annotate("text",x=1.15,y=0,label=venn_counts$aad_only,size=6)+annotate("text",x=-1.25,y=1.55,label=paste0("Uremic DEG\n(n=",length(u),")"),size=4,fontface="bold")+
 annotate("text",x=1.25,y=1.55,label=paste0("AAD DEG\n(n=",length(a),")"),size=4,fontface="bold")+
 scale_colour_manual(values=c(Uremic="#D55E00",AAD="#0072B2"),guide="none")+coord_fixed(xlim=c(-2.5,2.5),ylim=c(-1.7,2.0),clip="off")+
 labs(title="Human 1:1 ortholog DEG overlap",subtitle="Exact counts; circle areas are schematic",x=NULL,y=NULL)+theme_void(base_size=10)+theme(plot.title=element_text(face="bold",hjust=.5),plot.subtitle=element_text(hjust=.5))
save_pair(pvenn,"Figure2A_DEG_Venn",100,80)

# Reconstruct normalized matrices exactly as Phase 1 and draw top-30 DEG heatmaps.
m1<-fread(file.path(geo,"GSE146638_sample_metadata.tsv")); x1<-fread(file.path(root,"project","data","geo","processed","GSE146638_expression.tsv.gz"))
ids1<-x1[[1]]; mat1<-rowsum(as.matrix(x1[,-1]),group=ids1,reorder=FALSE,na.rm=TRUE); m1<-m1[include_primary=="yes"]
log1<-log2(mat1[,m1$expression_column,drop=FALSE]+1); d1<-fread(file.path(geo,"GSE146638_DEG.csv")); top1<-head(d1[order(-abs(t))]$gene_symbol,30); h1<-zscore(log1[top1,,drop=FALSE])

m2<-fread(file.path(geo,"GSE325608_sample_metadata.tsv")); x2<-fread(file.path(root,"project","data","geo","GSE325608_raw_counts.txt.gz"))
ids2<-x2[[1]]; mat2<-as.matrix(x2[,-1]); rownames(mat2)<-ids2; m2<-m2[include_primary=="yes"]
keep2<-rowSums(mat2[,m2$expression_column,drop=FALSE]>=10)>=5; m2$sex<-factor(m2$sex,levels=c("Male","Female")); m2$primary_group<-factor(m2$primary_group,levels=c("Control","AAD"))
design2<-model.matrix(~sex+primary_group,data=m2); v2<-voom(mat2[keep2,m2$expression_column,drop=FALSE],design2,plot=FALSE)
d2<-fread(file.path(geo,"GSE325608_DEG.csv")); top2<-head(d2[order(-abs(t))]$gene_symbol,30); h2<-zscore(v2$E[top2,,drop=FALSE])

heat_plot<-function(mat,meta,title){
 dd<-as.data.table(as.table(mat)); names(dd)<-c("gene","sample","z"); dd$gene<-factor(dd$gene,levels=rev(rownames(mat)))
 sm<-data.frame(sample=meta$expression_column,group=meta$primary_group); dd<-merge(dd,sm,by="sample",sort=FALSE); dd$sample<-factor(dd$sample,levels=meta$expression_column)
 ggplot(dd,aes(sample,gene,fill=z))+geom_tile()+scale_fill_gradient2(low="#0072B2",mid="white",high="#D55E00",midpoint=0,limits=c(-2.5,2.5),oob=scales::squish)+
  facet_grid(~group,scales="free_x",space="free_x")+labs(title=title,x="Samples",y=NULL,fill="Row z-score")+theme_pub+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),strip.background=element_blank(),strip.text=element_text(face="bold"),panel.spacing.x=unit(1,"mm"))
}
pheat1<-heat_plot(h1,m1,"GSE146638 top 30 DEGs")
pheat2<-heat_plot(h2,m2,"GSE325608 top 30 DEGs (sex-adjusted model)")
save_pair(pheat1,"GSE146638_top_DEG_heatmap",120,125); save_pair(pheat2,"GSE325608_top_DEG_heatmap",150,125)

# Complete GO/KEGG bubble overview: top five terms per dataset/database by BH-adjusted P.
read_enr<-function(path,label,db){x<-fread(path);x[,`:=`(dataset=label,database=db)];x}
enr<-rbindlist(list(read_enr(file.path(geo,"GOBP_UREMIC.csv"),"Uremic","GO BP"),read_enr(file.path(geo,"GOBP_AD.csv"),"AAD","GO BP"),read_enr(file.path(geo,"KEGG_UREMIC.csv"),"Uremic","KEGG"),read_enr(file.path(geo,"KEGG_AD.csv"),"AAD","KEGG")),fill=TRUE)
enr<-enr[is.finite(p.adjust)&p.adjust<.05]
top<-enr[order(p.adjust,-Count),head(.SD,5),by=.(dataset,database)]
fwrite(top[, .(dataset, database, ID, Description, Count, FoldEnrichment, pvalue, p.adjust)],
       file.path(fig,"Supplementary_Figure_S3a_ORA_bubble_source_data.tsv"),sep="\t")
top[,term:=wrap_text(Description,42)]; top[,term_panel:=paste(term,dataset,database,sep="___")]
top[,term_panel:=factor(term_panel,levels=rev(unique(term_panel)))]
top[,neglogq:=-log10(pmax(p.adjust,.Machine$double.xmin))]
pbubble<-ggplot(top,aes(FoldEnrichment,term_panel,size=Count,colour=neglogq))+geom_point(alpha=.9)+facet_grid(database~dataset,scales="free_y",space="free_y")+
 scale_y_discrete(labels=function(x)sub("___.*$","",x))+scale_colour_viridis_c(option="C")+scale_size_continuous(range=c(2.5,7))+
 labs(title="Animal-model GO/KEGG over-representation",subtitle="Five smallest BH-adjusted P values per model and database",x="Fold enrichment",y=NULL,colour=expression(-log[10](BH~q)),size="Genes")+
 theme_pub+theme(strip.background=element_blank(),strip.text=element_text(face="bold"),legend.position="right")
save_pair(pbubble,"Figure2B_GO_KEGG_bubble",183,170)

# GO bar overview uses fold enrichment, zero baseline and exact BH labels.
go<-rbindlist(list(read_enr(file.path(geo,"GOBP_UREMIC.csv"),"Uremic","GO BP"),read_enr(file.path(geo,"GOBP_AD.csv"),"AAD","GO BP")),fill=TRUE)
gobar<-go[is.finite(p.adjust)&p.adjust<.05][order(p.adjust,-Count),head(.SD,10),by=dataset]
fwrite(gobar[, .(dataset, database, ID, Description, Count, FoldEnrichment, pvalue, p.adjust)],
       file.path(fig,"Supplementary_Figure_S3b_GO_bar_source_data.tsv"),sep="\t")
gobar[,term:=wrap_text(Description,48)]; gobar[,term_panel:=factor(paste(term,dataset,sep="___"),levels=rev(unique(paste(term,dataset,sep="___"))))]
pbar<-ggplot(gobar,aes(FoldEnrichment,term_panel,fill=dataset))+geom_col(width=.72,colour="black",linewidth=.15)+facet_grid(dataset~.,scales="free_y",space="free_y")+
 scale_y_discrete(labels=function(x)sub("___.*$","",x))+scale_fill_manual(values=c(Uremic="#D55E00",AAD="#0072B2"),guide="none")+
 scale_x_continuous(expand=expansion(mult=c(0,.04)))+labs(title="GO Biological Process overview",subtitle="Ten smallest BH-adjusted P values per dataset",x="Fold enrichment",y=NULL)+theme_pub+
 theme(strip.background=element_blank(),strip.text=element_text(face="bold"))
save_pair(pbar,"Supplementary_GO_barplot",160,180)

writeLines(capture.output(sessionInfo()),file.path(fig,"extended_figure_sessionInfo.txt"))
cat("EXTENDED_FIGURES_COMPLETE\n")
