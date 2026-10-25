#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
.libPaths(c("project/Rlib-4.5", "/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library", .libPaths()))
suppressPackageStartupMessages({
  library(edgeR); library(limma); library(AnnotationDbi); library(org.Hs.eg.db)
  library(hgu133a2.db); library(ggplot2); library(pheatmap)
})

set.seed(20260809)
out <- "project/results/geo/human_validation"
fig <- file.path(out, "figures")
dir.create(fig, recursive=TRUE, showWarnings=FALSE)

go_ids <- c("GO:0030198","GO:0001503","GO:0006939","GO:0048659","GO:0048660","GO:0048661")
go_names <- c("extracellular matrix organization","ossification","smooth muscle contraction",
              "smooth muscle cell proliferation","regulation of smooth muscle cell proliferation",
              "positive regulation of smooth muscle cell proliferation")
names(go_names) <- go_ids
go_map <- AnnotationDbi::select(org.Hs.eg.db, keys=go_ids, keytype="GOALL",
                                columns=c("ENTREZID","ONTOLOGYALL"))
go_map <- unique(go_map[go_map$ONTOLOGYALL == "BP" & !is.na(go_map$ENTREZID), c("GOALL","ENTREZID")])
gene_sets <- split(go_map$ENTREZID, go_map$GOALL)[go_ids]
write.table(go_map, file.path(out,"frozen_six_pathway_membership.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

read_series_matrix <- function(path) {
  x <- readLines(gzfile(path), warn=FALSE)
  a <- grep("^!series_matrix_table_begin", x) + 1L
  b <- grep("^!series_matrix_table_end", x) - 1L
  con <- textConnection(x[a:b]); on.exit(close(con))
  read.delim(con, check.names=FALSE, quote="\"", row.names=1)
}

collapse_by_id <- function(mat, ids, rule=c("sum","maxvar")) {
  rule <- match.arg(rule); keep <- !is.na(ids) & nzchar(ids)
  mat <- mat[keep,,drop=FALSE]; ids <- as.character(ids[keep])
  if (rule == "sum") return(rowsum(mat, ids, reorder=FALSE))
  v <- apply(mat,1,var); ord <- order(ids, -v, rownames(mat)); first <- !duplicated(ids[ord])
  ans <- mat[ord[first],,drop=FALSE]; rownames(ans) <- ids[ord[first]]; ans
}

save_gg <- function(p, stem, width=7.1, height=5.0) {
  ggsave(file.path(fig,paste0(stem,".png")), p, width=width, height=height, dpi=300, bg="white")
  pdf(file.path(fig,paste0(stem,".pdf")), width=width, height=height, useDingbats=FALSE)
  print(p); dev.off()
}

plot_qc <- function(expr, meta, prefix) {
  vars <- apply(expr,1,var); use <- names(sort(vars,decreasing=TRUE))[seq_len(min(1000,length(vars)))]
  pc <- prcomp(t(expr[use,,drop=FALSE]), scale.=FALSE)
  pd <- data.frame(sample=rownames(pc$x), PC1=pc$x[,1], PC2=pc$x[,2], group=meta$group)
  ve <- 100*pc$sdev^2/sum(pc$sdev^2)
  write.csv(pd,file.path(out,paste0(prefix,"_PCA_coordinates.csv")),row.names=FALSE)
  p <- ggplot(pd,aes(PC1,PC2,color=group,shape=group,label=sample))+geom_point(size=3)+
    scale_color_manual(values=c(Control="#0072B2",Case="#D55E00"))+
    labs(x=sprintf("PC1 (%.1f%%)",ve[1]),y=sprintf("PC2 (%.1f%%)",ve[2]),title=paste(prefix,"PCA"))+
    theme_classic(base_size=10)+theme(legend.position="top")
  save_gg(p,paste0(prefix,"_PCA"),6.0,4.7)
  ann <- data.frame(Group=meta$group); rownames(ann) <- colnames(expr)
  cols <- list(Group=c(Control="#0072B2",Case="#D55E00"))
  png(file.path(fig,paste0(prefix,"_sample_correlation.png")),1800,1600,res=250)
  pheatmap(cor(expr,method="pearson"),annotation_col=ann,annotation_row=ann,
           annotation_colors=cols,border_color=NA,main=paste(prefix,"sample correlation"),fontsize=7)
  dev.off()
  pdf(file.path(fig,paste0(prefix,"_sample_correlation.pdf")),7.2,6.4)
  pheatmap(cor(expr,method="pearson"),annotation_col=ann,annotation_row=ann,
           annotation_colors=cols,border_color=NA,main=paste(prefix,"sample correlation"),fontsize=7)
  dev.off()
  pd
}

plot_de <- function(tab, expr, meta, prefix) {
  tab$status <- ifelse(tab$FDR<0.05 & tab$logFC>0.5,"Up",ifelse(tab$FDR<0.05 & tab$logFC< -0.5,"Down","Not significant"))
  p <- ggplot(tab,aes(logFC,-log10(pmax(P.Value,1e-300)),color=status))+
    geom_point(alpha=.55,size=1)+scale_color_manual(values=c(Down="#0072B2","Not significant"="#B7B7B7",Up="#D55E00"))+
    geom_vline(xintercept=c(-.5,.5),linetype=2,linewidth=.3)+geom_hline(yintercept=-log10(.05),linetype=3,linewidth=.3)+
    labs(x="log2 fold change",y="-log10(P)",title=paste(prefix,"differential expression"),color=NULL)+
    theme_classic(base_size=10)+theme(legend.position="top")
  save_gg(p,paste0(prefix,"_volcano"),6.2,4.8)
  top <- head(tab$ENTREZID[order(abs(tab$t),decreasing=TRUE)],40)
  z <- t(scale(t(expr[top,,drop=FALSE]))); z[z>3]<-3; z[z< -3]<- -3
  ann <- data.frame(Group=meta$group); rownames(ann)<-colnames(expr)
  cols <- list(Group=c(Control="#0072B2",Case="#D55E00"))
  png(file.path(fig,paste0(prefix,"_top40_heatmap.png")),1700,1900,res=250)
  pheatmap(z,annotation_col=ann,annotation_colors=cols,show_colnames=FALSE,border_color=NA,
           main=paste(prefix,"top 40 genes by |moderated t|"),fontsize_row=6)
  dev.off()
}

gsea_es <- function(stats, members) {
  hit <- names(stats) %in% members; nh <- sum(hit); n <- length(stats)
  if(nh<5 || nh>=n) return(NA_real_)
  w <- abs(stats); ph <- ifelse(hit,w/sum(w[hit]),0); pm <- ifelse(!hit,1/(n-nh),0)
  run <- cumsum(ph-pm); mx <- max(run); mn <- min(run)
  if(abs(mx)>=abs(mn)) mx else mn
}

run_gsea <- function(stats, prefix, nperm=5000L) {
  stats <- sort(stats[is.finite(stats)],decreasing=TRUE); ids <- names(stats); n <- length(stats)
  rows <- vector("list",length(go_ids)); curves <- vector("list",length(go_ids))
  for(i in seq_along(go_ids)) {
    gs <- intersect(gene_sets[[go_ids[i]]],ids); es <- gsea_es(stats,gs); k <- length(gs)
    null <- replicate(nperm,gsea_es(stats,sample(ids,k)))
    den <- if(es>=0) mean(null[null>=0]) else abs(mean(null[null<0])); nes <- es/den
    p <- if(es>=0) (1+sum(null>=es))/(1+sum(null>=0)) else (1+sum(null<=es))/(1+sum(null<0))
    hit <- ids %in% gs; w <- abs(stats); running <- cumsum(ifelse(hit,w/sum(w[hit]),1/(n-k)*-1))
    peak <- if(es>=0) which.max(running) else which.min(running)
    leading <- if(es>=0) ids[which(hit & seq_along(ids)<=peak)] else ids[which(hit & seq_along(ids)>=peak)]
    rows[[i]] <- data.frame(dataset=prefix,GO_ID=go_ids[i],pathway=go_names[i],set_size=k,ES=es,NES=nes,P_value=p,
                            leading_edge=paste(leading,collapse="/"))
    curves[[i]] <- data.frame(rank=seq_len(n),running=running,GO_ID=go_ids[i],pathway=go_names[i])
  }
  ans <- do.call(rbind,rows); ans$BH_q <- p.adjust(ans$P_value,"BH")
  write.csv(ans,file.path(out,paste0(prefix,"_six_pathway_GSEA.csv")),row.names=FALSE)
  cv <- do.call(rbind,curves); lab <- paste0(ans$pathway,"\nNES=",sprintf("%.2f",ans$NES),", q=",sprintf("%.3g",ans$BH_q)); names(lab)<-ans$pathway
  cv$facet <- factor(lab[cv$pathway],levels=lab)
  p <- ggplot(cv,aes(rank,running))+geom_hline(yintercept=0,color="grey70",linewidth=.25)+
    geom_line(color="#009E73",linewidth=.45)+facet_wrap(~facet,scales="free_y",ncol=2)+
    labs(x="Rank in moderated-t gene list",y="Running enrichment score",title=paste(prefix,"frozen six-pathway GSEA"))+
    theme_classic(base_size=8)+theme(strip.text=element_text(size=7),plot.title=element_text(size=10))
  save_gg(p,paste0(prefix,"_six_pathway_GSEA"),7.2,7.0)
  ans
}

# GSE153434: raw integer counts, gene-level sum before testing.
aad <- read.delim(gzfile("inputs/GSE153434_all.counts.txt.gz"),check.names=FALSE)
sample_cols <- c(paste0("Normal",c(1,10,2:9)),paste0("Patient",c(10,1:9)))
cnt <- as.matrix(aad[,sample_cols]); storage.mode(cnt)<-"numeric"
ent <- AnnotationDbi::mapIds(org.Hs.eg.db,keys=aad$Symbol,keytype="SYMBOL",column="ENTREZID",multiVals="first")
cnt <- collapse_by_id(cnt,ent,"sum")
meta_aad <- data.frame(sample=colnames(cnt),group=factor(c(rep("Control",10),rep("Case",10)),levels=c("Control","Case")))
rownames(meta_aad)<-meta_aad$sample
y <- DGEList(cnt); keep <- filterByExpr(y,group=meta_aad$group); y<-y[keep,,keep.lib.sizes=FALSE]; y<-calcNormFactors(y)
design <- model.matrix(~meta_aad$group); v<-voom(y,design,plot=FALSE); fit<-eBayes(lmFit(v,design))
tt <- topTable(fit,coef=2,number=Inf,sort.by="none"); tt$ENTREZID<-rownames(tt); tt$FDR<-tt$adj.P.Val
sym <- AnnotationDbi::mapIds(org.Hs.eg.db,keys=tt$ENTREZID,keytype="ENTREZID",column="SYMBOL",multiVals="first")
tt$Symbol<-unname(sym[tt$ENTREZID]); tt<-tt[,c("ENTREZID","Symbol","logFC","AveExpr","t","P.Value","FDR","B")]
write.csv(tt,file.path(out,"GSE153434_all_genes.csv"),row.names=FALSE)
write.table(meta_aad,file.path(out,"GSE153434_sample_metadata.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
plot_qc(v$E,meta_aad,"GSE153434"); plot_de(tt,v$E,meta_aad,"GSE153434")
aad_gsea_path <- file.path(out,"GSE153434_six_pathway_GSEA.csv")
g_aad <- if(file.exists(aad_gsea_path)) read.csv(aad_gsea_path) else run_gsea(setNames(tt$t,tt$ENTREZID),"GSE153434")

# GSE38750: GEO normalized matrix; highest-variance probe retained per Entrez gene.
mic <- read_series_matrix("inputs/GSE38750_series_matrix.txt.gz")
probe_map <- AnnotationDbi::select(hgu133a2.db,keys=rownames(mic),keytype="PROBEID",columns=c("ENTREZID","SYMBOL"))
probe_map <- probe_map[!duplicated(probe_map$PROBEID),]
idx <- match(rownames(mic),probe_map$PROBEID); ent2 <- probe_map$ENTREZID[idx]
expr <- collapse_by_id(as.matrix(mic),ent2,"maxvar")
meta_ckd <- data.frame(sample=colnames(expr),group=factor(c(rep("Control",19),rep("Case",15)),levels=c("Control","Case")))
rownames(meta_ckd)<-meta_ckd$sample
design2<-model.matrix(~meta_ckd$group); fit2<-eBayes(lmFit(expr,design2)); tt2<-topTable(fit2,coef=2,number=Inf,sort.by="none")
tt2$ENTREZID<-rownames(tt2); tt2$FDR<-tt2$adj.P.Val
sym2<-AnnotationDbi::mapIds(org.Hs.eg.db,keys=tt2$ENTREZID,keytype="ENTREZID",column="SYMBOL",multiVals="first")
tt2$Symbol<-unname(sym2[tt2$ENTREZID]); tt2<-tt2[,c("ENTREZID","Symbol","logFC","AveExpr","t","P.Value","FDR","B")]
write.csv(tt2,file.path(out,"GSE38750_all_genes.csv"),row.names=FALSE)
write.table(meta_ckd,file.path(out,"GSE38750_sample_metadata.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(probe_map,file.path(out,"GPL571_probe_mapping.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
plot_qc(expr,meta_ckd,"GSE38750"); plot_de(tt2,expr,meta_ckd,"GSE38750")
ckd_gsea_path <- file.path(out,"GSE38750_six_pathway_GSEA.csv")
g_ckd <- if(file.exists(ckd_gsea_path)) read.csv(ckd_gsea_path) else run_gsea(setNames(tt2$t,tt2$ENTREZID),"GSE38750")

# GSE38751: independent participants with same-site internal mammary arteries;
# GEO sample annotations place five normal controls first and five uremic cases last.
ima <- read_series_matrix("inputs/GSE38751_series_matrix.txt.gz")
idx3 <- match(rownames(ima),probe_map$PROBEID); ent3 <- probe_map$ENTREZID[idx3]
expr3 <- collapse_by_id(as.matrix(ima),ent3,"maxvar")
meta_ima <- data.frame(sample=colnames(expr3),group=factor(c(rep("Control",5),rep("Case",5)),levels=c("Control","Case")))
rownames(meta_ima)<-meta_ima$sample
design3<-model.matrix(~meta_ima$group); fit3<-eBayes(lmFit(expr3,design3)); tt3<-topTable(fit3,coef=2,number=Inf,sort.by="none")
tt3$ENTREZID<-rownames(tt3); tt3$FDR<-tt3$adj.P.Val
sym3<-AnnotationDbi::mapIds(org.Hs.eg.db,keys=tt3$ENTREZID,keytype="ENTREZID",column="SYMBOL",multiVals="first")
tt3$Symbol<-unname(sym3[tt3$ENTREZID]); tt3<-tt3[,c("ENTREZID","Symbol","logFC","AveExpr","t","P.Value","FDR","B")]
write.csv(tt3,file.path(out,"GSE38751_all_genes.csv"),row.names=FALSE)
write.table(meta_ima,file.path(out,"GSE38751_sample_metadata.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
plot_qc(expr3,meta_ima,"GSE38751"); plot_de(tt3,expr3,meta_ima,"GSE38751")
ima_gsea_path <- file.path(out,"GSE38751_six_pathway_GSEA.csv")
g_ima <- if(file.exists(ima_gsea_path)) read.csv(ima_gsea_path) else run_gsea(setNames(tt3$t,tt3$ENTREZID),"GSE38751")

allg <- rbind(g_aad,g_ckd,g_ima); write.csv(allg,file.path(out,"human_validation_six_pathway_GSEA.csv"),row.names=FALSE)
allg$label <- factor(allg$pathway,levels=rev(go_names)); allg$sig <- allg$BH_q<0.05
p <- ggplot(allg,aes(NES,label,shape=sig,color=dataset)) + geom_vline(xintercept=0,color="grey70")+
  geom_point(size=3)+scale_color_manual(values=c(GSE153434="#D55E00",GSE38750="#0072B2",GSE38751="#009E73"))+
  scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1),labels=c(`TRUE`="Supported",`FALSE`="Not supported"))+
  labs(x="Normalized enrichment score",y=NULL,color="Dataset",shape="Six-pathway BH q<0.05",
  title="Human arterial evaluation of six frozen pathways")+
  guides(color=guide_legend(order=1),shape=guide_legend(order=2))+
  theme_classic(base_size=9)+theme(legend.position="bottom",legend.box="vertical",legend.box.just="left")
save_gg(p,"human_validation_GSEA_summary",7.2,4.8)

summ <- data.frame(dataset=c("GSE153434","GSE38750","GSE38751"),tested_genes=c(nrow(tt),nrow(tt2),nrow(tt3)),
  DEG_up=c(sum(tt$FDR<.05&tt$logFC>.5),sum(tt2$FDR<.05&tt2$logFC>.5),sum(tt3$FDR<.05&tt3$logFC>.5)),
  DEG_down=c(sum(tt$FDR<.05&tt$logFC< -.5),sum(tt2$FDR<.05&tt2$logFC< -.5),sum(tt3$FDR<.05&tt3$logFC< -.5)),
  pathways_q05=c(sum(g_aad$BH_q<.05),sum(g_ckd$BH_q<.05),sum(g_ima$BH_q<.05)))
write.table(summ,file.path(out,"human_validation_summary.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(out,"human_validation_sessionInfo.txt"))
cat("HUMAN VALIDATION COMPLETE\n"); print(summ); print(allg[,c("dataset","GO_ID","NES","P_value","BH_q")])
