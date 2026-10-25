options(stringsAsFactors=FALSE)
.libPaths(c("project/Rlib-4.5", "/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library", .libPaths()))
suppressPackageStartupMessages({
  library(edgeR); library(limma); library(AnnotationDbi); library(org.Hs.eg.db)
  library(hgu133a2.db); library(glmnet); library(ggplot2); library(pheatmap); library(jsonlite)
})

set.seed(20260809)
out <- "project/results/phase4/analysis"; fig <- file.path(out,"figures")
dir.create(fig,recursive=TRUE,showWarnings=FALSE)

read_series_matrix <- function(path) {
  x <- readLines(gzfile(path),warn=FALSE)
  a <- grep("^!series_matrix_table_begin",x)+1L; b <- grep("^!series_matrix_table_end",x)-1L
  if(length(a)!=1L || length(b)!=1L || b<a) stop("invalid series matrix: ",path)
  con <- textConnection(x[a:b]); on.exit(close(con))
  z <- read.delim(con,check.names=FALSE,quote="\"",row.names=1)
  as.matrix(z)
}

collapse_symbol <- function(mat,symbols) {
  symbols <- trimws(as.character(symbols)); keep <- !is.na(symbols)&nzchar(symbols)&!grepl(" /// |//|;",symbols)
  mat <- mat[keep,,drop=FALSE]; symbols <- symbols[keep]
  v <- apply(mat,1,var,na.rm=TRUE); ord <- order(symbols,-v,seq_along(symbols)); take <- !duplicated(symbols[ord])
  ans <- mat[ord[take],,drop=FALSE]; rownames(ans)<-symbols[ord[take]]; ans
}

read_platform_map <- function(path,symbol_col) {
  z<-read.delim(path,header=FALSE,quote="",comment.char="",fill=TRUE,check.names=FALSE)
  setNames(trimws(as.character(z[[symbol_col]])),trimws(as.character(z[[1]])))
}

rank_scores <- function(expr,membership,dataset,group) {
  ranks <- apply(expr,2,function(x) rank(x,ties.method="average")/length(x))
  if(is.null(dim(ranks))) ranks<-matrix(ranks,ncol=ncol(expr),dimnames=list(rownames(expr),colnames(expr)))
  programs <- unique(membership[,c("program_id","program_name","go_id")])
  scores<-matrix(NA_real_,nrow=nrow(programs),ncol=ncol(expr),dimnames=list(programs$program_id,colnames(expr)))
  observed<-integer(nrow(programs))
  for(i in seq_len(nrow(programs))){
    genes<-intersect(membership$human_gene[membership$program_id==programs$program_id[i]],rownames(expr))
    observed[i]<-length(genes); if(length(genes)<10) stop(dataset,": immune program coverage <10: ",programs$program_id[i])
    scores[i,]<-colMeans(ranks[genes,,drop=FALSE])
  }
  design<-model.matrix(~group); fit<-eBayes(lmFit(scores,design)); tt<-topTable(fit,coef=2,number=Inf,sort.by="none")
  tt<-tt[programs$program_id,,drop=FALSE]
  pooled_sd<-sapply(seq_len(nrow(scores)),function(i){sqrt(((sum(group=="Control")-1)*var(scores[i,group=="Control"])+(sum(group=="Case")-1)*var(scores[i,group=="Case"]))/(length(group)-2))})
  contrast<-data.frame(dataset=dataset,program_id=programs$program_id,program_name=programs$program_name,go_id=programs$go_id,
    observed_gene_n=observed,mean_control=rowMeans(scores[,group=="Control",drop=FALSE]),mean_case=rowMeans(scores[,group=="Case",drop=FALSE]),
    mean_difference=tt$logFC,standardized_effect=tt$logFC/pooled_sd,p_value=tt$P.Value,q_value=p.adjust(tt$P.Value,"BH"),stringsAsFactors=FALSE)
  long<-data.frame(dataset=dataset,program_id=rep(rownames(scores),ncol(scores)),sample=rep(colnames(scores),each=nrow(scores)),
    group=rep(as.character(group),each=nrow(scores)),score=as.vector(scores),stringsAsFactors=FALSE)
  list(scores=scores,contrast=contrast,long=long)
}

auc <- function(y,s) {
  n1<-sum(y==1); n0<-sum(y==0); if(n1==0||n0==0)return(NA_real_)
  (sum(rank(s,ties.method="average")[y==1])-n1*(n1+1)/2)/(n1*n0)
}
average_precision <- function(y,s){o<-order(s,decreasing=TRUE); y<-y[o]; sum(cumsum(y)/seq_along(y)*y)/sum(y)}

make_inner_foldid <- function(y,k=5,seed=1){set.seed(seed); f<-integer(length(y)); for(v in 0:1){idx<-sample(which(y==v)); f[idx]<-rep(seq_len(k),length.out=length(idx))}; f}

nested_oof <- function(X,y,seed=20260809,alpha=.5){
  pred<-rep(NA_real_,nrow(X)); selected<-integer(nrow(X)); lambdas<-numeric(nrow(X))
  for(i in seq_len(nrow(X))){
    tr<-setdiff(seq_len(nrow(X)),i); te<-i; xtr<-X[tr,,drop=FALSE]; xte<-X[te,,drop=FALSE]
    sdv<-apply(xtr,2,sd); keep<-is.finite(sdv)&sdv>0
    xtr<-xtr[,keep,drop=FALSE]; xte<-xte[,keep,drop=FALSE]
    if(ncol(xtr)>100){top<-names(sort(apply(xtr,2,var),decreasing=TRUE))[seq_len(100)]; xtr<-xtr[,top,drop=FALSE]; xte<-xte[,top,drop=FALSE]}
    mu<-colMeans(xtr); ss<-apply(xtr,2,sd); xtr<-sweep(sweep(xtr,2,mu,"-"),2,ss,"/"); xte<-sweep(sweep(xte,2,mu,"-"),2,ss,"/")
    fid<-make_inner_foldid(y[tr],5,seed+i)
    cv<-cv.glmnet(xtr,y[tr],family="binomial",alpha=alpha,foldid=fid,nlambda=20,standardize=FALSE,type.measure="deviance",maxit=100000)
    lambdas[i]<-cv$lambda.1se; co<-as.matrix(coef(cv,s="lambda.1se")); selected[i]<-sum(co[-1,1]!=0)
    pred[i]<-as.numeric(predict(cv,newx=xte,s="lambda.1se",type="response"))
  }
  if(any(!is.finite(pred)))stop("non-finite OOF prediction")
  list(pred=pred,selected=selected,lambdas=lambdas)
}

fit_final <- function(X,y,seed=20260809,alpha=.5){
  sdv<-apply(X,2,sd); X<-X[,is.finite(sdv)&sdv>0,drop=FALSE]
  if(ncol(X)>100){top<-names(sort(apply(X,2,var),decreasing=TRUE))[seq_len(100)]; X<-X[,top,drop=FALSE]}
  mu<-colMeans(X); ss<-apply(X,2,sd)
  Z<-sweep(sweep(X,2,mu,"-"),2,ss,"/"); fid<-make_inner_foldid(y,5,seed)
  cv<-cv.glmnet(Z,y,family="binomial",alpha=alpha,foldid=fid,nlambda=40,standardize=FALSE,type.measure="deviance",maxit=100000)
  co<-as.matrix(coef(cv,s="lambda.1se")); beta<-co[-1,1]; names(beta)<-rownames(co)[-1]; beta<-beta[beta!=0]
  if(length(beta)==0){co<-as.matrix(coef(cv,s="lambda.min")); beta<-co[-1,1]; names(beta)<-rownames(co)[-1]; beta<-beta[beta!=0]}
  if(length(beta)==0)stop("no nonzero final signature")
  list(beta=beta,intercept=co[1,1],lambda_1se=cv$lambda.1se,lambda_min=cv$lambda.min,mean=mu,sd=ss)
}

signature_rank_score <- function(expr,beta){genes<-intersect(names(beta),rownames(expr)); if(!length(genes))return(rep(NA_real_,ncol(expr))); r<-apply(expr,2,function(x)rank(x,ties.method="average")/length(x)); w<-beta[genes]; as.numeric(crossprod(w,r[genes,,drop=FALSE])/sum(abs(w)))}

bootstrap_auc_ci <- function(y,s,B=2000,seed=1){set.seed(seed); vals<-rep(NA_real_,B); for(b in seq_len(B)){idx<-c(sample(which(y==0),sum(y==0),TRUE),sample(which(y==1),sum(y==1),TRUE));vals[b]<-auc(y[idx],s[idx])}; quantile(vals,c(.025,.975),na.rm=TRUE)}

# Inputs and primary expression matrices
membership<-read.delim("project/results/phase4/immune_program_membership.tsv",check.names=FALSE)
candidates<-read.delim("project/results/phase4/renal_origin_candidates.tsv",check.names=FALSE)

rat<-read.delim(gzfile("project/data/geo/processed/GSE146638_expression.tsv.gz"),check.names=FALSE)
rat_meta<-read.delim("project/results/geo/GSE146638_sample_metadata.tsv",check.names=FALSE)
rat_meta<-rat_meta[rat_meta$include_primary=="yes",]; rat_mat<-as.matrix(rat[match(unique(rat$gene_symbol),rat$gene_symbol),rat_meta$expression_column]); storage.mode(rat_mat)<-"numeric"; rownames(rat_mat)<-unique(rat$gene_symbol)
rat_map<-read.csv("project/results/geo/ortholog_rat_human_1to1.csv"); mi<-match(rownames(rat_mat),rat_map$source_symbol); rat_h<-collapse_symbol(rat_mat,rat_map$human_symbol[mi])
y_rat<-DGEList(rat_h); y_rat<-calcNormFactors(y_rat); rat_expr<-cpm(y_rat,log=TRUE,prior.count=2); rat_group<-factor(ifelse(rat_meta$primary_group=="Uremic","Case","Control"),levels=c("Control","Case"))

aad<-read.delim(gzfile("inputs/GSE153434_all.counts.txt.gz"),check.names=FALSE)
aad_cols<-c(paste0("Normal",c(1,10,2:9)),paste0("Patient",c(10,1:9))); aad_cnt<-as.matrix(aad[,aad_cols]); storage.mode(aad_cnt)<-"numeric"
aad_expr0<-collapse_symbol(aad_cnt,aad$Symbol); aad_group<-factor(c(rep("Control",10),rep("Case",10)),levels=c("Control","Case"))
y_aad<-DGEList(aad_expr0); keep<-filterByExpr(y_aad,group=aad_group); y_aad<-calcNormFactors(y_aad[keep,,keep.lib.sizes=FALSE]); aad_expr<-voom(y_aad,model.matrix(~aad_group),plot=FALSE)$E

mic<-read_series_matrix("inputs/GSE38750_series_matrix.txt.gz")
pm<-AnnotationDbi::select(hgu133a2.db,keys=rownames(mic),keytype="PROBEID",columns=c("SYMBOL")); pm<-pm[!duplicated(pm$PROBEID),]; ck_expr<-collapse_symbol(mic,pm$SYMBOL[match(rownames(mic),pm$PROBEID)])
ck_group<-factor(c(rep("Control",19),rep("Case",15)),levels=c("Control","Case"))

# Immune analysis
imm1<-rank_scores(rat_expr,membership,"GSE146638",rat_group); imm2<-rank_scores(aad_expr,membership,"GSE153434",aad_group); imm3<-rank_scores(ck_expr,membership,"GSE38750",ck_group)
imm_con<-rbind(imm1$contrast,imm2$contrast,imm3$contrast); write.table(imm_con,file.path(out,"immune_program_contrasts.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
imm_long<-rbind(imm1$long,imm2$long,imm3$long); write.table(imm_long,file.path(out,"immune_program_scores.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
hm<-do.call(cbind,lapply(list(imm1$scores,imm2$scores,imm3$scores),function(z)t(scale(t(z)))))
ann<-data.frame(Dataset=c(rep("GSE146638",ncol(imm1$scores)),rep("GSE153434",ncol(imm2$scores)),rep("GSE38750",ncol(imm3$scores))),Group=c(as.character(rat_group),as.character(aad_group),as.character(ck_group))); rownames(ann)<-colnames(hm)
png(file.path(fig,"immune_program_heatmap.png"),width=2400,height=1300,res=300); pheatmap(hm,cluster_cols=FALSE,annotation_col=ann,fontsize=7,main="Frozen immune activation programs (within-dataset z scores)"); dev.off()
pimm<-ggplot(imm_con,aes(standardized_effect,program_name,color=dataset))+geom_vline(xintercept=0,color="grey70")+geom_point(position=position_dodge(width=.5),size=2.5)+scale_color_manual(values=c(GSE146638="#009E73",GSE153434="#D55E00",GSE38750="#0072B2"))+labs(x="Standardized case-control difference",y=NULL,title="Immune activation-program effects",subtitle="Scores are transcriptional programs, not cell fractions")+theme_classic(base_size=9)+theme(legend.position="bottom")
ggsave(file.path(fig,"immune_effect_forest.png"),pimm,width=7.2,height=4.8,dpi=300,bg="white")

# Primary renal-origin model
genes<-intersect(candidates$human_gene,rownames(aad_expr)); X<-t(aad_expr[genes,,drop=FALSE]); y<-as.integer(aad_group=="Case")
oof<-nested_oof(X,y); oof_auc<-auc(y,oof$pred); oof_ci<-bootstrap_auc_ci(y,oof$pred,2000,20260809); oof_ap<-average_precision(y,oof$pred)
oof_tab<-data.frame(sample=rownames(X),group=as.character(aad_group),truth=y,probability=oof$pred,outer_selected_n=oof$selected,outer_lambda=oof$lambdas)
write.table(oof_tab,file.path(out,"primary_oof_predictions.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

# Full-pipeline label permutations
perm_auc<-unlist(parallel::mclapply(seq_len(1000),function(b){set.seed(20260809+b); yp<-sample(y); auc(yp,nested_oof(X,yp,seed=20260809+b)$pred)},mc.cores=4,mc.set.seed=FALSE))
perm_p<-(1+sum(perm_auc>=oof_auc))/(1001); write.table(data.frame(iteration=seq_len(1000),auc=perm_auc),file.path(out,"permutation_auc.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

# Final signature and 500 bootstrap stability fits
final<-fit_final(X,y); coef_tab<-data.frame(gene=names(final$beta),coefficient=as.numeric(final$beta),renal_direction=candidates$renal_direction[match(names(final$beta),candidates$human_gene)])
write.table(coef_tab,file.path(out,"locked_signature_coefficients.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
stab_genes<-parallel::mclapply(seq_len(500),function(b){set.seed(20270000+b); idx<-c(sample(which(y==0),sum(y==0),TRUE),sample(which(y==1),sum(y==1),TRUE)); f<-try(fit_final(X[idx,,drop=FALSE],y[idx],seed=20270000+b),silent=TRUE); if(inherits(f,"try-error"))character() else names(f$beta)},mc.cores=4,mc.set.seed=FALSE)
stab<-setNames(integer(ncol(X)),colnames(X)); for(g in stab_genes)stab[g]<-stab[g]+1
stab_tab<-data.frame(gene=names(stab),selection_frequency=as.numeric(stab)/500); stab_tab<-stab_tab[order(-stab_tab$selection_frequency,stab_tab$gene),]; write.table(stab_tab,file.path(out,"feature_stability.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

# Size-matched random-gene negative controls using the same nested pipeline (100 repetitions)
noncand<-setdiff(colnames(t(aad_expr)),genes)
random_auc<-unlist(parallel::mclapply(seq_len(100),function(b){set.seed(20300000+b); g<-sample(noncand,length(genes),replace=FALSE); auc(y,nested_oof(t(aad_expr[g,,drop=FALSE]),y,seed=20300000+b)$pred)},mc.cores=4,mc.set.seed=FALSE))
write.table(data.frame(iteration=seq_len(100),auc=random_auc),file.path(out,"random_gene_control_auc.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

# External matrices in locked order
load_external<-function(dataset){
  if(dataset=="GSE52093"){m<-read_series_matrix("project/data/geo/phase4/raw/GSE52093_series_matrix.txt.gz"); mp<-read_platform_map("project/data/geo/phase4/platform/GPL10558-tbl-1.txt",13); m<-log2(pmax(m,1)); m<-normalizeBetweenArrays(m,method="quantile"); e<-collapse_symbol(m,mp[rownames(m)]); grp<-factor(c(rep("Case",4),rep("Control",3),rep("Case",3),rep("Control",2)),levels=c("Control","Case"))}
  else if(dataset=="GSE98770"){m<-read_series_matrix("project/data/geo/phase4/raw/GSE98770-GPL14550_series_matrix.txt.gz"); mp<-read_platform_map("project/data/geo/phase4/platform/GPL14550-tbl-1.txt",7); e<-collapse_symbol(m,mp[rownames(m)]); grp<-factor(c(rep("Case",6),rep("Control",5)),levels=c("Control","Case"))}
  else if(dataset=="GSE190635"){m<-read_series_matrix("project/data/geo/phase4/raw/GSE190635_series_matrix.txt.gz"); mp<-read_platform_map("project/data/geo/phase4/platform/GPL570-tbl-1.txt",11); e<-collapse_symbol(m,mp[rownames(m)]); grp<-factor(c(rep("Case",4),rep("Control",4)),levels=c("Control","Case"))}
  else stop("unknown external dataset")
  list(expr=e,group=grp)
}

external_order<-c("GSE52093","GSE98770","GSE147026","GSE190635"); ext_rows<-list(); ext_scores<-list()
for(ds in external_order){
  if(ds=="GSE147026"){ext_rows[[ds]]<-data.frame(dataset=ds,status="NOT_READY",case_n=4,control_n=4,signature_gene_n=NA,signature_coverage=NA,auc=NA,auc_ci_low=NA,auc_ci_high=NA,standardized_effect=NA,effect_se=NA,reason="metadata-only series matrix; zero expression rows"); next}
  z<-load_external(ds); score<-signature_rank_score(z$expr,final$beta); yy<-as.integer(z$group=="Case"); aa<-auc(yy,score); ci<-bootstrap_auc_ci(yy,score,2000,20260809+match(ds,external_order)); ng<-length(intersect(names(final$beta),rownames(z$expr))); cov<-ng/length(final$beta)
  mc<-mean(score[yy==0]); mt<-mean(score[yy==1]); sp<-sqrt(((sum(yy==0)-1)*var(score[yy==0])+(sum(yy==1)-1)*var(score[yy==1]))/(length(yy)-2)); d<-(mt-mc)/sp; se<-sqrt(length(yy)/(sum(yy==0)*sum(yy==1))+d^2/(2*(length(yy)-2)))
  status<-if(cov>=.5)"EVALUATED" else "INELIGIBLE_LOW_COVERAGE"
  ext_rows[[ds]]<-data.frame(dataset=ds,status=status,case_n=sum(yy==1),control_n=sum(yy==0),signature_gene_n=ng,signature_coverage=cov,auc=aa,auc_ci_low=ci[1],auc_ci_high=ci[2],standardized_effect=d,effect_se=se,reason="locked rank-weighted signature; no refit")
  ext_scores[[ds]]<-data.frame(dataset=ds,sample=colnames(z$expr),group=as.character(z$group),score=score)
}
ext<-do.call(rbind,ext_rows); rownames(ext)<-NULL; write.table(ext,file.path(out,"external_validation_results.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="")
write.table(do.call(rbind,ext_scores),file.path(out,"external_validation_scores.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

# Model figures
roc_df<-data.frame(threshold=sort(unique(c(-Inf,oof$pred,Inf)))); roc_df$tpr<-sapply(roc_df$threshold,function(t)mean(oof$pred[y==1]>=t)); roc_df$fpr<-sapply(roc_df$threshold,function(t)mean(oof$pred[y==0]>=t)); roc_df<-roc_df[order(roc_df$fpr,roc_df$tpr),]
proc<-ggplot(roc_df,aes(fpr,tpr))+geom_abline(slope=1,intercept=0,lty=2,color="grey60")+geom_step(color="#D55E00",linewidth=.8)+coord_equal()+labs(x="False-positive rate",y="True-positive rate",title=sprintf("Nested OOF ROC: AUC %.3f (95%% CI %.3f–%.3f)",oof_auc,oof_ci[1],oof_ci[2]),subtitle=sprintf("Full-pipeline permutation P=%.4f",perm_p))+theme_classic(base_size=10); ggsave(file.path(fig,"primary_oof_roc.png"),proc,width=5.5,height=5,dpi=300,bg="white")
ord<-order(oof$pred,decreasing=TRUE); prdf<-data.frame(recall=cumsum(y[ord])/sum(y),precision=cumsum(y[ord])/seq_along(y)); ppr<-ggplot(prdf,aes(recall,precision))+geom_step(color="#0072B2",linewidth=.8)+geom_hline(yintercept=mean(y),lty=2,color="grey60")+labs(x="Recall",y="Precision",title=sprintf("Nested OOF precision–recall: AP %.3f",oof_ap))+theme_classic(base_size=10); ggsave(file.path(fig,"primary_oof_pr.png"),ppr,width=5.5,height=5,dpi=300,bg="white")
cal<-data.frame(prob=oof$pred,truth=y,bin=cut(oof$pred,breaks=quantile(oof$pred,probs=seq(0,1,.2)),include.lowest=TRUE,duplicates="drop")); cal2<-aggregate(cbind(prob,truth)~bin,cal,mean); pcal<-ggplot(cal2,aes(prob,truth))+geom_abline(slope=1,intercept=0,lty=2,color="grey60")+geom_point(size=3,color="#009E73")+geom_line(color="#009E73")+coord_equal(xlim=c(0,1),ylim=c(0,1))+labs(x="Mean predicted probability",y="Observed case fraction",title="Exploratory calibration of nested OOF predictions")+theme_classic(base_size=10); ggsave(file.path(fig,"primary_calibration.png"),pcal,width=5.5,height=5,dpi=300,bg="white")
pst<-ggplot(head(stab_tab,30),aes(reorder(gene,selection_frequency),selection_frequency))+geom_col(fill="#0072B2")+coord_flip()+labs(x=NULL,y="Bootstrap selection frequency",title="Top renal-origin feature stability (500 resamples)")+theme_classic(base_size=9); ggsave(file.path(fig,"feature_stability.png"),pst,width=6.2,height=7,dpi=300,bg="white")
pcf<-ggplot(coef_tab,aes(reorder(gene,coefficient),coefficient,fill=coefficient>0))+geom_col()+coord_flip()+scale_fill_manual(values=c("TRUE"="#D55E00","FALSE"="#0072B2"),guide="none")+labs(x=NULL,y="Elastic-net coefficient",title="Locked renal-origin signature")+theme_classic(base_size=9); ggsave(file.path(fig,"locked_signature_coefficients.png"),pcf,width=6.2,height=max(4,0.22*nrow(coef_tab)+1.5),dpi=300,bg="white",limitsize=FALSE)
pext<-ggplot(ext[ext$status=="EVALUATED",],aes(auc,dataset))+geom_vline(xintercept=.5,lty=2,color="grey60")+geom_errorbarh(aes(xmin=auc_ci_low,xmax=auc_ci_high),height=.15)+geom_point(size=3,color="#CC79A7")+coord_cartesian(xlim=c(0,1))+labs(x="Locked-signature ROC AUC (bootstrap 95% CI)",y=NULL,title="External AAD validation in frozen order",subtitle="No cohort-specific refitting")+theme_classic(base_size=10); ggsave(file.path(fig,"external_validation_forest.png"),pext,width=6.5,height=4,dpi=300,bg="white")

summary<-list(primary_candidate_n=length(genes),oof_prediction_n=length(oof$pred),oof_auc=unname(oof_auc),oof_auc_ci=unname(oof_ci),average_precision=unname(oof_ap),permutations=1000,permutation_p=unname(perm_p),random_gene_controls=100,random_gene_auc_95th=unname(quantile(random_auc,.95)),renal_model_exceeds_random_95th=unname(oof_auc>quantile(random_auc,.95)),stability_resamples=500,locked_signature_n=length(final$beta),external_evaluated_n=sum(ext$status=="EVALUATED"),interpretation="prevalent tissue discrimination only; not incident-risk prediction")
write_json(summary,file.path(out,"model_summary.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,"phase4_sessionInfo.txt"))
cat("PHASE4 ANALYSIS COMPLETE\n"); print(summary); print(imm_con); print(ext)
