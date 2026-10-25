# Layout-only revision: frozen PCA, correlations, volcano and heatmap values.
# Contract: descriptive tissue-expression QC, not independent renal effects.
# Quantitative grid, R only, 183 x 240 mm, unchanged panels and source values.
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
root <- normalizePath(commandArgs(trailingOnly=TRUE)[1])
src <- file.path(root,'CLAIM_CITATION_AUDIT/recovered/AD-CKD/project/results/figures')
out <- file.path(root,'CLAIM_CITATION_AUDIT/aad_revision_candidate/figures')
theme_pub <- theme_classic(base_size=7.2)+theme(plot.title=element_text(face='bold',size=8.5),plot.subtitle=element_text(size=6.4),axis.title=element_text(face='bold'),plot.margin=margin(4,5,4,4))
embedded_pdf <- function(filename,width,height,bg,...) {
  stopifnot(capabilities('aqua'))
  grDevices::quartz(type='pdf',file=filename,width=width,height=height,bg=bg)
}
for (suffix in c('a','b')) {
  prefix <- paste0('Supplementary_Figure_S4',suffix)
  read_panel <- function(panel) read.delim(file.path(src,paste0(prefix,'_source_data_',panel,'.tsv')),check.names=FALSE)
  pd<-read_panel('pca');cd<-read_panel('cor');de<-read_panel('volcano');hd<-read_panel('heat')
  confounded<-suffix=='b'
  groups<-c('Control',if(confounded)'Uremic' else 'AAD')
  pd$group<-factor(pd$group,levels=groups)
  cd$display_x<-factor(cd$display_x,levels=pd$display);cd$display_y<-factor(cd$display_y,levels=rev(pd$display))
  hd$gene<-factor(hd$gene,levels=rev(unique(hd$gene)))
  hd$display<-factor(hd$display,levels=pd$display);hd$group<-factor(hd$group,levels=groups)
  de$status<-factor(de$status,levels=c('Down','Not DEG','Up'))
  de$neglogq<--log10(pmax(de$FDR,.Machine$double.xmin))
  stopifnot(nrow(pd)==if(confounded)34 else 20, length(levels(hd$gene))==40)
  p1<-ggplot(pd,aes(PC1,PC2,colour=group,shape=group))+geom_point(size=2)+scale_colour_manual(values=c(Control='#0072B2',AAD='#D55E00',Uremic='#CC79A7'))+labs(tag='A',title='PCA',subtitle=if(confounded)'1,000 most variable genes; uremic iliac versus control renal artery' else '1,000 most variable analysed genes; no PCA-based exclusion',x=sprintf('PC1 (%.1f%%)',pd$variance_PC1_pct[1]),y=sprintf('PC2 (%.1f%%)',pd$variance_PC2_pct[1]),colour=NULL,shape=NULL)+theme_pub+theme(legend.position='bottom')
  p2<-ggplot(cd,aes(display_x,display_y,fill=correlation))+geom_tile()+scale_fill_viridis_c(option='C',limits=c(min(cd$correlation),1),name='Pearson r')+labs(tag='B',title='Sample correlation',subtitle='Metadata order; display only',x=NULL,y=NULL)+theme_pub+theme(axis.text.x=element_text(angle=90,hjust=1,size=5.1),axis.text.y=element_text(size=5.1),axis.ticks=element_blank())
  p3<-ggplot(de,aes(logFC,neglogq,colour=status))+geom_point(size=.45,alpha=.55)+geom_vline(xintercept=c(-.5,.5),linetype='dashed',colour='grey55')+geom_hline(yintercept=-log10(.05),linetype='dashed',colour='grey55')+scale_colour_manual(values=c(Down='#0072B2','Not DEG'='grey78',Up='#D55E00'))+labs(tag='C',title='Differential-expression overview',subtitle='All tested genes; BH q<0.05 and abs(log2FC)>0.5',x='log2 fold change',y=expression(-log[10](BH~q)),colour=NULL)+theme_pub+theme(legend.position='bottom')
  p4<-ggplot(hd,aes(display,gene,fill=row_z))+geom_tile()+scale_fill_viridis_c(option='E',direction=-1,limits=c(-2.5,2.5),breaks=c(-2.5,0,2.5),oob=scales::squish,name='Row z-score')+facet_grid(~group,scales='free_x',space='free_x')+labs(tag='D',title='Selected-gene heatmap',subtitle='Top 40 by absolute moderated t; selected from this contrast\nLight = negative; dark = positive; colour limits at +/-2.5',x='Samples',y=NULL)+theme_pub+theme(axis.text.x=element_blank(),axis.text.y=element_text(size=5.5),axis.ticks.x=element_blank(),strip.background=element_blank(),strip.text=element_text(face='bold'))
  page<-((p1|p2)/(p3|p4))+plot_layout(heights=c(1,1.5))+plot_annotation(title=pd$dataset[1],subtitle=if(confounded)'Uremic iliac artery versus control renal artery; disease and arterial bed are completely confounded' else 'Ascending-aorta tissue: 10 controls and 10 aortic-dissection cases',theme=theme(plot.title=element_text(face='bold',size=10),plot.subtitle=element_text(size=7.2,face=if(confounded)'bold' else 'plain',colour=if(confounded)'#A33A2B' else 'black',margin=margin(b=4)))) & theme(plot.tag=element_text(face='bold',size=9),plot.tag.position=c(0,1))
  stem<-paste0(prefix,if(confounded)'_human_uremic_expression_QC' else '_human_AAD_expression_QC')
  ggsave(file.path(out,paste0(stem,'.pdf')),page,width=183,height=240,units='mm',device=embedded_pdf,bg='white')
  ggsave(file.path(out,paste0(stem,'.png')),page,width=183,height=240,units='mm',dpi=600,bg='white')
  for(panel in c('pca','cor','volcano','heat')) {
    name<-paste0(prefix,'_source_data_',panel,'.tsv')
    if(file.exists(file.path(out,name))) {
      stopifnot(unname(tools::md5sum(file.path(src,name)))==unname(tools::md5sum(file.path(out,name))))
    } else stopifnot(file.copy(file.path(src,name),file.path(out,name)))
  }
}
cat('S4 frozen layout export complete; no analysis rerun.\n')
