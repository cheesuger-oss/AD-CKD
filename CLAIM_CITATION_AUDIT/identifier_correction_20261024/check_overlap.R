#!/usr/bin/env Rscript
i <- normalizePath('CLAIM_CITATION_AUDIT/identifier_correction_20261024')
rat <- read.csv(file.path(i,'results/corrected_DEG.csv'))
mouse <- read.csv(file.path(i,'frozen/results/geo/GSE325608_all_genes.csv'))
rm <- read.csv(file.path(i,'ortholog_rat_human_1to1_corrected.csv'))
mm <- read.csv(file.path(i,'frozen/data/geo/processed/mouse_human_biomart_all.csv'))
mm <- subset(mm,orthology_type=='ortholog_one2one' & orthology_confidence==1)
collapse <- function(d,m) {
  d <- merge(subset(d,adj.P.Val<.05 & abs(logFC)>.5),m,by.x='gene_symbol',by.y='source_symbol')
  aggregate(logFC~human_symbol,d,function(x) x[which.max(abs(x))])
}
x <- merge(collapse(rat,rm),collapse(mouse,mm),by='human_symbol',suffixes=c('_renal','_ad'))
old <- read.delim(file.path(dirname(i),'recovered/AD-CKD/project/results/figures/Figure1_animal_transcriptomic_discovery_source_data_D_overlap.tsv'))
old <- old[match(x$human_symbol,old$human_symbol),]
stopifnot(nrow(x)==549,identical(x$human_symbol,old$human_symbol))
d <- max(abs(x$logFC_renal-old$logFC_renal),abs(x$logFC_ad-old$logFC_ad))
stopifnot(d<1e-9,sum(sign(x$logFC_renal)==sign(x$logFC_ad))==478)
write.csv(data.frame(status='PASS',overlap=549,concordant=478,discordant=71,max_logFC_delta=d),file.path(i,'overlap_verification.csv'),row.names=FALSE)
print(read.csv(file.path(i,'overlap_verification.csv')))
