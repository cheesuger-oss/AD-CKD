library(grid)
library(jsonlite)
library(digest)
out <- 'visual_candidates_20261025'
read <- function(s) read.delim(file.path('figures',s),check.names=FALSE)
stem <- 'Figure2_pathway_convergence_human_support'
files <- paste0(stem,'_source_data_',c('A_ORA','B_animal_GSEA','C_human_GSEA','D_design'),'.tsv')
src <- c(file.path('figures',files),'figures/Figure5_current_synthesis_source_data.tsv',paste0('figures/',c(stem,'Figure5_current_synthesis'),'.pdf'))
hash <- function(p) digest(file=p,algo='sha256')
before <- setNames(vapply(src,hash,''),src)
a<-read(files[1]); b<-read(files[2]); c<-read(files[3]); d<-read(files[4]); s<-read('Figure5_current_synthesis_source_data.tsv')
stopifnot(nrow(a)==12,nrow(b)==12,nrow(c)==18,nrow(d)==6,nrow(s)==5)
ids<-c('GO:0030198','GO:0001503','GO:0006939','GO:0048659','GO:0048660','GO:0048661')
labs<-c('Extracellular matrix\norganization','Ossification','Smooth muscle\ncontraction','Smooth muscle cell\nproliferation','Regulation of smooth\nmuscle cell proliferation','Positive regulation of smooth\nmuscle cell proliferation')
dat<-rbind(data.frame(id=b$pathway_id,cohort=b$model,NES=b$NES,q=b$p.adjust,passed=b$passed),data.frame(id=c$pathway_id,cohort=c$dataset,NES=c$NES,q=c$BH_q,passed=c$passed))
cohorts<-c('Renal dysfunction','AAD lesions','GSE153434','GSE38750','GSE38751')
stopifnot(!anyDuplicated(paste(dat$id,dat$cohort)),all(dat$passed==(dat$q<.05)))
write.table(dat,file.path(out,'Figure2_display_values.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
ink<-'#243443'; blue<-'#0072B2'; orange<-'#D55E00'; grey<-'#63717A'
tx<-function(t,x,y,size=9,bold=FALSE,col=ink,just='left') grid.text(t,x=unit(x,'mm'),y=unit(y,'mm'),just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col,lineheight=1.12))
rect<-function(x,y,w,h,fill,border=NA) grid.rect(unit(x,'mm'),unit(y,'mm'),unit(w,'mm'),unit(h,'mm'),just=c('left','bottom'),gp=gpar(fill=fill,col=border))
line<-function(x1,y1,x2,y2,col='#D8E0E5') grid.lines(unit(c(x1,x2),'mm'),unit(c(y1,y2),'mm'),gp=gpar(col=col,lwd=.6))
wrap<-function(t,n) paste(strwrap(t,width=n),collapse='\n')
f2<-function(){
 grid.newpage();tx('Shared pathways, context-dependent enrichment',8,201,14,TRUE)
 tx('a  Full-ranked gene-set enrichment',8,190,10,TRUE)
 tx('Animal models',92,181,9,TRUE,just='centre');tx('Human arterial tissue',162,181,9,TRUE,just='centre')
 xs<-c(78,106,137,165,193);ys<-seq(153,83,length.out=6)
 heads<-c('Uraemic\ncalcification','AAD\nlesions','AD aorta\nGSE153434','Uraemia / site\nGSE38750','Uraemia / same site\nGSE38751')
 for(j in 1:5)tx(heads[j],xs[j],170,8,TRUE,just='centre')
 line(121,77,121,181)
 for(i in 1:6){rect(7,ys[i]-7,198,13,if(i%%2)'#F3F6F8' else '#FFFFFF');tx(labs[i],9,ys[i],8.5)
  for(j in 1:5){v<-dat[dat$id==ids[i]&dat$cohort==cohorts[j],];stopifnot(nrow(v)==1)
   grid.circle(unit(xs[j],'mm'),unit(ys[i]+2,'mm'),r=unit(2.5,'mm'),gp=gpar(fill=if(v$NES>0)orange else blue,col=NA))
   tx(sprintf('%+.2f%s',v$NES,if(v$passed)' *' else ''),xs[j],ys[i]-3.5,8.5,v$passed,just='centre')
  }}
 tx('Orange: positive NES   Blue: negative NES   * BH q < 0.05 within the original analysis',8,70,8)
 tx('Signed values preserve direction in grayscale. Animal and human correction families differ.',8,65,8,col=grey)
 tx('b  Animal ORA (BH q)',8,55,10,TRUE);tx('c  Human sampling design',112,55,10,TRUE)
 tx('Calcification',59,48,8,TRUE,just='centre');tx('AAD',90,48,8,TRUE,just='centre')
 short<-c('Matrix organization','Ossification','SM contraction','SM proliferation','Regulation of SM prolif.','Positive regulation')
 for(i in 1:6){y<-43-(i-1)*5.4;tx(short[i],8,y,8)
  for(j in 1:2){v<-a[a$pathway_id==ids[i]&a$model==cohorts[j],];tx(formatC(v$q,format='e',digits=2),c(59,90)[j],y,8,just='centre')}}
 tx('AD: ascending aorta, 10 vs 10\nUraemia: iliac vs renal artery, 15 vs 19\nUraemia: internal mammary artery, 5 vs 5',112,40,8.5)
 tx('GSE38750: disease and site confounded.\nGSE38751: same study, not independent replication.\nAnimal ORA and GSEA are distinct analyses.',112,21,8,col=grey)
 tx('GO smooth-muscle contraction is distinct from the KEGG pathway. No analysis or thresholds changed.',8,7,7.7,col=grey)
}
f5<-function(){
 grid.newpage();tx('An evidence map of shared arterial remodelling',8,181,14,TRUE)
 tx('Experimental overlap is broader than human support',8,172,10,col=grey)
 tx('EVIDENCE LAYER',8,160,8,TRUE);tx('WHAT THE RESULTS SHOW',49,160,8,TRUE);tx('SCOPE OF INTERPRETATION',133,160,8,TRUE)
 names<-c('Animal\ndiscovery','Human\narterial tissue','Genetic\nanalyses','Immune &\nclassifier checks')
 for(i in 1:4){top<-153-(i-1)*31
  rect(7,top-29,196,29,if(i%%2)'#F1F5F8' else '#FAFBFC');rect(7,top-29,1.2,29,if(i<=2)blue else '#8D9BA6')
  tx(sprintf('0%d',i),11,top-5,10,TRUE,col=if(i<=2)blue else grey);tx(names[i],11,top-16,10,TRUE)
  tx(wrap(s$evidence[i],48),49,top-14,8.2)
  tx(wrap(s$boundary[i],40),133,top-14,8.2,col=grey)
 }
 line(8,25,202,25);tx('INTERPRETATION',8,19,8,TRUE,col=blue)
 tx(wrap(s$evidence[5],110),49,19,8.5,TRUE)
 tx(s$boundary[5],8,5,8,col=grey)
}
export<-function(fun,name,w,h){
 for(ext in c('pdf','svg','png')){
  p<-file.path(out,paste0(name,'.',ext))
  if(ext=='pdf')quartz(type='pdf',file=p,width=w/25.4,height=h/25.4,bg='white')
  if(ext=='svg')svglite::svglite(p,width=w/25.4,height=h/25.4)
  if(ext=='png')ragg::agg_png(p,width=w,height=h,units='mm',res=200)
  fun();dev.off()
 }
}
export(f2,'Figure2_candidate',210,210);export(f5,'Figure5_candidate',210,190)
for(k in c(2,5)){
 old<-magick::image_read(file.path('figures',paste0(if(k==2)stem else 'Figure5_current_synthesis','.png')))
 new<-magick::image_read(file.path(out,paste0('Figure',k,'_candidate.png')))
 pair<-lapply(list(old,new),function(x)magick::image_extent(magick::image_resize(x,'1450x1450'), '1500x1500',gravity='center',color='white'))
 pair[[1]]<-magick::image_annotate(pair[[1]],'CURRENT',size=30,gravity='northwest',location='+25+10',color=ink)
 pair[[2]]<-magick::image_annotate(pair[[2]],'CANDIDATE - NOT REPLACED',size=30,gravity='northwest',location='+25+10',color=ink)
 magick::image_write(magick::image_append(do.call(base::c,pair)),file.path(out,paste0('Figure',k,'comparison.png')))
 magick::image_write(magick::image_convert(new,colorspace='gray'),file.path(out,paste0('Figure',k,'_grayscale.png')))
}
stopifnot(identical(before,setNames(vapply(src,hash,''),src)))
write_json(list(status='CANDIDATE_PENDING_VISUAL_REVIEW',source_sha256=as.list(before),rows=list(ORA=12,NES=30,design=6,synthesis=5),source_unchanged=TRUE,no_analysis_rerun=TRUE),file.path(out,'BUILD.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
