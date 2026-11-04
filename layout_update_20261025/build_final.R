# Display-only adaptation of the author-approved design. No analyses are run.
args<-commandArgs(trailingOnly=FALSE)
script<-sub('^--file=','',args[grepl('^--file=',args)][1]);setwd(dirname(normalizePath(script)))
dir.create('outputs',showWarnings=FALSE)
txt<-paste(readLines('approved_design.R'),collapse='\n')
txt<-strsplit(txt,"export(f2,'Figure2_candidate'",fixed=TRUE)[[1]][1]
txt<-sub("out <- 'visual_candidates_20261025'","out <- 'outputs'",txt,fixed=TRUE)
txt<-sub("src <- c(file.path('figures',files),'figures/Figure5_current_synthesis_source_data.tsv',paste0('figures/',c(stem,'Figure5_current_synthesis'),'.pdf'))","src <- c(file.path('figures',files),'figures/Figure5_current_synthesis_source_data.tsv')",txt,fixed=TRUE)
eval(parse(text=txt))
sc<-1;xo<-0;yo<-0
X<-function(x)unit((x+xo)*sc,'mm');Y<-function(y)unit((y+yo)*sc,'mm')
tx<-function(t,x,y,size=9,bold=FALSE,col=ink,just='left'){
 if(t %in% c('Shared pathways, context-dependent enrichment','An evidence map of shared arterial remodelling','Experimental overlap is broader than human support'))return(invisible(NULL))
 t<-sub('^a  ','A  ',t);t<-sub('^b  ','B  ',t);t<-sub('^c  ','C  ',t)
 if(grepl('GO smooth-muscle contraction is distinct',t))t<-'GO smooth-muscle contraction is distinct from the KEGG pathway.'
 grid.text(t,X(x),Y(y),just=just,gp=gpar(fontfamily='Arial',fontsize=max(7,size*sc),fontface=if(bold)'bold' else 'plain',col=col,lineheight=1.12))
}
rect<-function(x,y,w,h,fill,border=NA)grid.rect(X(x),Y(y),unit(w*sc,'mm'),unit(h*sc,'mm'),just=c('left','bottom'),gp=gpar(fill=fill,col=border))
line<-function(x1,y1,x2,y2,col='#D8E0E5')grid.lines(X(c(x1,x2)),Y(c(y1,y2)),gp=gpar(col=col,lwd=.6))
grid.circle<-function(x,y,r,gp)grid::grid.circle(X(convertX(x,'mm',valueOnly=TRUE)),Y(convertY(y,'mm',valueOnly=TRUE)),unit(convertWidth(r,'mm',valueOnly=TRUE)*sc,'mm'),gp=gp)
# Reuse the approved hero-panel instructions, without clipping a composite.
body2<-paste(deparse(body(f2)),collapse='\n')
start<-regexpr('tx("a  Full-ranked',body2,fixed=TRUE)[1]
end<-regexpr('tx("b  Animal ORA',body2,fixed=TRUE)[1]
stopifnot(start>0,end>start)
f2A<-eval(parse(text=paste0('function(){grid.newpage();',substr(body2,start,end-1),'}')))
f2B<-function(){grid.newpage();tx('B  Animal ORA (BH q)',8,55,10,TRUE)
 tx('Calcification',59,48,8,TRUE,just='centre');tx('AAD',90,48,8,TRUE,just='centre')
 short<-c('Matrix organization','Ossification','SM contraction','SM proliferation','Regulation of SM prolif.','Positive regulation')
 for(i in 1:6){y<-43-(i-1)*5.4;tx(short[i],8,y,8)
 for(j in 1:2){v<-a[a$pathway_id==ids[i]&a$model==cohorts[j],];tx(formatC(v$q,format='e',digits=2),c(59,90)[j],y,8,just='centre')}}
 tx('SM: smooth muscle; ORA: over-representation analysis',8,8,7)
}
f2C<-function(){grid.newpage();tx('C  Human sampling design',112,55,10,TRUE)
 tx('AD: ascending aorta, 10 vs 10\nUraemia: iliac vs renal artery, 15 vs 19\nUraemia: internal mammary artery, 5 vs 5',112,40,8.5)
 tx('GSE38750: disease and site confounded.\nGSE38751: same study, not independent replication.\nAnimal ORA and GSEA are distinct analyses.',112,21,8,col=grey)
}
save<-function(fun,name,w,h,scale=1,x=0,y=0){
 sc<<-scale;xo<<-x;yo<<-y
 for(ext in c('pdf','svg','png')){
 p<-file.path('outputs',paste0(name,'.',ext))
 if(ext=='pdf')quartz(type='pdf',file=p,width=w/25.4,height=h/25.4,bg='white')
 if(ext=='svg')svglite::svglite(p,width=w/25.4,height=h/25.4)
 if(ext=='png')ragg::agg_png(p,width=w,height=h,units='mm',res=300)
 fun();dev.off()
 }
}
save(f2,'Figure2_composite_review',190,178,190/210)
save(f2A,'Figure2_A',190,123,190/210,y=-60)
save(f2B,'Figure2_B',140,76,140/110)
save(f2C,'Figure2_C',140,76,140/110,x=-105)
save(f5,'Figure5',190,151,190/210)
for(n in c('Figure2_A','Figure2_B','Figure2_C','Figure5')){
 img<-magick::image_read(file.path('outputs',paste0(n,'.png')))
 magick::image_write(magick::image_convert(img,colorspace='gray'),file.path('outputs',paste0(n,'_gray.png')))
}
stopifnot(identical(before,setNames(vapply(src,hash,''),src)))
write_json(list(source_sha256=as.list(before),source_unchanged=TRUE,no_analysis_rerun=TRUE,
 panel_mapping=list(A='Previous B + C; 30 original NES values',B='Previous A; 12 original ORA q values',C='Previous D; 6 original sampling rows'),
 device='R Quartz PDF; SVG svglite; PNG ragg',dimensions_mm=list(Figure2_A=c(190,123),Figure2_B=c(140,76),Figure2_C=c(140,76),Figure5=c(190,151)),minimum_text_pt=7), 'outputs/BUILD.json',pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),'outputs/sessionInfo.txt')
