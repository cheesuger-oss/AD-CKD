# Layout-only version at the journal's stated 130 mm readable width.
library(grid)
root <- normalizePath('CLAIM_CITATION_AUDIT/identifier_correction_20261024')
out <- file.path(root,'journal_figure_candidate','figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
draw <- function() {
 grid.newpage()
 navy <- '#17365D'
 grid.text('Cross-species evidence suggests a bounded\narterial-remodelling overlap',.5,.94,
  gp=gpar(fontfamily='Helvetica',fontsize=9,fontface='bold',col=navy,lineheight=1.15))
 heads <- c('Animal discovery','Preliminary human support','Boundary tests')
 bodies <- c(
 '549 shared genes\n478 concordant; 71 discordant\n154 shared pathways\nUraemic calcification and dissection-or-aneurysm models',
 'Ossification: 3 of 3 datasets\nMatrix and smooth-muscle signals varied\n5, 10 and 15 cases across datasets',
 'No renal trait passed BH correction\nGenetic correlation was non-estimable\nClassifier signal lacked specificity')
 ys <- c(.755,.535,.315)
 for(i in 1:3) {
  grid.roundrect(.5,ys[i],.92,.185,r=unit(1.5,'mm'),
   gp=gpar(fill=c('#DCE8F5','#E3F0E8','#F3EEE5')[i],col='#64748B',lwd=.6))
  grid.text(heads[i],.065,ys[i]+.064,just='left',gp=gpar(fontfamily='Helvetica',fontsize=8,fontface='bold',col=navy))
  grid.text(bodies[i],.065,ys[i]+.032,just=c('left','top'),
   gp=gpar(fontfamily='Helvetica',fontsize=7,lineheight=1.2,col='#222222'))
  if(i<3)grid.lines(x=c(.5,.5),y=c(ys[i]-.096,ys[i]-.121),
   arrow=arrow(type='closed',length=unit(1.5,'mm')),gp=gpar(col='#2F6DA5',fill='#2F6DA5',lwd=.7))
 }
 grid.roundrect(.5,.105,.92,.17,r=unit(1.5,'mm'),gp=gpar(fill=navy,col=navy))
 grid.text('Most reproducible human signal:\nossification-related expression',.5,.145,
  gp=gpar(fontfamily='Helvetica',fontsize=8,fontface='bold',col='white',lineheight=1.1))
 grid.text('The overlap guides focused studies but does not establish\ncausality or nominate treatment targets',.5,.067,
  gp=gpar(fontfamily='Helvetica',fontsize=7,col='white',lineheight=1.15))
}
pdf <- file.path(out,'GraphicalAbstract_R1.pdf')
# Longest side remains horizontal: 130 x 120 mm.
grDevices::quartz(type='pdf',file=pdf,width=130/25.4,height=120/25.4,family='Helvetica')
draw();dev.off()
stopifnot(system2('/opt/homebrew/bin/pdftocairo',c('-svg',shQuote(pdf),shQuote(file.path(out,'GraphicalAbstract_R1.svg'))))==0)
stopifnot(system2('/opt/homebrew/bin/pdftoppm',c('-png','-singlefile','-r','600',shQuote(pdf),shQuote(file.path(out,'GraphicalAbstract_R1'))))==0)
