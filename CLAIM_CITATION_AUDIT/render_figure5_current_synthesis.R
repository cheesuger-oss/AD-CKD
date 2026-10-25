#!/usr/bin/env Rscript
# Text-only synchronization to the approved three-dataset human results.
# R grid is retained; equal panels do not encode quantitative evidence strength.
library(grid)
out <- "CLAIM_CITATION_AUDIT/aad_revision_candidate/figures"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
stages <- data.frame(
  stage=c("Animal discovery", "Human arterial tissue", "Genetic analyses", "Immune and classifier checks", "Study interpretation"),
  evidence=c(
    "549 shared differentially expressed genes and 155 shared pathways; matrix and smooth-muscle programmes in uraemic calcification and dissection-or-aneurysm models.",
    "Ossification passed correction in all three datasets. Matrix organization passed in two; smooth-muscle findings varied by dataset.",
    "No primary MR estimate passed BH correction. Genetic correlation was not estimable; colocalization was prior-sensitive and reverse MR was instrument-limited.",
    "No immune programme passed correction in either human immune-analysis dataset. Random-gene controls matched internal classification; external performance varied.",
    "Shared arterial-remodelling expression programmes across experimental models, with narrower, context-dependent human support strongest for ossification-related expression."
  ),
  boundary=c(
    "478 shared genes changed concordantly; 71 changed in opposite directions. Animal overlap is not dissection-specific.",
    "GSE38750 confounds disease with arterial site. GSE38751 is a small same-study, same-site comparison, not independent replication.",
    "Smaller effects and shared genetic architecture remain unresolved. The fitted regional model limits the colocalization inference.",
    "Programme scores are not cell fractions. Retrospective tissue classification does not test future disease risk.",
    "A focused basis for mechanistic studies; no causal mechanism, therapeutic target or treatment effect is established."
  ), stringsAsFactors=FALSE)
write.table(stages,file.path(out,"Figure5_current_synthesis_source_data.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
wrap <- function(x,w) paste(strwrap(x,width=w),collapse="\n")
draw <- function() {
 grid.newpage()
 grid.text("Shared arterial-remodelling programmes across evidence layers",x=.04,y=.96,just="left",gp=gpar(fontfamily="sans",fontsize=11,fontface="bold"))
 for(i in 1:5) {
  y <- .84-(i-1)*.17
  grid.rect(x=.5,y=y,width=.92,height=.155,gp=gpar(fill=if(i==5) "#EAF2F5" else "white",col="#D1D5DB"))
  grid.text(stages$stage[i],x=.055,y=y+.058,just="left",gp=gpar(fontfamily="sans",fontsize=8.2,fontface="bold",col="#235A73"))
  grid.text(wrap(stages$evidence[i],67),x=.055,y=y+.029,just=c("left","top"),gp=gpar(fontfamily="sans",fontsize=7.1,lineheight=1.15))
  grid.text(wrap(stages$boundary[i],43),x=.625,y=y+.029,just=c("left","top"),gp=gpar(fontfamily="sans",fontsize=6.7,lineheight=1.15,col="#444444"))
 }
 grid.text("Conceptual synthesis; equal box dimensions do not encode effect size or evidential weight.",x=.04,y=.035,just="left",gp=gpar(fontfamily="sans",fontsize=6.2,col="#555555"))
}
stopifnot(capabilities("aqua"))
grDevices::quartz(type="pdf",file=file.path(out,"Figure5_current_synthesis.pdf"),width=183/25.4,height=175/25.4)
draw(); dev.off()
png(file.path(out,"Figure5_current_synthesis.png"),width=183,height=175,units="mm",res=600,type="quartz",bg="white")
draw(); dev.off()
