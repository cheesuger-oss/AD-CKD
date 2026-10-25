# Correction-only graphical abstract: original R/grid composition and evidence boundaries.
root <- normalizePath('CLAIM_CITATION_AUDIT/identifier_correction_20261024')
txt <- readLines('CLAIM_CITATION_AUDIT/assemble_graphical_abstract_AAD.R',warn=FALSE)
txt <- gsub('155 shared pathways','154 shared pathways',txt,fixed=TRUE)
ix <- grep('^out_dir <-',txt)
stopifnot(length(ix)==1)
txt[ix] <- sprintf('out_dir <- %s',encodeString(file.path(root,'figures'),quote='"'))
eval(parse(text=txt),envir=new.env())
