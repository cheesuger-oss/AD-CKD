args <- commandArgs(trailingOnly = TRUE)
.libPaths(c("project/rlib", .libPaths()))
suppressPackageStartupMessages(library(coloc))
eligible_window <- function(n) is.finite(n) && n >= 100
sanitize <- function(x) gsub("[[:space:]]+", " ", trimws(x))
drop_duplicate_rsids <- function(d) d[!(duplicated(d$rsid) | duplicated(d$rsid, fromLast=TRUE)), , drop=FALSE]

if (length(args) == 1 && args[1] == "--selftest") {
  stopifnot(!eligible_window(99), eligible_window(100))
  stopifnot(sanitize("a\nb\tc") == "a b c")
  dx <- data.frame(rsid=c("rs1","rs1","rs2"), value=1:3)
  stopifnot(identical(drop_duplicate_rsids(dx)$rsid, "rs2"))
  set.seed(20260809)
  n <- 500; maf <- runif(n, 0.05, 0.5); v <- rep(0.01^2, n)
  b1 <- rnorm(n, 0, sqrt(v)); b2 <- rnorm(n, 0, sqrt(v)); b1[100] <- 0.15; b2[100] <- 0.16
  d1 <- list(beta=b1, varbeta=v, MAF=maf, N=100000, type="quant", snp=paste0("rs",seq_len(n)))
  d2 <- list(beta=b2, varbeta=v, MAF=maf, N=100000, type="quant", snp=paste0("rs",seq_len(n)))
  z <- coloc.abf(d1, d2, p1=1e-4, p2=1e-4, p12=1e-5)
  stopifnot(z$summary[["PP.H4.abf"]] > 0.8)
  cat("SELFTEST PASS\n"); quit(status=0)
}

if (length(args) != 2) stop("usage: run_coloc_all_windows.R input_dir output.tsv")
indir <- args[1]; outfile <- args[2]
traits <- c("egfrcrea", "egfrcys", "cystatin_c", "bun", "ckd", "uacr")
p12s <- c(1e-6, 5e-6, 1e-5, 5e-5)
out <- list(); k <- 0

for (trait in traits) {
  path <- file.path(indir, paste0("aad__", trait, "_harmonized.tsv.gz"))
  d <- read.delim(gzfile(path), check.names=FALSE)
  for (wid in unique(d$window_id)) {
    x <- d[d$window_id == wid, , drop=FALSE]
    n_before_dedup <- nrow(x); x <- drop_duplicate_rsids(x); n_dup_excluded <- n_before_dedup - nrow(x)
    if (!eligible_window(nrow(x))) next
    d1 <- list(beta=x$renal_beta, varbeta=x$renal_varbeta, MAF=pmin(x$renal_freq, 1-x$renal_freq),
               N=median(x$renal_n), type=if (trait == "ckd") "cc" else "quant", snp=x$rsid)
    if (trait == "ckd") d1$s <- 41395 / 480698
    d2 <- list(beta=x$aad_beta, varbeta=x$aad_varbeta, MAF=pmin(x$aad_freq, 1-x$aad_freq),
               N=519870, s=1262/519870, type="cc", snp=x$rsid)
    for (prior in p12s) {
      k <- k + 1; warning_text <- ""; error_text <- ""
      z <- tryCatch(withCallingHandlers({capture.output(result <- coloc.abf(d1, d2, p1=1e-4, p2=1e-4, p12=prior)); result},
                  warning=function(w) { warning_text <<- paste(warning_text, conditionMessage(w), sep=" | "); invokeRestart("muffleWarning") }),
                    error=function(e) { error_text <<- conditionMessage(e); NULL })
      if (is.null(z)) {
        out[[k]] <- data.frame(trait=trait, window_id=wid, nsnps=nrow(x), duplicate_rows_excluded=n_dup_excluded, p12=prior,
          PP.H0=NA, PP.H1=NA, PP.H2=NA, PP.H3=NA, PP.H4=NA, H4_over_H3H4=NA,
          status="ERROR", warning=sanitize(warning_text), error=sanitize(error_text))
      } else {
        s <- z$summary; h3 <- unname(s[["PP.H3.abf"]]); h4 <- unname(s[["PP.H4.abf"]])
        out[[k]] <- data.frame(trait=trait, window_id=wid, nsnps=unname(s[["nsnps"]]), duplicate_rows_excluded=n_dup_excluded, p12=prior,
          PP.H0=unname(s[["PP.H0.abf"]]), PP.H1=unname(s[["PP.H1.abf"]]),
          PP.H2=unname(s[["PP.H2.abf"]]), PP.H3=h3, PP.H4=h4,
          H4_over_H3H4=ifelse(h3+h4>0, h4/(h3+h4), NA), status="OK", warning=sanitize(warning_text), error="")
      }
    }
  }
  cat(trait, length(unique(d$window_id)), "windows\n")
}
write.table(do.call(rbind, out), outfile, sep="\t", row.names=FALSE, quote=FALSE)
