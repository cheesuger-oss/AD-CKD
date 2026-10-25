cran_mirror <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 2026-08-07 probe: common mainland mirrors did not expose Bioconductor 3.22.
# Keep TUNA for CRAN and use the verified Posit mirror only for Bioconductor.
bioc_mirror <- "https://bioconductor.posit.co/packages"

bioc_version <- if (getRversion() >= "4.5.0") {
  "3.22"
} else if (getRversion() >= "4.4.0") {
  "3.20"
} else {
  "3.18"
}

repositories <- c(
  BioCsoft = sprintf("%s/%s/bioc", bioc_mirror, bioc_version),
  BioCann = sprintf("%s/%s/data/annotation", bioc_mirror, bioc_version),
  BioCexp = sprintf("%s/%s/data/experiment", bioc_mirror, bioc_version),
  CRAN = cran_mirror
)

options(
  repos = repositories,
  BioC_mirror = bioc_mirror,
  timeout = 1800
)

required_cran <- c(
  "renv", "ggplot2", "pheatmap", "stringr", "dplyr", "data.table",
  "igraph", "MendelianRandomization"
)

required_bioc <- c(
  "GEOquery", "DESeq2", "limma", "clusterProfiler", "org.Hs.eg.db", "biomaRt"
)

required_github <- c(
  "MRCIEU/TwoSampleMR",
  "rondolab/MR-PRESSO"
)

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = cran_mirror)
}

if (!file.exists("renv.lock")) {
  renv::init(bare = TRUE)
}

renv::install(required_cran, repos = cran_mirror)

install.packages(
  required_bioc,
  lib = .libPaths()[1],
  repos = repositories,
  dependencies = NA
)
renv::install(required_github)
renv::snapshot(prompt = FALSE)

all_packages <- c(required_cran, required_bioc, "TwoSampleMR", "MRPRESSO")
stopifnot(all(vapply(all_packages, requireNamespace, logical(1), quietly = TRUE)))

message("Phase 0 R dependencies verified: ", paste(all_packages, collapse = ", "))
