project_library <- file.path(
  getwd(), "renv", "library", "macos", "R-4.5", "aarch64-apple-darwin20"
)
bootstrap_library <- file.path(getwd(), "project", "Rlib-4.5")

if (!dir.exists(project_library)) {
  stop("Phase 0 R dependency gate failed: project library is missing: ", project_library)
}

.libPaths(c(project_library, bootstrap_library, .libPaths()))

required <- c(
  "renv", "GEOquery", "DESeq2", "limma", "clusterProfiler",
  "org.Hs.eg.db", "biomaRt", "TwoSampleMR", "MendelianRandomization",
  "MRPRESSO", "ggplot2", "pheatmap", "stringr", "dplyr",
  "data.table", "igraph"
)

versions <- vapply(
  required,
  function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) return("MISSING")
    as.character(utils::packageVersion(pkg))
  },
  character(1)
)

report <- data.frame(package = required, version = unname(versions))
print(report, row.names = FALSE)

missing <- report$package[report$version == "MISSING"]
if (length(missing)) {
  stop("Phase 0 R dependency gate failed: ", paste(missing, collapse = ", "))
}

if (!file.exists("renv.lock")) {
  stop("Phase 0 R dependency gate failed: renv.lock is missing")
}

message("Phase 0 R dependency gate passed")
