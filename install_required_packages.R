repo_dir <- "/rsrch9/home/genetics/vanloolab/data/MPNST/phase1/260810_zenodo/MPNST-Zenodo"
lib_dir <- file.path(repo_dir, "library", "MPNST_phase1_zenodo")
summary_file <- file.path(repo_dir, "code", "mpnst_phase_1", "package_install_summary.tsv")

dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_dir, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))
dep_types <- c("Depends", "Imports", "LinkingTo")

cran_pkgs <- c(
  "ape", "dendextend", "ggplot2", "ggrepel", "gridExtra", "jpeg",
  "lubridate", "magrittr", "nnls", "patchwork", "pbmcapply",
  "phangorn", "pheatmap", "quadprog", "readxl", "Seurat",
  "SeuratObject", "stringr", "tidyverse", "umap", "wrMisc"
)

bioc_pkgs <- c(
  "Biobase", "BiocGenerics", "CARD", "ComplexHeatmap", "EpiDISH",
  "GenomeInfoDb", "GenomicRanges", "ggtree", "IRanges",
  "limma", "MatrixGenerics", "MuSiC", "S4Vectors",
  "SingleCellExperiment", "SummarizedExperiment", "tidytree",
  "TOAST", "treeio"
)

install_if_missing <- function(pkgs, installer) {
  installed_now <- rownames(installed.packages(lib.loc = .libPaths()))
  missing <- setdiff(pkgs, installed_now)
  if (length(missing) > 0) {
    installer(missing)
  }
}

lock_dirs <- Sys.glob(file.path(lib_dir, "00LOCK*"))
if (length(lock_dirs) > 0) {
  unlink(lock_dirs, recursive = TRUE, force = TRUE)
}

start_time <- Sys.time()

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib_dir, dependencies = dep_types)
}

suppressPackageStartupMessages(library(BiocManager, quietly = TRUE))

install_if_missing(
  bioc_pkgs,
  function(pkgs) BiocManager::install(pkgs, lib = lib_dir, ask = FALSE, update = FALSE, dependencies = dep_types)
)

install_if_missing(
  cran_pkgs,
  function(pkgs) install.packages(pkgs, lib = lib_dir, dependencies = dep_types)
)

end_time <- Sys.time()
installed_final <- rownames(installed.packages(lib.loc = .libPaths()))
missing_final <- setdiff(c(cran_pkgs, bioc_pkgs, "BiocManager"), installed_final)

summary_tbl <- data.frame(
  start_time = format(start_time, "%Y-%m-%d %H:%M:%S %Z"),
  end_time = format(end_time, "%Y-%m-%d %H:%M:%S %Z"),
  elapsed_seconds = round(as.numeric(difftime(end_time, start_time, units = "secs")), 2),
  library_dir = normalizePath(lib_dir, mustWork = FALSE),
  missing_after_install = if (length(missing_final) > 0) paste(missing_final, collapse = ",") else "",
  stringsAsFactors = FALSE
)

write.table(summary_tbl, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Library directory:", normalizePath(lib_dir, mustWork = FALSE), "\n")
cat("Elapsed seconds:", summary_tbl$elapsed_seconds, "\n")
if (length(missing_final) > 0) {
  cat("Missing after install:", paste(missing_final, collapse = ", "), "\n")
} else {
  cat("All requested packages are installed.\n")
}
