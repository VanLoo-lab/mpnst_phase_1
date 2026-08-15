# Example code for package installation

# Environment information:
# - Linux server environment (RHEL 8)
# - R 4.5.2
# - Bioconductor 3.22
# - Repo-local Zenodo data already pre-downloaded at checkout time
# - Local project library: ./library/MPNST_phase1_zenodo
# - Total installation time: ~ 1 hr

# Notes:
# - Uses a repo-local library under ./library/MPNST_phase1_zenodo
# - Installs CRAN packages and Bioconductor packages into that library
# - For environments where Bioconductor does not provide MuSiC/CARD,
#   install MuSiC from GitHub and CARD from GitHub after patching CARD's
#   src/Makevars to request C++14 instead of C++11.
# - On some systems, the runtime ICU libraries used by xml2/igraph are
#   visible only when LD_LIBRARY_PATH includes the Miniforge library dir.

repo_dir <- getwd()
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
  if (length(missing) > 0) installer(missing)
}

cleanup_lockdirs <- function(lib_dir) {
  lock_dirs <- Sys.glob(file.path(lib_dir, "00LOCK*"))
  if (length(lock_dirs) > 0) {
    unlink(lock_dirs, recursive = TRUE, force = TRUE)
  }
}

install_github_if_missing <- function(pkg, repo, lib_dir) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools", lib = lib_dir, dependencies = dep_types)
    }
    suppressPackageStartupMessages(library(devtools, quietly = TRUE))
    devtools::install_github(repo, lib = lib_dir, upgrade = "never", dependencies = TRUE)
  }
}

install_card_from_local_source <- function(card_src_dir, lib_dir) {
  stopifnot(dir.exists(card_src_dir))
  makevars <- file.path(card_src_dir, "src", "Makevars")
  makevars_win <- file.path(card_src_dir, "src", "Makevars.win")
  for (f in c(makevars, makevars_win)) {
    if (file.exists(f)) {
      x <- readLines(f, warn = FALSE)
      x <- gsub("CXX_STD = CXX11", "CXX_STD = CXX14", x, fixed = TRUE)
      writeLines(x, f)
    }
  }
  system2("R", c("CMD", "INSTALL", paste0("--library=", lib_dir), card_src_dir), stdout = TRUE, stderr = TRUE)
}

cleanup_lockdirs(lib_dir)
start_time <- Sys.time()

# Core installer: prefer Bioconductor/CRAN first.
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

# Optional fallback path for MuSiC and CARD when Bioconductor does not expose them.
install_github_if_missing("MuSiC", "xuranw/MuSiC", lib_dir)

if (!requireNamespace("CARD", quietly = TRUE)) {
  card_src_dir <- file.path(repo_dir, "code", "mpnst_phase_1", "vendor", "CARD-master")
  if (dir.exists(card_src_dir)) {
    install_card_from_local_source(card_src_dir, lib_dir)
  } else {
    install_github_if_missing("CARD", "YMa-lab/CARD", lib_dir)
  }
}

end_time <- Sys.time()
installed_final <- rownames(installed.packages(lib.loc = .libPaths()))
missing_final <- setdiff(c(cran_pkgs, bioc_pkgs, "BiocManager", "devtools"), installed_final)

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
