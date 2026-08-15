# Chromosomal instability shapes spatial and temporal phenotypic diversity in a malignant peripheral nerve sheath tumour
This repository contains code for the figures in the manuscript "Chromosomal instability shapes spatial and temporal phenotypic diversity in a malignant peripheral nerve sheath tumour".

The code for each main figure has its own folder:
- **figure_1**
- **figure_2**
- **figure_3**
- **figure_4**

## Installation

All scripts are written in R. Each figure folder contains a `sessionInfo.txt` listing the exact R version and package versions used to generate that figure. Install matching versions, e.g. with `renv` or `remotes::install_version()`:

```r
sessionInfo <- readLines("figure_1/sessionInfo.txt")  # inspect required package versions
```

Alternatively, install the latest CRAN/Bioconductor versions of the packages listed under "attached packages" / "other attached packages" and only pin specific versions if you hit compatibility issues.

## Running the code

1. Download the processed data from Zenodo: https://zenodo.org/records/19653314
2. In each script, set `zenodo.dir` to the local path where you downloaded the Zenodo data (it should contain `data/` and `results/` subfolders).
3. Run the scripts in a folder (e.g. `figure_1/figure_1B.R`) individually in R/RStudio; each script writes its output figure to `results/figure_N/` within the Zenodo directory. The total estimated run time is less than 30 minutes (Apple Silicon Mac, R 4.5.1, Zenodo data pre-downloaded locally).

Raw sequencing data (EGA) are not required to reproduce the figures — the Zenodo data already contains the processed inputs.
