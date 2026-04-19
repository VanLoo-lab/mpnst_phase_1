################################################################################
### Generate Figure 3.f
### meta program for snRNA data
################################################################################

# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"
# Define input directory with data and output directory to save the figure
input.dir <- paste0(zenodo.dir, "data/seurat_objects")
output.dir <- paste0(zenodo.dir, "results/figure_3/")
#dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
setwd(output.dir)
### =========================
### 0. Load libraries and data
### =========================

library(ggplot2)
library(Seurat)
library(SeuratObject)

hallmarks<-readRDS(paste0(input.dir,"/MPNST_hallmarks_cancer_states.rds"))

cancer_state = "Stress"
png(filename = paste0(output.dir,"figure3F_", cancer_state,".png"), width = 2000, height = 2000, res = 200)
FeaturePlot(hallmarks, features = paste0(cancer_state, "_RNA1"))
dev.off()

cancer_state = "Metal"
png(filename = paste0(output.dir,"figure3F_", cancer_state,".png"), width = 2000, height = 2000, res = 200)
FeaturePlot(hallmarks, features = paste0(cancer_state, "_RNA1"))
dev.off()
