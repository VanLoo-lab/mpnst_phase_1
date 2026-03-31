### title: Plot highly expressed genes

# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"

# Define input directory with data and output directory to save the figure
input.dir <- paste0(zenodo.dir, "data/spRNA/")
output.dir <- paste0(zenodo.dir, "results/figure_4/")
setwd(output.dir)


####################################################################################################################################
### Part 0: Load libraries and data
####################################################################################################################################

### Load libraries
library(tidyverse)
library(Seurat)
library(patchwork)
library(pbmcapply)
library(readxl)

### Load data and update seurat object for the purposes here
MPNST_sp_markers <- readRDS(paste0(input.dir,"/MPNST_sp_markers.rds"))
MPNST_sp_markers <- lapply(MPNST_sp_markers, UpdateSeuratObject)


####################################################################################################################################
### Part 1: Make varfeature plot with top genes
####################################################################################################################################

### Check for R4 samples subclone 4
i = 6
de_markers <- FindMarkers(MPNST_sp_markers[[i]], ident.1 = 4, only.pos = F)

### Pick the features
pos_markers <- de_markers[de_markers$avg_log2FC > 0.25,]
pos_markers <- pos_markers[!grepl("RPL|RPS", rownames(pos_markers)),]

neg_markers <- de_markers[de_markers$avg_log2FC < -0.25,]
neg_markers <- neg_markers[!grepl("RPL|RPS", rownames(neg_markers)),]

SpatialFeaturePlot(
  MPNST_sp_markers[[i]], 
  features = rownames(pos_markers)[1:4], 
  ncol = 2,
  pt.size.factor = 60,
  alpha = c(0.1, 1)
  )

### Make the plots
png(filename = paste0("figure_4E.png"), width = 2000, height = 2000, res = 200)
print(
  SpatialFeaturePlot(
    MPNST_sp_markers[[i]], 
    features = rownames(pos_markers)[1:4], 
    ncol = 2, 
    pt.size.factor = 60,
    alpha = c(0.1, 1)))
dev.off()

# png(filename = paste0("MPNST_varfeatureplot_neg_R4_4.png"), width = 2000, height = 2000, res = 200)
# print(
#   SpatialFeaturePlot(
#     MPNST_sp_markers[[i]], 
#     features = rownames(neg_markers)[1:4], 
#     ncol = 2, 
#     pt.size.factor = 60,
#     alpha = c(0.1, 1)))
# dev.off()
