################################################################################
### Generate Figure 3d
### Plot inferCNV clusters 8, 9, 10, 11
################################################################################

# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"
# Define input directory with data and output directory to save the figure
input.dir <- paste0(zenodo.dir, "data/snRNA")
output.dir <- paste0(zenodo.dir, "results/figure_3/")
#dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
setwd(output.dir)

### =========================
### 0. Load libraries and data
### =========================

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggrepel)
seurat_obj<-readRDS("/Users/ycheng3/Projects/MPNST_phase1/data/MPNST_C_updated.rds")

### =========================
### 1. keep only clusters 8, 9, 10, 11
### =========================
seurat_sub <- subset(
  seurat_obj,
  subset = inferCNV_clusters %in% c("8", "9", "10", "11")
)

seurat_sub$inferCNV_clusters <- factor(
  seurat_sub$inferCNV_clusters,
  levels = c("8", "9", "10", "11")
)
### =========================
### 2. plot UMAP
### =========================

cluster_colors_4 <- c(
  "8"  = "#BC80BD",
  "9"  = "#FFD700",
  "10" = "#00BFC4",
  "11" = "#FF1493"
)

p <- DimPlot(
  seurat_sub,
  reduction = "umap",
  group.by = "inferCNV_clusters",
  cols = cluster_colors_4,
  pt.size = 0.1
) +
  theme_classic()

png(file = paste0(output.dir, "/figure_3D.png"), width = 1500, height = 1350, res = 300)
print(p +
  coord_cartesian(
    xlim = c(-2, 6),
    ylim = c(6, 12.5)
  ))
dev.off()
