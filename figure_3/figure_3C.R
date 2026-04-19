################################################################################
### Generate Figure 3c
### inferCNV hclust heatmap
################################################################################

# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"
#zenodo.dir <- "/Users/ycheng3/Projects/MPNST_phase1/"
# Define input directory with data and output directory to save the figure
input.dir <- paste0(zenodo.dir, "data/snRNA/inferCNV")
output.dir <- paste0(zenodo.dir, "results/figure_3/")
#dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
setwd(output.dir)

### =========================
### 0. Load libraries and data
### =========================

library(ComplexHeatmap)
library(tidyverse)
MPNST_obs_merge <- readRDS(paste0(input.dir,"/MPNST_obs_merge.rds"))
MPNST_hclust_cor_ward <- readRDS(paste0(input.dir,"/MPNST_hclust_cor_ward.rds"))
infercnv_chr_probes <- readRDS(paste0(input.dir,"/infercnv_chr_probes.rds"))
### =========================
### 1. cut into K clusters from hierarchical clustering 
### =========================                               

K = 24

### =========================
### 2. Plot with trees cut
### =========================     
color_scale <- colorRampPalette(c("darkblue", "white", "darkred"))(100)
png(filename = paste0(output.dir,"figure_3C.png"), width = 4000, height = 4000, res = 200)
pheatmap(mat = MPNST_obs_merge, cluster_rows = MPNST_hclust_cor_ward, clustering_distance_rows = "correlation", 
         show_rownames = F, show_colnames = F, cutree_rows = K, cluster_cols = F, color = color_scale, fontsize = 14, main = "MPNST Heatmap")
dev.off()

### =========================
### 3. Plot Original Heatmap
### =========================  
#print("Original Heatmap")
#ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(MPNST_obs_merge))),
#                      col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D")),
#                     annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 22), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
#png(filename = paste0("/Users/ycheng3/Projects/MPNST_phase1/figures/MPNST_C_merge_orig.png"), width = 4000, height = 4000, res = 200)
#infercnv_heatmap(MPNST_obs_merge, probes = infercnv_chr_probes$cum.probes %>% set_names(nm = c(1:22, "X")), row_ann = ha_row)
#Heatmap(MPNST_obs_merge, cluster_rows = F, cluster_columns=F, show_column_names = F, show_row_names = F, 
#      col = color_scale, row_names_gp=grid::gpar(fontsize=6), name="MPNST Heatmap", use_raster = T)
#dev.off()


