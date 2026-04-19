################################################################################
### Generate Figure 3a
### Plot cell type UMAP from snRNA data
################################################################################

# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"
#zenodo.dir <- "/Users/ycheng3/Projects/MPNST_phase1/"
# Define input directory with data and output directory to save the figure
input.dir <- paste0(zenodo.dir, "data/snRNA/seurat_objects")
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

seurat_obj<-readRDS(paste0(input.dir,"/MPNST_C_updated.rds"))
                    
### =========================
### 1. better plotting for Figure 3a
### =========================

umap_df <- Embeddings(seurat_obj, "umap") %>%
  as.data.frame() %>%
  mutate(
    cell_type = seurat_obj$cell_type
  )
centers <- umap_df %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
# separate endothelial
centers_endo <- centers %>% dplyr::filter(cell_type == "Endothelial")
centers_other <- centers %>% dplyr::filter(cell_type != "Endothelial")
png(file = paste0(output.dir, "/figure_3A.png"), width = 2500, height = 2500, res = 300)
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.3, alpha = 0.7) +
  
  # normal clusters (with lines)
  geom_text_repel(
    data = centers_other,
    aes(x = UMAP_1, y = UMAP_2, label = cell_type),
    inherit.aes = FALSE,
    color = "black",
    size = 4,
    box.padding = 1,
    point.padding = 0.3,
    force = 5,
    nudge_x = 2,
    nudge_y = 2,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.4,
    segment.alpha = 0.6,
    max.overlaps = Inf
  ) +
  
  # endothelial (no line, pushed right)
  geom_text_repel(
    data = centers_endo,
    aes(x = UMAP_1, y = UMAP_2, label = cell_type),
    inherit.aes = FALSE,
    color = "black",
    size = 4,
    box.padding = 1,
    point.padding = 0.3,
    force = 5,
    nudge_x = 4,   # <-- push more to the right
    nudge_y = -3,
    segment.color = NA,  # <-- remove line
    max.overlaps = Inf
  ) +
  
  theme_classic(base_size = 14) +
  theme(legend.position = "none")
dev.off()
