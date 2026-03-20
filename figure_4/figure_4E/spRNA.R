
####################################################################################################################################
### Part 0: Libraries and preperation 
####################################################################################################################################

#options(bitmapType='cairo') #to solve plotting issue on CAMP

library(tidyverse)
library(Seurat) #previously had 3.1.4, now 3.2.1
library(patchwork)
library(pbmcapply)
library(readxl)

# samples = c("P", "R1", "R3", "P-b", "R2", "R4", "R5", "R5-b") #changes names for my own runs (previous labelling was Anneliens)
samples = c("Pa", "R1", "R3", "Pb", "R2", "R4", "R5a", "R5b")
names(samples) = c("VER683A1", "VER683A2", "VER683A3", "VER683A4", "VER683A5", "VER683A6", "VER683A7", "VER683A8")
# regions = c("Pa", "R1", "R3", "Pb", "R2", "R4", "R5a", "R5b") #names for barcodes can't have hyphens i think
slide = c(rep("V19N11-046",4), rep("V19N11-047",4))
area = rep(c("D1", "C1", "B1", "A1"), 2)

gg_colours <- hcl(h = seq(15, 375, length = length(samples) + 1), l = 65, c = 100)[1:length(samples)]
color_scale <- colorRampPalette(c("green", "darkred"))(100)
# barplot(1:length(samples), col = gg_colours)

# images = c("V19N11-046_D1_P_VER683A1.jpg", "V19N11-046_C1_.R1_VER683A2.jpg", "V19N11-046_B1_R3_VER683A3.jpg", "V19N11-046_A1_P-b_VER683A4.jpg", 
#            "V19N11-047_D1_R2_VER683A5.jpg", "V19N11-047_C1_R4_VER683A6.jpg", "V19N11-047_B1_R5_VER683A7.jpg", "V19N11-047_A1_R5-b_VER683A8.jpg")

scrna.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_RNA/results/C/")
input.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_spatial/data/")
system(paste0("mkdir -p ", input.dir))

rawrun1.dir = "/camp/project/proj-vanloo/analyses/averfaillie/10x_spatial/PM_20021/run1/"
rawrun2.dir = "/camp/project/proj-vanloo/analyses/averfaillie/10x_spatial/PM_20021/run2/"
rawrun3.dir = "/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_spatial/data/PM_20021/run3/201203_K00102_0542_AHK5YNBBXY/"

transcriptome.ref = "/camp/project/proj-vanloo/analyses/hyan/ref_files/SpaceRanger/refdata-gex-GRCh38-2020-A/"
# transcriptome.ref = "/srv/shared/vanloo/pipeline-files/human/references/SpaceRanger/refdata-cellranger-GRCh38-3.0.0"

output.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_spatial/results/separate/")
system(paste0("mkdir -p ", output.dir))
setwd(output.dir)

# system(paste0("cp -r ", rawrun1.dir, "images ", input.dir))
# system(paste0("cp -r ", rawrun1.dir, "slide_layout ", input.dir))



####################################################################################################################################
### Part 1: Create MPNST_sp object
####################################################################################################################################

if(F) {
  Load10X_Spatial_hy <- function(data.dir, assay = 'Spatial', slice = 'slice1', filter.matrix = TRUE, to.upper = FALSE, ...) {
    data <- Read10X(data.dir)
    if (to.upper) {
      rownames(x = data) <- toupper(x = rownames(x = data))
    }
    object <- CreateSeuratObject(counts = data, assay = assay)
    image <- Read10X_Image(image.dir = file.path(data.dir, '../spatial'), filter.matrix = filter.matrix)
    image <- image[Cells(x = object)]
    DefaultAssay(object = image) <- assay
    object[[slice]] <- image
    return(object)
  }
  
  MPNST_list <- list()
  for (i in 1:length(samples)) {
    MPNST_list[[i]] <- Load10X_Spatial_hy(data.dir = paste0(input.dir, "run1_2_3_man_align/",samples[i],"/outs/filtered_feature_bc_matrix"))
  }
  
  saveRDS(MPNST_list, "MPNST_list.rds")
  MPNST_list <- readRDS("MPNST_list.rds")
  
  MPNST_SCT <- pbmclapply(MPNST_list, function(x) {
    SCTransform(x, assay = "Spatial", verbose = FALSE)
  }, mc.cores = 8)
  
  saveRDS(MPNST_SCT, "MPNST_SCT.rds")
  MPNST_SCT <- readRDS("MPNST_SCT.rds")
  
  #Check certain gene expression
  # SpatialFeaturePlot(MPNST_SCT[[1]], features = c("TTN", "CD163", "PTPRC", "VCAN"), pt.size.factor = 1.2, alpha = c(0.1,1))
  
  #Perform dimension retuction
  Seurat_pipe <- function(s) {
    s <- RunPCA(s, assay = "SCT", verbose = FALSE)
    s <- FindNeighbors(s, reduction = "pca", dims = 1:30)
    s <- FindClusters(s, verbose = FALSE)
    s <- RunUMAP(s, reduction = "pca", dims = 1:30)
    return(s)
  }
  
  MPNST_sp <- pbmclapply(MPNST_SCT, Seurat_pipe, mc.cores = 8)
  saveRDS(MPNST_sp, "MPNST_sp.rds")
  MPNST_sp <- readRDS("MPNST_sp.rds")
  
  
  #Check variable features by spatial patterning
  MPNST_sp_markers <- pbmclapply(MPNST_sp, function(x) {
    FindSpatiallyVariableFeatures(x, assay = "SCT", features = VariableFeatures(x)[1:1000], selection.method = "markvariogram")
  }, mc.cores = 8)
  
  saveRDS(MPNST_sp_markers, "MPNST_sp_markers.rds")
  MPNST_sp_markers <- readRDS("MPNST_sp_markers.rds")
}

####################################################################################################################################
### Part 2: Make varfeature plot with top genes (Figure 4E)
####################################################################################################################################

setwd("~/Documents/GitHub/MPNST-Zenodo/figure_4/results/figure_4E")

### Read in the spatial transcriptomic data (seurat object)
MPNST_sp_markers <- readRDS("~/Documents/GitHub/MPNST-Zenodo/figure_1/data/10X_spatial/MPNST_sp_markers.rds")
MPNST_sp_markers <- lapply(MPNST_sp_markers, UpdateSeuratObject)

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
png(filename = paste0("MPNST_varfeatureplot_pos_R4_4.png"), width = 2000, height = 2000, res = 200)
print(
  SpatialFeaturePlot(
    MPNST_sp_markers[[i]], 
    features = rownames(pos_markers)[1:4], 
    ncol = 2, 
    pt.size.factor = 60,
    alpha = c(0.1, 1)))
dev.off()

png(filename = paste0("MPNST_varfeatureplot_neg_R4_4.png"), width = 2000, height = 2000, res = 200)
print(
  SpatialFeaturePlot(
    MPNST_sp_markers[[i]], 
    features = rownames(neg_markers)[1:4], 
    ncol = 2, 
    pt.size.factor = 60,
    alpha = c(0.1, 1)))
dev.off()
