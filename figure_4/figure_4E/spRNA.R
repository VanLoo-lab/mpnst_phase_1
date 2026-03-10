options(bitmapType='cairo') #to solve plotting issue on CAMP
library(tidyverse, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(Seurat, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/") #previously had 3.1.4, now 3.2.1
library(patchwork, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(pbmcapply, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(readxl, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")

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
###21/09/20 Run spaceranger on MPNST
####################################################################################################################################
#Running spaceranger on top upped samples (run on instance)
#nohup Rscript spRNA.R &> spRNA_p.out &
if (F) {
  for (i in 2:length(samples)) {
  # i = 1
    system(paste0("mkdir -p ", input.dir, "run1"))
    setwd(paste0(input.dir, "run1"))
    system(paste0("spaceranger count --id=", samples[i],
           " --transcriptome=", transcriptome.ref,
           " --fastqs=", rawrun1.dir, "fastq",
           " --sample=", names(samples)[i],
           " --image=", input.dir, "images/", slide[i],"_",area[i],"_",(samples)[i],"_",names(samples)[i],".jpg",
           " --slide=", slide[i],
           " --area=", area[i],
           " --slidefile=", input.dir, "slide_layout/", slide[i], ".gpr",
           " --localcores=6"))
  }
  #/srv/sw/eb/software/SpaceRanger/1.0.0/bin/spatial_rna/count --id=R1_slidelayout --transcriptome=/srv/shared/vanloo/pipeline-files/human/references/SpaceRanger/refdata-cellranger-GRCh38-3.0.0 --fastqs=fastq --sample=VER683A2 --image=images/V19N11-046_C1_.R1_VER683A2.jpg --slide=V19N11-046 --area=C1 --slidefile=slide_layout/V19N11-046.gpr --localcores=6
}

#Use fastqs from both sequencing runs
if (F) {
  for (i in 1:length(samples)) {
    # i = 1
    system(paste0("mkdir -p ", input.dir, "run1_2")) #Reran with 10X supplied transcriptome
    setwd(paste0(input.dir, "run1_2_new"))
    system(paste0("spaceranger count --id=", samples[i],
                  " --transcriptome=", transcriptome.ref,
                  " --fastqs=", rawrun1.dir, "fastq", ",",rawrun2.dir, "fastq",
                  " --sample=", names(samples)[i],
                  " --image=", input.dir, "images/", slide[i],"_",area[i],"_",(samples)[i],"_",names(samples)[i],".jpg",
                  " --slide=", slide[i],
                  " --area=", area[i],
                  " --slidefile=", input.dir, "slide_layout/", slide[i], ".gpr",
                  " --localcores=8"))
  }
}

#Use fastqs from both sequencing runs AND manually aligned
if (F) {
  for (i in 1:length(samples)) {
    # i = 1
    system(paste0("mkdir -p ", input.dir, "run1_2_man_align")) #Reran with 10X supplied transcriptome
    setwd(paste0(input.dir, "run1_2_man_align"))
    system(paste0("spaceranger count --id=", samples[i],
                  " --transcriptome=", transcriptome.ref,
                  " --fastqs=", rawrun1.dir, "fastq", ",",rawrun2.dir, "fastq",
                  " --sample=", names(samples)[i],
                  " --image=", input.dir, "images/", slide[i], "_", area[i], "_", samples[i], "_", names(samples)[i], ".jpg",
                  " --slide=", slide[i],
                  " --area=", area[i],
                  " --loupe-alignment=", input.dir, "images/", slide[i], "-", area[i], ".json",
                  " --slidefile=", input.dir, "slide_layout/", slide[i], ".gpr",
                  " --localcores=8"))
  }
}

#Run with spaceranger version 1.1.0
if (F) {
  for (i in 1:length(samples)) {
    # i = 1
    system(paste0("mkdir -p ", input.dir, "run1_2_man_align")) #Reran with 10X supplied transcriptome
    setwd(paste0(input.dir, "run1_2_man_align"))
    system(paste0("/srv/shared/vanloo/home/hyan/sw/Spaceranger/spaceranger-1.1.0/bin/spaceranger count --id=", samples[i],
                  " --transcriptome=", transcriptome.ref,
                  " --fastqs=", rawrun1.dir, "fastq", ",",rawrun2.dir, "fastq",
                  " --sample=", names(samples)[i],
                  " --image=", input.dir, "images/", slide[i], "_", area[i], "_", samples[i], "_", names(samples)[i], ".jpg",
                  " --slide=", slide[i],
                  " --area=", area[i],
                  " --loupe-alignment=", input.dir, "images/", slide[i], "-", area[i], ".json",
                  " --slidefile=", input.dir, "slide_layout/", slide[i], ".gpr",
                  " --localcores=8"))
  }
}

#Run with spaceranger version 1.1.0 and all 3 sequencing runs
if (F) {
  for (i in 1:length(samples)) {
    # i = 1
    system(paste0("mkdir -p ", input.dir, "run1_2_3_man_align")) #Reran with 10X supplied transcriptome
    setwd(paste0(input.dir, "run1_2_3_man_align"))
    system(paste0("/srv/shared/vanloo/home/hyan/sw/Spaceranger/spaceranger-1.1.0/bin/spaceranger count --id=", samples[i],
                  " --transcriptome=", transcriptome.ref,
                  " --fastqs=", rawrun1.dir, "fastq", ",", rawrun2.dir, "fastq", ",", rawrun3.dir, "fastq",
                  " --sample=", names(samples)[i],
                  " --image=", input.dir, "images/", slide[i], "_", area[i], "_", samples[i], "_", names(samples)[i], ".jpg",
                  " --slide=", slide[i],
                  " --area=", area[i],
                  " --loupe-alignment=", input.dir, "images/", slide[i], "-", area[i], ".json",
                  " --slidefile=", input.dir, "slide_layout/", slide[i], ".gpr",
                  " --localcores=8"))
  }
}

####################################################################################################################################
###24/09/20 Run on regions separately on MPNST
####################################################################################################################################
if (F) {
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
  
  for (i in 1:length(samples)) {
    plot1 <- VlnPlot(MPNST_list[[i]], features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
    plot2 <- SpatialFeaturePlot(MPNST_list[[i]], features = "nCount_Spatial", pt.size.factor = 1.2) + theme(legend.position = "right")
    png(filename = paste0("MPNST_nCount_spatial_",samples[i],".png"), width = 4000, height = 2000, res = 200)
    print(wrap_plots(plot1, plot2))
    dev.off()
  }
  
  for (i in 1:length(samples)) {
    plot1 <- VlnPlot(MPNST_list[[i]], features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
    plot2 <- SpatialFeaturePlot(MPNST_list[[i]], features = "nFeature_Spatial", pt.size.factor = 1.2) + theme(legend.position = "right")
    png(filename = paste0("MPNST_nFeature_spatial_",samples[i],".png"), width = 4000, height = 2000, res = 200)
    print(wrap_plots(plot1, plot2))
    dev.off()
  }
  
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
  
  for (i in 1:length(samples)) {
    # png(filename = paste0("MPNST_dimplot_",samples[i],".png"), width = 4000, height = 2000, res = 200)
    # clusters <- length(levels(MPNST_sp[[i]]))
    # plot_colours <- hcl(h = seq(15, 375, length = length(levels(MPNST_sp[[i]])) + 1), l = 65, c = 100)[1:length(levels(MPNST_sp[[i]]))]
    # p1 <- DimPlot(MPNST_sp[[i]], reduction = "umap", label = TRUE)
    # p2 <- SpatialDimPlot(MPNST_sp[[i]], cols = NULL, label = TRUE, label.size = 5, pt.size.factor = 1.2) + scale_fill_manual(values = plot_colours)
    # print(p1 + p2)
    # dev.off()
    
    pdf(file = paste0("MPNST_dimplot_",samples[i],".pdf"), width = 14, height = 7)
    clusters <- length(levels(MPNST_sp[[i]]))
    plot_colours <- hcl(h = seq(15, 375, length = length(levels(MPNST_sp[[i]])) + 1), l = 65, c = 100)[1:length(levels(MPNST_sp[[i]]))]
    p1 <- DimPlot(MPNST_sp[[i]], reduction = "umap", label = TRUE)
    p2 <- SpatialDimPlot(MPNST_sp[[i]], cols = NULL, label = TRUE, label.size = 5, pt.size.factor = 1.2, stroke = NA) + scale_fill_manual(values = plot_colours) +
      labs(fill = "Cluster") + guides(fill = guide_legend(override.aes = list(size=3)))
    print(p1 + p2)
    dev.off()
    
    pdf(file = paste0("MPNST_spatialdimplot_",samples[i],".pdf"), width = 7, height = 7)
    clusters <- length(levels(MPNST_sp[[i]]))
    plot_colours <- hcl(h = seq(15, 375, length = length(levels(MPNST_sp[[i]])) + 1), l = 65, c = 100)[1:length(levels(MPNST_sp[[i]]))]
    p2 <- SpatialDimPlot(MPNST_sp[[i]], cols = NULL, label = TRUE, label.size = 8, pt.size.factor = 1, stroke = NA) + scale_fill_manual(values = plot_colours) +
      labs(fill = "Cluster") + theme(legend.position = "none", legend.text=element_text(size=20)) + guides(colour = guide_legend(override.aes = list(size=8)))
    print(p2)
    dev.off()
  }
  
  #Check marker genes for scRNA data
  for (i in 1:length(samples)) {
    DefaultAssay(MPNST_sp[[i]]) <- "Spatial"
    png(filename = paste0("MPNST_map_by_gene_",samples[i],".png"), width = 4000, height = 4000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp[[i]], features = c("F13A1", "CD163", #macrophage
                                                         "PTPRB", "VWF", #endothelial
                                                         "TTN", "NEB", #Skeletal muscle
                                                         "PTPRC", "IKZF1", "THEMIS"), ncol = 3, alpha = c(0.1,1)) +  #T cell
            scale_color_gradientn(colours = c('green', 'darkred')))
    dev.off()
    DefaultAssay(MPNST_sp[[i]]) <- "SCT"
  }
  
  #Check marker genes for scRNA data
  for (i in 1:length(samples)) {
    DefaultAssay(MPNST_sp[[i]]) <- "Spatial"
    png(filename = paste0("MPNST_map_by_gene_",samples[i],".png"), width = 4000, height = 4000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp[[i]], features = c("F13A1", "CD163", #macrophage
                                                         "PTPRB", "VWF", #endothelial
                                                         "TTN", "NEB", #Skeletal muscle
                                                         "PTPRC", "IKZF1", "THEMIS"), ncol = 3, alpha = c(0.1,1)) +  #T cell
            scale_color_gradientn(colours = c('green', 'darkred')))
    dev.off()
    DefaultAssay(MPNST_sp[[i]]) <- "SCT"
  }
  
  #Check FOS expression for scRNA data
  for (i in 1:length(samples)) {
    DefaultAssay(MPNST_sp[[i]]) <- "Spatial"
    png(filename = paste0("MPNST_map_by_gene_FOS_",samples[i],".png"), width = 2000, height = 2000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp[[i]], features = c("FOS"), alpha = c(0.1,1)) +
            scale_color_gradientn(colours = c('green', 'darkred')))
    dev.off()
    DefaultAssay(MPNST_sp[[i]]) <- "SCT"
  }
  
  #Check homo del genes (Chr 2 and 18)
  for (i in 1:length(samples)) {
    DefaultAssay(MPNST_sp[[i]]) <- "Spatial"
    png(filename = paste0("MPNST_map_by_homo_del_",samples[i],".png"), width = 4000, height = 4000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp[[i]], features = c("SH3YL1", "ACP1",
                                                         "MAPK4", "MRO",
                                                         "ENO1"), ncol = 3, alpha = c(0.1,1)))
    dev.off()
    DefaultAssay(MPNST_sp[[i]]) <- "SCT"
  }
  
  # #Check variable features by clusters (not run yet)
  # de_markers <- FindMarkers(brain, ident.1 = 4, ident.2 = 6)
  # SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
  
  #Check variable features by spatial patterning
  MPNST_sp_markers <- pbmclapply(MPNST_sp, function(x) {
    FindSpatiallyVariableFeatures(x, assay = "SCT", features = VariableFeatures(x)[1:1000], selection.method = "markvariogram")
  }, mc.cores = 8)
  
  saveRDS(MPNST_sp_markers, "MPNST_sp_markers.rds")
  MPNST_sp_markers <- readRDS("MPNST_sp_markers.rds")
  
  for (i in 1:length(samples)) {
    MPNST_top_markers <- head(SpatiallyVariableFeatures(MPNST_sp_markers[[i]], selection.method = "markvariogram"), 9)
    png(filename = paste0("MPNST_varfeatureplot_",samples[i],".png"), width = 4000, height = 4000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp_markers[[i]], features = MPNST_top_markers, ncol = 3, alpha = c(0.1, 1)))
    dev.off()
  }
  
  #Check for R4 samples subclone 4
  i = 6
  de_markers <- FindMarkers(MPNST_sp_markers[[i]], ident.1 = 4, only.pos = F)
  pos_markers <- de_markers[de_markers$avg_logFC > 0.25,]
  pos_markers <- pos_markers[!grepl("RPL|RPS", rownames(pos_markers)),]
  neg_markers <- de_markers[de_markers$avg_logFC < -0.25,]
  neg_markers <- neg_markers[!grepl("RPL|RPS", rownames(neg_markers)),]
  png(filename = paste0("MPNST_varfeatureplot_pos_R4_4.png"), width = 2000, height = 2000, res = 200)
  print(SpatialFeaturePlot(MPNST_sp_markers[[i]], features = rownames(pos_markers)[1:4], ncol = 2, alpha = c(0.1, 1)))
  dev.off()
  
  png(filename = paste0("MPNST_varfeatureplot_neg_R4_4.png"), width = 2000, height = 2000, res = 200)
  print(SpatialFeaturePlot(MPNST_sp_markers[[i]], features = rownames(neg_markers)[1:4], ncol = 2, alpha = c(0.1, 1)))
  dev.off()
}

####################################################################################################################################
###Integrate with scRNA data (NOT RERUN)
####################################################################################################################################
if (F) {
  MPNST_non_CC <- readRDS(paste0(scrna.dir, "MPNST_C_non_CC.rds"))
  new.cluster.ids <- c("Malignant_R5", "Malignant_R2", "Malignant_R4_2", "Malignant_R3", "Malignant_R1_1",
                       "Malignant_P", "Malignant_R4_1", "Monocyte", "Malignant_R1_2", "Malignant_Ribo",
                       "Malignant_G2MS_R2R3", "Malignant_G2MS_R1R5", "Malignant_G2MS_R4", "T_cell", "Malignant_R4_3",
                       "Endothelial", "Skeletal_Muscle")
  names(new.cluster.ids) <- levels(MPNST_non_CC)
  MPNST_non_CC <- RenameIdents(MPNST_non_CC, new.cluster.ids)
  MPNST_non_CC$cell_type <- MPNST_non_CC@active.ident
  
  #Load in cluster markers
  MPNST_non_CC_cluster_markers <- readRDS(paste0(scrna.dir, "MPNST_C_non_CC_cluster_markers.rds"))
  
  ##Make table with more info on top 10 genes
  n = 100
  system(paste0("mkdir -p ", output.dir, n, "_genes_module"))
  setwd(paste0(output.dir, n, "_genes_module"))
  top_n <- as_tibble(MPNST_non_CC_cluster_markers) %>% group_by(cluster) %>% arrange(cluster, -avg_logFC) %>% top_n(n = n, wt = avg_logFC)

  RNA_modules <- lapply(1:length(levels(MPNST_non_CC)), function(c) {
    top_n %>% filter(cluster == c-1) %>% pull(gene)
  })
  
  for (c in 1:length(levels(MPNST_non_CC))) {
    print(paste0("Adding module score for ", levels(MPNST_non_CC)[c], " cluster"))
    MPNST_non_CC <- AddModuleScore(MPNST_non_CC, features = RNA_modules[c], ctrl = 100, name = levels(MPNST_non_CC)[c])
  }
  # saveRDS(MPNST_non_CC, "MPNST_C_non_CC_module_scores.rds")
  #Note module names shouldn't have spaces so renamed to T.cell and Skeletal.Muscle
  
  for (c in new.cluster.ids) {
    DefaultAssay(MPNST_non_CC) <- "RNA"
    png(filename = paste0("MPNST_C_non_CC_module_",c,".png"), width = 2000, height = 2000, res = 200)
    print(FeaturePlot(MPNST_non_CC, features = paste0(c,1)))
    dev.off()
  }

  #Add module score to spatial
  for (i in 1:length(samples)) {
    print(samples[i])
    for (c in 1:length(levels(MPNST_non_CC))) {
      print(paste0("Adding module score for ", levels(MPNST_non_CC)[c], " cluster"))
      MPNST_sp_markers[[i]] <- AddModuleScore(MPNST_sp_markers[[i]], features = RNA_modules[c], ctrl = 100, name = levels(MPNST_non_CC)[c])
    }
  }
  
  for (i in 1:length(samples)) {
    print(samples[i])
    DefaultAssay(MPNST_sp_markers[[i]]) <- "Spatial"
    
    png(filename = paste0("MPNST_dimplot_",samples[i],"_module_1_to_9.png"), width = 4000, height = 2000, res = 200)
    print(FeaturePlot(MPNST_sp_markers[[i]], features = paste0(new.cluster.ids[1:9],1), label = TRUE, ncol = 3))
    dev.off()

    png(filename = paste0("MPNST_dimplot_",samples[i],"_module_10_to_17.png"), width = 4000, height = 2000, res = 200)
    print(FeaturePlot(MPNST_sp_markers[[i]], features = paste0(new.cluster.ids[10:17],1), label = TRUE, ncol = 3))
    dev.off()
    
    png(filename = paste0("MPNST_spatial_",samples[i],"_module_1_to_9.png"), width = 4000, height = 2000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp_markers[[i]], features = paste0(new.cluster.ids[1:9],1), ncol = 3, alpha = c(0.1,1)) + theme(legend.position = "right"))
    dev.off()
    
    png(filename = paste0("MPNST_spatial_",samples[i],"_module_10_to_17.png"), width = 4000, height = 2000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp_markers[[i]], features = paste0(new.cluster.ids[10:17],1), ncol = 3, alpha = c(0.1,1)) + theme(legend.position = "right"))
    dev.off()
    DefaultAssay(MPNST_sp_markers[[i]]) <- "SCT"
  }
  setwd(paste0(output.dir))
}  

#Try modules between regions and tumour cells.
if (F) {
  MPNST_non_CC <- readRDS(paste0(scrna.dir, "MPNST_C_non_CC.rds"))
  new.cluster.ids <- c("Malignant_R5", "Malignant_R2", "Malignant_R4_2", "Malignant_R3", "Malignant_R1_1",
                       "Malignant_P", "Malignant_R4_1", "Monocyte", "Malignant_R1_2", "Malignant_Ribo",
                       "Malignant_G2MS_R2R3", "Malignant_G2MS_R1R5", "Malignant_G2MS_R4", "T_cell", "Malignant_R4_3",
                       "Endothelial", "Skeletal_Muscle")
  names(new.cluster.ids) <- levels(MPNST_non_CC)
  MPNST_non_CC <- RenameIdents(MPNST_non_CC, new.cluster.ids)
  MPNST_non_CC$cell_type <- MPNST_non_CC@active.ident
  
  library(plyr)
  coarse.cluster.ids <- MPNST_non_CC@active.ident
  coarse.cluster.ids <- revalue(coarse.cluster.ids, c("Malignant_R1_1"="Malignant_R1", "Malignant_R1_2"="Malignant_R1", 
                                                      "Malignant_R4_1"="Malignant_R4", "Malignant_R4_2"="Malignant_R4", "Malignant_R4_3"="Malignant_R4",
                                                      "Malignant_G2MS_R2R3"="Malignant_G2MS", "Malignant_G2MS_R1R5"="Malignant_G2MS", "Malignant_G2MS_R4"="Malignant_G2MS"))
  MPNST_non_CC$coarse <- coarse.cluster.ids
  Idents(MPNST_non_CC) <- MPNST_non_CC$coarse
  png(filename = paste0("MPNST_C_non_CC_umap_coarse_raw.png"), width = 2000, height = 1000, res = 200)
  DimPlot(MPNST_non_CC, reduction = "umap")
  dev.off()
  
  #Find cluster markers
  MPNST_non_CC_cluster_markers <- FindAllMarkers(MPNST_non_CC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  saveRDS(MPNST_non_CC_cluster_markers, file = paste0("MPNST_C_non_CC_coarse_cluster_markers.rds"))
  MPNST_non_CC_cluster_markers <- readRDS(paste0("MPNST_C_non_CC_corase_cluster_markers.rds"))
  
  ##Make table with more info on top 10 genes
  annotations <- readRDS("/srv/shared/vanloo/home/hyan/mpnst/10X_RNA/results/annotations.rds")
  top10 <- as.tbl(MPNST_non_CC_cluster_markers) %>% group_by(cluster) %>% arrange(cluster, -avg_logFC) %>% top_n(n = 10, wt = avg_logFC)
  MPNST_non_CC_cluster_top10_markers <- left_join(top10, annotations[, c(1:2, 3, 9)], by = c("gene" = "gene_name"))
  MPNST_non_CC_cluster_top10_markers <- MPNST_non_CC_cluster_top10_markers[, c(6:8, 1:5, 9:10)]
  write.csv(MPNST_non_CC_cluster_top10_markers, paste0("MPNST_C_non_CC_coarse_top10.csv"), quote = F)
  
  ##Make table with more info on top 10 genes
  n = 100
  system(paste0("mkdir -p ", output.dir, n, "_coarse_genes_module"))
  setwd(paste0(output.dir, n, "_coarse_genes_module"))
  top_n <- as_tibble(MPNST_non_CC_cluster_markers) %>% group_by(cluster) %>% arrange(cluster, -avg_logFC) %>% top_n(n = n, wt = avg_logFC)
  
  RNA_modules <- lapply(1:length(levels(MPNST_non_CC)), function(c) {
    top_n %>% filter(cluster == levels(MPNST_non_CC)[c]) %>% pull(gene)
  })
  
  for (c in 1:length(levels(MPNST_non_CC))) {
    print(paste0("Adding module score for ", levels(MPNST_non_CC)[c], " cluster"))
    MPNST_non_CC <- AddModuleScore(MPNST_non_CC, features = RNA_modules[c], ctrl = 100, name = levels(MPNST_non_CC)[c])
  }
  # saveRDS(MPNST_non_CC, "MPNST_C_non_CC_module_scores.rds")
  #Note module names shouldn't have spaces so renamed to T.cell and Skeletal.Muscle
  
  for (c in levels(coarse.cluster.ids)) {
    DefaultAssay(MPNST_non_CC) <- "RNA"
    png(filename = paste0("MPNST_C_non_CC_coarse_module_",c,".png"), width = 2000, height = 2000, res = 200)
    print(FeaturePlot(MPNST_non_CC, features = paste0(c,1)))
    dev.off()
  }
  
  png(filename = paste0("MPNST_C_non_CC_coarse_module_1_to_6.png"), width = 2000, height = 2000, res = 200)
  print(FeaturePlot(MPNST_non_CC, features = paste0(levels(coarse.cluster.ids)[1:6],1)), ncol = 2)
  dev.off()
  
  png(filename = paste0("MPNST_C_non_CC_coarse_module_7_to_12.png"), width = 2000, height = 2000, res = 200)
  print(FeaturePlot(MPNST_non_CC, features = paste0(levels(coarse.cluster.ids)[7:12],1)), ncol = 2)
  dev.off()
  
  #Add module score to spatial
  for (i in 1:length(samples)) {
    print(samples[i])
    for (c in 1:length(levels(MPNST_non_CC))) {
      print(paste0("Adding module score for ", levels(MPNST_non_CC)[c], " cluster"))
      MPNST_sp_markers[[i]] <- AddModuleScore(MPNST_sp_markers[[i]], features = RNA_modules[c], ctrl = 100, name = levels(MPNST_non_CC)[c])
    }
  }
  
  for (i in 1:length(samples)) {
    print(samples[i])
    DefaultAssay(MPNST_sp_markers[[i]]) <- "Spatial"
    
    png(filename = paste0("MPNST_dimplot_",samples[i],"_coarse_module_1_to_6.png"), width = 4000, height = 2000, res = 200)
    print(FeaturePlot(MPNST_sp_markers[[i]], features = paste0(levels(coarse.cluster.ids)[1:6],1), label = TRUE, ncol = 2))
    dev.off()
    
    png(filename = paste0("MPNST_dimplot_",samples[i],"_coarse_module_7_to_12.png"), width = 4000, height = 2000, res = 200)
    print(FeaturePlot(MPNST_sp_markers[[i]], features = paste0(levels(coarse.cluster.ids)[7:12],1), label = TRUE, ncol = 2))
    dev.off()
    
    png(filename = paste0("MPNST_spatial_",samples[i],"_coarse_module_1_to_6.png"), width = 4000, height = 2000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp_markers[[i]], features = paste0(levels(coarse.cluster.ids)[1:6],1), ncol = 2, alpha = c(0.1,1)) + theme(legend.position = "right"))
    dev.off()
    
    png(filename = paste0("MPNST_spatial_",samples[i],"_coarse_module_7_to_12.png"), width = 4000, height = 2000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp_markers[[i]], features = paste0(levels(coarse.cluster.ids)[7:12],1), ncol = 2, alpha = c(0.1,1)) + theme(legend.position = "right"))
    dev.off()
    DefaultAssay(MPNST_sp_markers[[i]]) <- "SCT"
  }
  setwd(paste0(output.dir))
}  

#Try Seurat's assign
if (F) {
  # DefaultAssay(MPNST_non_CC) <- "RNA"
  # png(filename = paste0("MPNST_C_non_CC_umap_SUZ12.png"), width = 2000, height = 2000, res = 200)
  # FeaturePlot(MPNST_non_CC, features = c("SUZ12"))
  # dev.off()

  #Perform anchoring  
  #Output for "P"   # Performing PCA on the provided reference using 1967 features as input.
  # Found 5894 anchors  # Retained 3577 anchors
  # MPNST_predictions <- pbmclapply(MPNST_sp, function(x) {
  # 
  # }, mc.cores = 8)
  # 
  MPNST_predictions <- list()
  MPNST_non_CC_anchors <- list()
  for (i in 1:length(samples)) {
    MPNST_non_CC_anchors[[i]] <- FindTransferAnchors(reference = MPNST_non_CC, query = MPNST_sp_markers[[i]], normalization.method = "SCT")
    MPNST_predictions[[i]] <- TransferData(anchorset = MPNST_non_CC_anchors[[i]], refdata = MPNST_non_CC$cell_type, prediction.assay = TRUE, 
                 weight.reduction = MPNST_sp_markers[[i]][["pca"]])
  }
  ###***NOTE***Seurat does like underscores so names changed to - e.g "Malignant_R5"
  saveRDS(MPNST_predictions, "MPNST_sp_predictions.rds")
  saveRDS(MPNST_non_CC_anchors, "MPNST_non_CC_anchors.rds")
  
  for (i in 1:length(samples)) {
    MPNST_sp_markers[[i]][["predictions"]] <- MPNST_predictions[[i]]

    DefaultAssay(MPNST_sp_markers[[i]]) <- "predictions"
    png(filename = paste0("MPNST_cell_type_",samples[i],"_pt1.png"), width = 4000, height = 4000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp_markers[[i]], features = c("Malignant-R5", "Malignant-R2", "Malignant-R4-2",
                                                   "Malignant-R3", "Malignant-R1-1", "Malignant-P",
                                                   "Malignant-R4-1", "Monocyte", "Malignant-R1-2"), ncol = 3, crop = TRUE))
    dev.off()
    
    png(filename = paste0("MPNST_cell_type_",samples[i],"_pt2.png"), width = 4000, height = 4000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp_markers[[i]], features = c("Malignant-Ribo", "Malignant-G2MS-R2R3", "Malignant-G2MS-R1R5",
                                                   "Malignant-G2MS-R4", "T cell", "Malignant-R4-3",
                                                   "Endothelial", "Skeletal Muscle"), ncol = 3, crop = TRUE))
    dev.off()
  }
}

####################################################################################################################################
###Integrate with harmony scRNA data (NOT RERUN)
####################################################################################################################################
if (F) {
  system(paste0("mkdir -p ", output.dir, "Harmony"))
  MPNST_non_CC <- readRDS(paste0(scrna.dir, "../C_harmony/MPNST_C_non_CC_har.rds"))
  new.cluster.ids <- c("Tumour_1","Tumour_2","Tumour_G2MS","Tumour_3","Tumour_4",
                       "Tumour_5","Macrophage","Tumour_6","Tumour_7","Tumour_Ribo",
                       "Tumour_8","T_Cell","Endothelial", "Skeletal_Muscle")
  names(new.cluster.ids) <- levels(MPNST_non_CC)
  MPNST_non_CC <- RenameIdents(MPNST_non_CC, new.cluster.ids)
  MPNST_non_CC$cell_type <- MPNST_non_CC@active.ident
  
  #Load in cluster markers
  MPNST_non_CC_cluster_markers <- readRDS(paste0(scrna.dir, "../C_harmony/MPNST_C_non_CC_cluster_markers.rds"))
  ##Make table with more info on top 10 genes
  top10 <- as_tibble(MPNST_non_CC_cluster_markers) %>% group_by(cluster) %>% arrange(cluster, -avg_logFC) %>% top_n(n = 10, wt = avg_logFC)
  
  RNA_modules <- lapply(1:length(levels(MPNST_non_CC)), function(c) {
    top10 %>% filter(cluster == c-1) %>% pull(gene)
  })
  
  for (c in 1:length(levels(MPNST_non_CC))) {
    print(paste0("Adding module score for ", levels(MPNST_non_CC)[c], " cluster"))
    MPNST_non_CC <- AddModuleScore(MPNST_non_CC, features = RNA_modules[c], ctrl = 100, name = levels(MPNST_non_CC)[c])
  }
  saveRDS(MPNST_non_CC, "Harmony/MPNST_C_non_CC_module_scores.rds")
  #Note module names shouldn't have spaces so renamed to T.cell and Skeletal.Muscle
  
  for (c in new.cluster.ids) {
    DefaultAssay(MPNST_non_CC) <- "RNA"
    png(filename = paste0("Harmony/MPNST_C_non_CC_module_",c,".png"), width = 2000, height = 2000, res = 200)
    print(FeaturePlot(MPNST_non_CC, features = paste0(c,1)))
    dev.off()
  }
  
  #Add module score to spatial
  for (i in 1:length(samples)) {
    print(samples[i])
    for (c in 1:length(levels(MPNST_non_CC))) {
      print(paste0("Adding module score for ", levels(MPNST_non_CC)[c], " cluster"))
      MPNST_sp_markers[[i]] <- AddModuleScore(MPNST_sp_markers[[i]], features = RNA_modules[c], ctrl = 100, name = levels(MPNST_non_CC)[c])
    }
  }
  
  for (i in 1:length(samples)) {
    print(samples[i])
    DefaultAssay(MPNST_sp_markers[[i]]) <- "Spatial"
    
    png(filename = paste0("Harmony/MPNST_dimplot_",samples[i],"_module_1_to_9.png"), width = 4000, height = 2000, res = 200)
    print(FeaturePlot(MPNST_sp_markers[[i]], features = paste0(new.cluster.ids[1:9],1), label = TRUE, ncol = 3))
    dev.off()
    
    png(filename = paste0("Harmony/MPNST_dimplot_",samples[i],"_module_10_to_14.png"), width = 4000, height = 2000, res = 200)
    print(FeaturePlot(MPNST_sp_markers[[i]], features = paste0(new.cluster.ids[10:14],1), label = TRUE, ncol = 3))
    dev.off()
    
    png(filename = paste0("Harmony/MPNST_spatial_",samples[i],"_module_1_to_9.png"), width = 2000, height = 2000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp_markers[[i]], features = paste0(new.cluster.ids[1:9],1), ncol = 3, alpha = c(0.1,1)) + theme(legend.position = "right"))
    dev.off()
    
    png(filename = paste0("Harmony/MPNST_spatial_",samples[i],"_module_10_to_14.png"), width = 2000, height = 2000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp_markers[[i]], features = paste0(new.cluster.ids[10:14],1), ncol = 3, alpha = c(0.1,1)) + theme(legend.position = "right"))
    dev.off()
    DefaultAssay(MPNST_sp_markers[[i]]) <- "SCT"
  }
}

####################################################################################################################################
###13/10/20 Try to add genotype data
####################################################################################################################################
if (F) {
  MPNST_sp <- readRDS("MPNST_sp.rds")
  both_haplo_count_combined_by_cell <- readRDS("/srv/shared/vanloo/home/hyan/mpnst/10X_spatial/results/Genotype_SNPs/both_haplo_count_combined_by_cell.rds")
  
  for (i in 1:length(samples)) {
    slide_barcodes <- enframe(MPNST_sp[[i]]@assays[["Spatial"]]@counts@Dimnames[[2]], name = NULL, value = "Barcode")
    haplo_count_pct <- both_haplo_count_combined_by_cell %>% filter(str_detect(Barcode, regions[i])) %>% 
      mutate(Haplo_Percent = Haplo_1_Count/(Haplo_1_Count+Haplo_2_Count)) %>% 
      mutate(Barcode = paste0(gsub(paste0(regions[i], "_"),"",Barcode), "-1")) %>%
      select(Barcode, Haplo_Percent)
    haplo_count_pct_barcodes <- left_join(slide_barcodes, haplo_count_pct, by = "Barcode")
    MPNST_sp[[i]] <- AddMetaData(MPNST_sp[[i]], pull(haplo_count_pct_barcodes, Haplo_Percent), col.name = "Haplo_Percent")
    
    png(filename = paste0("MPNST_haplo1_pct_",samples[i],".png"), width = 4000, height = 4000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp[[i]], features = "Haplo_Percent", alpha = c(0.1,1)) + ggtitle("Percent of all counts which are Haplotype 1"))
    dev.off()
    # Tried to change limits but didn't work for SpatialFeaturePlot
    # FeaturePlot(MPNST_sp[[1]], features = "TTN") + scale_color_gradientn(colours = c("red","green"), limits = c(0,1))
  }  
  
  for (i in 1:length(samples)) {
    slide_barcodes <- enframe(MPNST_sp[[i]]@assays[["Spatial"]]@counts@Dimnames[[2]], name = NULL, value = "Barcode")
    haplo_count_ratio <- both_haplo_count_combined_by_cell %>% filter(str_detect(Barcode, regions[i])) %>% 
      mutate(Haplo_Ratio = Haplo_1_Count/Haplo_2_Count) %>% 
      mutate(Barcode = paste0(gsub(paste0(regions[i], "_"),"",Barcode), "-1")) %>%
      select(Barcode, Haplo_Ratio)
    haplo_count_ratio_barcodes <- left_join(slide_barcodes, haplo_count_ratio, by = "Barcode")
    MPNST_sp[[i]] <- AddMetaData(MPNST_sp[[i]], pull(haplo_count_ratio_barcodes, Haplo_Ratio), col.name = "Haplo_Ratio")
    
    png(filename = paste0("MPNST_haplo_ratio_",samples[i],".png"), width = 4000, height = 4000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp[[i]], features = "Haplo_Ratio", alpha = c(0.1,1)) + ggtitle("Percent of all counts which are Haplotype 1"))
    dev.off()
    # Tried to change limits but didn't work for SpatialFeaturePlot
    # FeaturePlot(MPNST_sp[[1]], features = "TTN") + scale_color_gradientn(colours = c("red","green"), limits = c(0,1))
  }  
}

####################################################################################################################################
###03/10/22 Try to add hallmarks/cancer states data
####################################################################################################################################
if (!file.exists("MPNST_sp_cancer_states.rds")) {
  MPNST_sp <- readRDS("MPNST_sp.rds")
 
  #Add Cancer cell states from Barkley et.al Nature Genetics 2022
  cancer_states_csv <- read.delim("../../../10X_RNA/input/41588_2022_1141_MOESM6_ESM.csv", sep = ",", header = F, stringsAsFactors = F)
  cancer_states <- lapply(1:ncol(cancer_states_csv), function(state){
    return(cancer_states_csv[4:nrow(cancer_states_csv),state][cancer_states_csv[4:nrow(cancer_states_csv),state] != ""])
  })
  names(cancer_states) <- cancer_states_csv[3,]
  MPNST_sp_states <- MPNST_sp
  
  for (s in 1:length(MPNST_sp_states)) {
    # for (i in 1:length(cancer_states)) {
    #   MPNST_sp_states[[s]] <- AddModuleScore(MPNST_sp_states[[s]], features = cancer_states[i], assay = "Spatial", nbin = 20, name = paste0(names(cancer_states)[i],"_SP"), seed = 1)
    #   MPNST_sp_states[[s]] <- AddModuleScore(MPNST_sp_states[[s]], features = cancer_states[i], assay = "SCT", nbin = 20, name = paste0(names(cancer_states)[i],"_SCT"), seed = 1)
    # }
    #Plot cancer states
    png(filename = paste0("MPNST_spatial_cancer_states_RNA_", samples[s], ".png"), width = 4000, height = 4000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp_states[[s]], features = paste0(names(cancer_states), "_SP1"), ncol = 4))
    dev.off()
    
    png(filename = paste0("MPNST_spatial_cancer_states_SCT_", samples[s], ".png"), width = 4000, height = 4000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp_states[[s]], features = paste0(names(cancer_states), "_SCT1"), ncol = 4))
    dev.off()
  }
  
  #Add hallmarks of ITH from Gavish et al. from Nature 2023
  hallmarks_ITH_xls <- read_excel(path = paste0("../../../10X_RNA/input/41586_2023_6130_MOESM6_ESM.xlsx"), sheet = "Cancer MPs", col_names = T)
  hallmarks_ITH <- lapply(colnames(hallmarks_ITH_xls), function(state){
    return(hallmarks_ITH_xls %>% pull(state))
  })
  names(hallmarks_ITH) <- gsub("Cylce", "Cycle", gsub(" " , "_", gsub("  ", " ", gsub("\\(|\\)|", "", gsub("-|/", "_", gsub(" - ", "_", colnames(hallmarks_ITH_xls)))))))
  
  for (s in 1:length(MPNST_sp_states)) {
    for (i in 1:length(hallmarks_ITH)) {
      print(paste0("Adding module score for ", names(hallmarks_ITH)[i]))
      MPNST_sp_states[[s]] <- AddModuleScore(MPNST_sp_states[[s]], features = hallmarks_ITH[i], assay = "Spatial", nbin = 20, name = paste0(names(hallmarks_ITH)[i],"_SP"), seed = 1)
      MPNST_sp_states[[s]] <- AddModuleScore(MPNST_sp_states[[s]], features = hallmarks_ITH[i], assay = "SCT", nbin = 20, name = paste0(names(hallmarks_ITH)[i],"_SCT"), seed = 1)
    }
    #Plot hallmarks of ITH
    png(filename = paste0("MPNST_spatial_hallmarks_ITH_RNA_", samples[s], ".png"), width = 3000, height = 2000, res = 100)
    print(SpatialFeaturePlot(MPNST_sp_states[[s]], features = paste0(names(hallmarks_ITH), "_SP1")[1:40], ncol = 8))
    dev.off()
    
    png(filename = paste0("MPNST_spatial_hallmarks_ITH_SCT_", samples[s], ".png"), width = 3000, height = 2000, res = 100)
    print(SpatialFeaturePlot(MPNST_sp_states[[s]], features = paste0(names(hallmarks_ITH), "_SCT1")[1:40], ncol = 8))
    dev.off()
  }
  saveRDS(MPNST_sp_states, "MPNST_sp_cancer_states.rds") 
} else {
  MPNST_sp_states <- readRDS("MPNST_sp_cancer_states.rds")
}

#Pick region and plot
sample = samples[1]
cancer_state = "Stress"
png(filename = paste0("MPNST_spatial_",cancer_state,"_RNA_", sample, ".png"), width = 4000, height = 4000, res = 200)
SpatialFeaturePlot(MPNST_sp_states[[which(samples == sample)]], features = paste0(names(cancer_state), "_SCT1"), ncol = 4)
dev.off()

####################################################################################################################################
###28/08/23 Check chr10 expression for R4_4 subclone
####################################################################################################################################
if (F) {
  gene_positions <- read.delim("../../../../ref_files/inferCNV/gencode_v21_gen_pos.complete.clean.txt", header = F, stringsAsFactors = F)
  colnames(gene_positions) <- c("gene", "chr", "start", "end")
  
  MPNST_sp <- readRDS("MPNST_sp.rds")
  
  MPNST_sp_count_mtx <- MPNST_sp[[6]]@assays[["Spatial"]]@counts %>% as.matrix() %>% t()
  MPNST_SCT_count_mtx <- MPNST_sp[[6]]@assays[["SCT"]]@counts %>% as.matrix() %>% t()
  
  chr_10_genes <- colnames(MPNST_SCT_count_mtx)[colnames(MPNST_SCT_count_mtx) %in% (gene_positions %>% filter(chr == "chr10") %>% pull(gene))]
  chr_11_genes <- colnames(MPNST_SCT_count_mtx)[colnames(MPNST_SCT_count_mtx) %in% (gene_positions %>% filter(chr == "chr11") %>% pull(gene))]
  total_exp <- rowSums(MPNST_SCT_count_mtx)
  chr10_exp <- rowSums(MPNST_SCT_count_mtx[,chr_10_genes])/total_exp*100
  chr11_exp <- rowSums(MPNST_SCT_count_mtx[,chr_11_genes])/total_exp*100
  MPNST_sp[[6]] <- AddMetaData(MPNST_sp[[6]], chr10_exp, col.name = "chr10_exp")
  MPNST_sp[[6]] <- AddMetaData(MPNST_sp[[6]], chr11_exp, col.name = "chr11_exp")
  
  # png(filename = paste0("MPNST_R4_chr10_exp.png"), width = 2000, height = 2000, res = 200)
  pdf(file = paste0("MPNST_R4_chr10_exp.pdf"), width = 7, height = 7)
  print(SpatialFeaturePlot(MPNST_sp[[6]], features = "chr10_exp", alpha = c(0.1, 1)))
  dev.off()
  # png(filename = paste0("MPNST_R4_chr11_exp.png"), width = 2000, height = 2000, res = 200)
  pdf(file = paste0("MPNST_R4_chr11_exp.pdf"), width = 7, height = 7)
  print(SpatialFeaturePlot(MPNST_sp[[6]], features = "chr11_exp", alpha = c(0.1, 1)))
  dev.off()
}

####################################################################################################################################
###14/10/20 Try to add STARCH data (NOT RERUN)
####################################################################################################################################
if (F) {
  MPNST_sp <- readRDS("MPNST_sp.rds")
  starch.output.dir = paste0("/srv/shared/vanloo/home/hyan/mpnst/10X_spatial/results/STARCH/", samples)
  
  for (i in 1:length(samples)) {
    starch.labels = read.csv(paste0(starch.output.dir[i], "/labels_", samples[i], ".csv"), col.names = c("Spot", "Cluster"))
    MPNST_sp[[i]] <- AddMetaData(MPNST_sp[[i]], starch.labels$Cluster, col.name = "Starch_Cluster")
    
    png(filename = paste0("MPNST_starch_cluster_",samples[i],".png"), width = 4000, height = 4000, res = 200)
    print(SpatialFeaturePlot(MPNST_sp[[i]], features = "Starch_Cluster", alpha = c(0.1,1)))
    dev.off()
  }
}
