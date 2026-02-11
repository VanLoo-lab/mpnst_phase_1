options(bitmapType='cairo') #to solve plotting issue on CAMP
library(tidyverse, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(ggpubr, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(gridExtra, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(VariantAnnotation, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(liftOver, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(copynumber, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(pheatmap, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(pbapply, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(pbmcapply, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(dendextend, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(grid, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(gridExtra, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")

source("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/code/10X_DLP_metadata.R")

bulk_samples = c("VER236A1", "VER236A2", "VER236A3", "VER236A4", "VER236A5", "VER236A6")
names(bulk_samples) <- samples

chrom_numbers <- c(1:23)
chrom_names <- c(1:22,"X")

output.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/results/Co-occurence_SNVs/")
system(paste0("mkdir -p ", output.dir))
setwd(output.dir)

if (F) {
  #Load in 10X and DLP genotyped SNVs
  bulk_mt_10X_counts_SNV_vec_cell <- readRDS("../../../10X_DNA/results/Genotype_SNVs/MPNST_scDNA_mt_counts_SNV_vec_cell.rds")
  names(bulk_mt_10X_counts_SNV_vec_cell) <- unlist(barcodes)
  
  bulk_mt_DLP_counts_SNV_vec_cell <- readRDS("../../../DLP_plus/results/Genotype_SNVs_CICC/MPNST_scDNA_mt_counts_SNV_vec_cell.rds")
  names(bulk_mt_DLP_counts_SNV_vec_cell) <- all_DLP_barcodes
  
  bulk_mt_counts_SNV_vec_cell <- c(bulk_mt_10X_counts_SNV_vec_cell, bulk_mt_DLP_counts_SNV_vec_cell)
  
  mt_counts_SNV_vec_cell <- lapply(1:length(bulk_mt_counts_SNV_vec_cell), function(c) {
    return(bulk_mt_counts_SNV_vec_cell[[c]] %>% arrange(chr,pos) %>% pull(Present))
  })
  
  #Calculate co-occurence matrix
  mt_counts_SNV_mtx_cell <- matrix(rep(0, length(mt_counts_SNV_vec_cell[[1]])^2), 
                                   nrow = length(mt_counts_SNV_vec_cell[[1]]), ncol = length(mt_counts_SNV_vec_cell[[1]]))
  for (i in 1:length(mt_counts_SNV_vec_cell)) { #for each cell
    if (i %% 100 == 0) {
      print(i)
    }
    mt_counts_SNV_mtx_cell <- mt_counts_SNV_mtx_cell + mt_counts_SNV_vec_cell[[i]] %*% t(mt_counts_SNV_vec_cell[[i]]) #add matrix together
  }
  
  saveRDS(mt_counts_SNV_mtx_cell, "MPNST_scDNA_mt_counts_SNV_mtx_cell.rds")
}

#Plot
if (F) {
  #Reorder clusters
  mt_counts_SNV_mtx_cell <- readRDS("MPNST_scDNA_mt_counts_SNV_mtx_cell.rds")
  SNVs_clusters <- readRDS("../../../DLP_plus/results/Genotype_SNVs_CICC/MPNST_SNVs_clusters.rds")
  rownames(mt_counts_SNV_mtx_cell) <- paste0(do.call(rbind, SNVs_clusters) %>% arrange(chr, pos) %>% pull(chr), "_", do.call(rbind, SNVs_clusters) %>% arrange(chr, pos) %>% pull(pos))
  colnames(mt_counts_SNV_mtx_cell) <- paste0(do.call(rbind, SNVs_clusters) %>% arrange(chr, pos) %>% pull(chr), "_", do.call(rbind, SNVs_clusters) %>% arrange(chr, pos) %>% pull(pos))
  
  SNVs_clusters_reorder <- SNVs_clusters
  # SNVs_clusters_reorder <- SNVs_clusters[c(2,6,16,15,17,4,7,3,1,10,12,13,5,11,9,8,14)]
  mt_counts_SNV_mtx_cell <- mt_counts_SNV_mtx_cell[paste0(do.call(rbind, SNVs_clusters_reorder)$chr, "_", do.call(rbind, SNVs_clusters_reorder)$pos), 
                                                   paste0(do.call(rbind, SNVs_clusters_reorder)$chr, "_", do.call(rbind, SNVs_clusters_reorder)$pos)]
  
  #Filter to set threshold
  n_cells = 1
  mt_counts_SNV_mtx_filter <- mt_counts_SNV_mtx_cell
  mt_counts_SNV_mtx_filter[mt_counts_SNV_mtx_filter>n_cells] <- n_cells
  mt_counts_SNV_mtx_filter <- 1-(mt_counts_SNV_mtx_filter/max(mt_counts_SNV_mtx_filter))
  
  # cluster_names <- cluster_names_orig[clusters_reorder]
  cluster_names <- names(SNVs_clusters)
  cluster_colours_reorder <- c("black", "#F8766D", "pink", "coral", "#B79F00", "gold", "#F564E3", "mediumpurple1",
                               "darkslategray4", "darkgreen", "#00BA38", "chartreuse", "#619CFF", "deepskyblue", "#00BFC4", "cyan")
  
  add_box <- T
  if (add_box) {
    add_box_suffix <- "_add_box"
  } else {
    add_box_suffix <- ""
  }
  # png(filename = paste0("MPNST_scDNA_SNVs_co-clustering_raster_filter_",n_cells,"_reorder_no_artefact",add_box_suffix,".png"), width = 4000, height = 4000, res = 200)
  pdf(file = paste0("MPNST_scDNA_SNVs_co-clustering_raster_filter_",n_cells,"_reorder_no_artefact",add_box_suffix,".pdf"), width = 7, height = 7)
  plot(0:1, 0:1, type="n", xlim = c(-0.05,1.05), ylim = c(0,1.05), axes = F, xlab = NA, ylab = NA)
  rasterImage(mt_counts_SNV_mtx_filter, xleft = 0, ybottom = 0, xright = 1, ytop = 1, interpolate = F)
  for (r in 1:length(SNVs_clusters_reorder)) {
    if (r == 9) { #adjusts label position for a small cluster
      adj_x = 0.002
      adj_y = 0.005
    } else {
      adj_x = 0.005
      adj_y = 0
    }
    #Plot boxes
    rect(xleft = (sum(unlist(lapply(SNVs_clusters_reorder[1:r], nrow)))-nrow(SNVs_clusters_reorder[[r]]))/ncol(mt_counts_SNV_mtx_filter),
         ybottom = 1.01,
         xright = sum(unlist(lapply(SNVs_clusters_reorder[1:r], nrow)))/ncol(mt_counts_SNV_mtx_filter),
         ytop = 1.02,
         col = cluster_colours_reorder[r])
    rect(xleft = -0.02,
         ybottom = 1-(sum(unlist(lapply(SNVs_clusters_reorder[1:r], nrow)))/ncol(mt_counts_SNV_mtx_filter)),
         xright = -0.01,
         ytop = 1-(sum(unlist(lapply(SNVs_clusters_reorder[1:r], nrow)))-nrow(SNVs_clusters_reorder[[r]]))/ncol(mt_counts_SNV_mtx_filter),
         col = cluster_colours_reorder[r])
    
    #Plot horizontal labels
    text(x = (sum(unlist(lapply(SNVs_clusters_reorder[1:r], nrow)))-(nrow(SNVs_clusters_reorder[[r]])/2))/ncol(mt_counts_SNV_mtx_filter)+adj_x,
         y = 1.05,
         labels = cluster_names[r],
         adj = c(0.3,-0.1),
         srt = 90,
         cex = 1.2)    
    #Plot vertical labels
    text(y = 1-(sum(unlist(lapply(SNVs_clusters_reorder[1:r], nrow)))-(nrow(SNVs_clusters_reorder[[r]])/2))/ncol(mt_counts_SNV_mtx_filter)+adj_y,
         x = -0.05,
         labels = cluster_names[r],
         adj = c(0.7,0.5),
         cex = 1.2)
    if (add_box) {
      rect(xleft = 0,
           ybottom = 1-(sum(unlist(lapply(SNVs_clusters_reorder[1:1], nrow)))/ncol(mt_counts_SNV_mtx_filter)),
           xright = 1,
           ytop = 1,
           border = "black",
           lwd = 3)
      rect(xleft = 0,
           ybottom = 0,
           xright = sum(unlist(lapply(SNVs_clusters_reorder[1:1], nrow)))/ncol(mt_counts_SNV_mtx_filter),
           ytop = 1,
           border = "black",
           lwd = 3)
      rect(xleft = (sum(unlist(lapply(SNVs_clusters_reorder[1:1], nrow))))/ncol(mt_counts_SNV_mtx_filter),
           ybottom = 1-(sum(unlist(lapply(SNVs_clusters_reorder[1:1], nrow)))/ncol(mt_counts_SNV_mtx_filter)),
           xright = sum(unlist(lapply(SNVs_clusters_reorder[1:3], nrow)))/ncol(mt_counts_SNV_mtx_filter),
           ytop = 1-(sum(unlist(lapply(SNVs_clusters_reorder[1:3], nrow))))/ncol(mt_counts_SNV_mtx_filter),
           border = "#F8766D",
           lwd = 3)
      rect(xleft = (sum(unlist(lapply(SNVs_clusters_reorder[1:3], nrow))))/ncol(mt_counts_SNV_mtx_filter),
           ybottom = 1-(sum(unlist(lapply(SNVs_clusters_reorder[1:3], nrow)))/ncol(mt_counts_SNV_mtx_filter)),
           xright = sum(unlist(lapply(SNVs_clusters_reorder[1:8], nrow)))/ncol(mt_counts_SNV_mtx_filter),
           ytop = 1-(sum(unlist(lapply(SNVs_clusters_reorder[1:8], nrow))))/ncol(mt_counts_SNV_mtx_filter),
           border = "coral",
           lwd = 3)
      rect(xleft = (sum(unlist(lapply(SNVs_clusters_reorder[1:8], nrow))))/ncol(mt_counts_SNV_mtx_filter),
           ybottom = 1-(sum(unlist(lapply(SNVs_clusters_reorder[1:8], nrow)))/ncol(mt_counts_SNV_mtx_filter)),
           xright = 1,
           ytop = 0,
           border = "darkslategray4",
           lwd = 3)
    }
  }
  dev.off()
}