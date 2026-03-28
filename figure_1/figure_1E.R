### title: Plot the copy number profiles of single cells. 

# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"

# Define input directory with data and output directory to save the figure
input.dir <- paste0(zenodo.dir, "data/scDNA/")
output.dir <- paste0(zenodo.dir, "results/figure_1/")
setwd(output.dir)

####################################################################################################################################
### Part 0: Load libraries and data
####################################################################################################################################

### Load libraries
library(tidyverse)
library(ComplexHeatmap)
library(grid)
library(gridExtra)
library(magrittr)

### Load cluster ID information from CN clusters (k-means)
K <- 22
scDNA_CN_cluster_ids <- readRDS(paste0(input.dir, "CN_profiles/MPNST_all_k_means_K",K,"_clusters.rds"))

### Load ASCAT.sc CN profiles and probe information
gamma = 5
asCN_chr_probes <- readRDS(paste0(input.dir, "CN_profiles//MPNST_all_asCN_chr_probes_mpcf_", gamma,".rds"))

scDNA_asCN_mtx <- readRDS(paste0(input.dir,"CN_profiles/MPNST_DLP_asfree_asCN_mtx.rds"))

####################################################################################################################################
### Part 1: Prepare metadata for 10X and DLP
####################################################################################################################################

K = 22
run <- "22_kmeans_10X_DLP"
remove_kmeans2 = T
scDNA_CN_cluster_ids <- scDNA_CN_cluster_ids
scDNA_normal_ids <- scDNA_CN_cluster_ids[[1]]
scDNA_CN_cluster_orig_names <- (1:length(scDNA_CN_cluster_ids))[-1] #remove normal
scDNA_CN_cluster_ids <- scDNA_CN_cluster_ids[-1] #remove normal
scDNA_CN_cluster_orig_names <- scDNA_CN_cluster_orig_names[lengths(scDNA_CN_cluster_ids) > 5] #remove small clusters
scDNA_CN_cluster_ids <- scDNA_CN_cluster_ids[lengths(scDNA_CN_cluster_ids) > 5] #remove small clusters
if (remove_kmeans2) {
  scDNA_CN_cluster_orig_names <- scDNA_CN_cluster_orig_names[-1] #remove k means cluster 2
  scDNA_CN_cluster_ids <- scDNA_CN_cluster_ids[-1] #remove k means cluster 2
  scDNA_CN_clusters_origins <- c("R4", "P", "R1", "R4", 
                                 "R2", "R4", "R5", "R5", "R3",
                                 "R2", "R2", "P", "R5", "R1",
                                 "R1", "R1", "R5", "P")
} else {
  scDNA_CN_clusters_origins <- c("R1", "R4", "P", "R1", "R4", 
                                 "R2", "R4", "R5", "R5", "R3",
                                 "R2", "R2", "P", "R5", "R1",
                                 "R1", "R1", "R5", "P")
}

names(scDNA_CN_cluster_ids) <- ave(scDNA_CN_clusters_origins, scDNA_CN_clusters_origins, FUN = function(i) paste0(i, "_", seq_along(i)))

####################################################################################################################################
### Part 2: A heatmap function to plot allele-specific CN profiles of single cells
####################################################################################################################################

sc_asCN_heatmap3 <- function(CN_mtx, hclust = F, km = NULL, row_split = NULL, probes, row_ann = NULL, column_title = "Single Cell Copy Number Heatmap", title = "Allele-specific copy number state", 
                             colour_scheme = c("royalblue3", "skyblue2", "grey80", "white", "gold1", "khaki1", "darkorange3", "darkorange1", "orange", "red4", "red", "orangered2", "purple4")) {
  show_raw_dend = T
  if (!is.null(km)) {
    hclust = T
    show_raw_dend = F
  } #Cluster each kmeans cluster but don't show dendrogram
  CN_mtx[CN_mtx<0] <- 0 #set min CN value
  CN_mtx[CN_mtx>60] <- 60 #set max CN value
  names(colour_scheme) <- c(0, 10, 20, 21, 30, 31, 40, 41, 42, 50, 51, 52, 60)#Colours
  
  ### The actual change is in defining CN_labels
  #CN_labels <- c("0+0 ", "1+0 ", "2+0 ", "1+1 ", "3+0 ", "2+1 ", "4+0 ", "3+1 ", "2+2 ", "5+0 ", "4+1 ", "3+2 ", ">6 ")[which(names(colour_scheme) %in% sort(unique(as.numeric(CN_mtx))))] #Get labels in legend
  CN_at <- intersect(as.numeric(names(colour_scheme)),
                     sort(unique(as.numeric(CN_mtx))))
  
  CN_label_map <- c(
    "0"="0+0 ", "10"="1+0 ", "20"="2+0 ", "21"="1+1 ",
    "30"="3+0 ", "31"="2+1 ", "40"="4+0 ", "41"="3+1 ",
    "42"="2+2 ", "50"="5+0 ", "51"="4+1 ", "52"="3+2 ",
    "60"=">6 "
  )
  CN_labels <- unname(CN_label_map[as.character(CN_at)])
  
  ha_column = HeatmapAnnotation(odd = anno_empty(border = F),
                                even = anno_empty(border = F)) #Chromosome annotation at bottom
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(matrix = CN_mtx, cluster_rows = hclust, row_km = km, row_split = row_split, cluster_row_slices = FALSE, show_row_dend = show_raw_dend, use_raster = FALSE,
                                               row_title_rot = 0, cluster_columns = F, show_row_names = F, show_column_names = F, row_title_gp = grid::gpar(fontsize = 24),
                                               heatmap_legend_param = list(at = CN_at, labels = CN_labels, title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 20), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), border = "black", nrow=1),
                                               bottom_annotation = ha_column, left_annotation = row_ann, col = colour_scheme, column_title = column_title, name = title),
                       heatmap_legend_side="bottom", annotation_legend_side="bottom")
  probes <- c(0, probes) #Adds posiiton zero for start of first chromosome
  
  decorate_heatmap_body(heatmap = title, {
    for (k in 1:(length(probes))) {
      grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1), gp=gpar(col="black", lty = 1, lwd = 1.5))
    }
  }) #Chr lines
  if (!is.null(row_split)) {
    for (i in 2:length(unique(row_split))) {
      decorate_heatmap_body(heatmap = title, row_slice = i, {
        for (k in 1:(length(probes))) {
          grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
        }
      }) #Chr lines
    }
  }
  if (!is.null(km)) {
    for (i in 2:K) {
      decorate_heatmap_body(heatmap = title, row_slice = i, {
        for (k in 1:(length(probes))) {
          grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
        }
      }) #Chr lines
    }
  }
  decorate_annotation(annotation = "odd", {
    for (k in 2:(length(probes))) {
      if (k%%2 == 0) {grid.text(names(probes)[k], x=((probes[k-1]+((probes[k]-probes[k-1])/2))/ncol(CN_mtx)) + 0.001, y=0.5, just = "left", gp=gpar(cex = 3))}
    }
  }) #Chr numbers
  decorate_annotation(annotation = "even", {
    for (k in 2:(length(probes))) {
      if (k%%2 == 1) {grid.text(names(probes)[k], x=((probes[k-1]+((probes[k]-probes[k-1])/2))/ncol(CN_mtx)) + 0.001, y=0.5, just = "left", gp=gpar(cex = 3))}
    }
  }) #Chr numbers
}

####################################################################################################################################
### Part 3: Plot the large heatmap based on asCN profiles
####################################################################################################################################

### We require:
# Objects: (1) scDNA_asCN_mtx (2) scDNA_CN_cluster_ids (3) asCN_chr_probes
# Function: sc_asCN_heatmap3

### Create row annotation for region and tech

named_clusters <- unlist(lapply(1:length(scDNA_CN_cluster_ids), function(i) return(rep(names(scDNA_CN_cluster_ids)[i], length(scDNA_CN_cluster_ids[[i]]))))) %>% set_names(unlist(scDNA_CN_cluster_ids))

ha_row = rowAnnotation(
  df = data.frame(Region = gsub("_.*", "", rownames(scDNA_asCN_mtx[names(named_clusters),])), 
                  Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_asCN_mtx[names(named_clusters),]))),
  col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"), 
             Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")),
  annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 22), 
                                               grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1), 
                                 Tech = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 22), grid_height = unit(0.8, "cm"), 
                                             grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)

### Create the plot
png(filename = paste0("figure_1E.png"), width = 4000, height = 4000, res = 200)
print(sc_asCN_heatmap3(CN_mtx = scDNA_asCN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, row_ann = ha_row, column_title = NA, probes = asCN_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X"))))
dev.off()


