### title: Cell type deconvolution of spatial transcriptomics data

# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"

# Define input directory with data and output directory to save the figure
input.dir <- paste0(zenodo.dir, "data/")
output.dir <- paste0(zenodo.dir, "results/figure_4/")
setwd(output.dir)


####################################################################################################################################
### Part 0: Load libraries and data
####################################################################################################################################

### Load libraries
library(tidyverse)
library(grid)
library(gridExtra)
library(Seurat)
library(ComplexHeatmap)

### Load scDNA data
scDNA.dir = paste0(input.dir, "scDNA/CN_profiles/")

all_probes <- readRDS(paste0(scDNA.dir, "MPNST_all_probes_mpcf_5.rds"))
all_chr_probes <- readRDS(paste0(scDNA.dir, "MPNST_all_chr_probes_mpcf_5.rds"))

### Load spRNA data
spRNA.dir <- paste0(input.dir, "spRNA/")

MPNST_sp_markers <- readRDS(paste0(spRNA.dir, "MPNST_sp_markers.rds"))  # spRNA markers
R1_ref_CN <- readRDS(paste0(spRNA.dir, "MPNST_R1_ref_CN.rds"))          # Load reference CN profile

### Load inferCNV results
inferCNV.dir <- paste0(input.dir,"spRNA/inferCNV/")
# loaded per sample in for loop below

### Load metadata
#source( paste0(CODEDIR, "metadata/MPNST_common_functions.R") )
#source( paste0(CODEDIR, "metadata/10X_DLP_metadata_minimal.R") )
#source( paste0(CODEDIR, "metadata/Visium_metadata_minimal.R") )

####################################################################################################################################
### Part 1: Define metadata
####################################################################################################################################

### Define color pallets
color_scale <- colorRampPalette(c("darkblue", "grey95", "darkred"))(100)

scRNA_mtx_breaks = seq(0.8,1.2, by = 0.004)
scRNA_adj_mtx_breaks = seq(-0.2,0.2, by = 0.004)

scDNA_CN_colour <- c("#00008B", "#7B7BC3", "#FFFFFF", "#FFCCCC", "#FF8080", "#FF3333", "#E60000", "#B31000", "#8B0000", "#660000", "#330000")
scDNA_CN_breaks = seq(-0.5,10.5, by = 1)

### Define function to discretize
disc_mtx <- function(mtx, breaks, labels) {
  d_mtx <- matrix(as.matrix(cut(mtx, breaks = breaks, labels = labels)) %>% as.numeric(), nrow(mtx))
  rownames(d_mtx) <- rownames(mtx)
  colnames(d_mtx) <- colnames(mtx)
  return(d_mtx)
}

### Define sample names
visium_samples = c("Pa", "R1", "R3", "Pb", "R2", "R4", "R5a", "R5b")
names(visium_samples) = c("VER683A1", "VER683A2", "VER683A3", "VER683A4", "VER683A5", "VER683A6", "VER683A7", "VER683A8")
# sample = "C"
color_scale <- colorRampPalette(c("darkblue", "white", "darkred"))(100)

combos <- combn(visium_samples[c(2:8,1)], 2)
pairs <- paste0(combos[1,], "_", combos[2,])

### Rename spRNA clusters to include sample name and "SP" suffix to avoid confusion with scDNA clusters in the following analysis
for (s in visium_samples) {
  MPNST_sp_markers[[which(visium_samples == s)]] <- RenameCells(MPNST_sp_markers[[which(visium_samples == s)]],
                                                                new.names = gsub(paste0("(................)-1"), paste0(s, "_\\1","SP"),colnames(MPNST_sp_markers[[which(visium_samples == s)]])))
  sp.cluster.ids <- paste0(s, "_",0:(length(levels(MPNST_sp_markers[[which(visium_samples == s)]]))-1))
  names(sp.cluster.ids) <- levels(MPNST_sp_markers[[which(visium_samples == s)]])
  MPNST_sp_markers[[which(visium_samples == s)]] <- RenameIdents(MPNST_sp_markers[[which(visium_samples == s)]], sp.cluster.ids)
}

### Define parameters for copy number gain and loss
discretize_breaks <- c(-Inf,-0.2,-0.05,0.05,0.2,Inf)
discretize_labels <- c(-2,-1,0,1,2)
discretize_colors = c("darkblue", "royalblue1", "grey95", "firebrick3", "darkred")

### Define parameters for clustering
rm_small_cluster <- 0 #Set to 50 or 0
adj_cluster <- F
high_res <- F
single_spot <- F
kd_centroid <- T
suffix <- paste0(ifelse(rm_small_cluster > 0, paste0("_rm_sub", rm_small_cluster), ""), ifelse(adj_cluster, "_adj", ""), ifelse(high_res, "_res", ""), ifelse(single_spot, "_single", ""))

if (high_res) {
  for (s in visium_samples) {
    MPNST_sp_markers[[which(visium_samples == s)]] <- FindClusters(MPNST_sp_markers[[which(visium_samples == s)]], resolution = 1.2, verbose = FALSE)
  }
}

####################################################################################################################################
### Part 2: Function to create complex heatmap
####################################################################################################################################

sc_totCN_heatmap2 <- function(CN_mtx, hclust = F, km = NULL, row_split = NULL, probes, row_ann = NULL, max_CN = 10, column_title = "Single Cell Copy Number Heatmap", title = "Total copy number state", 
                              colour_scheme = c("#00008B", "#7B7BC3", "#FFFFFF", "#FFCCCC", "#FF8080", "#FF3333", "#E60000", "#B31000", "#8B0000", "#660000", "#330000")) {
  show_raw_dend = T
  if (!is.null(km)) {
    hclust = T
    show_raw_dend = F
  } #Cluster each kmeans cluster but don't show dendrogram
  CN_mtx[CN_mtx<0] <- 0 #set min CN value
  CN_mtx[CN_mtx>max_CN] <- max_CN #set max CN value
  names(colour_scheme) <- 0:max_CN #Colours
  #CN_labels <- paste0(0:max_CN, " ")[which(names(colour_scheme) %in% sort(unique(as.numeric(CN_mtx))))] #Get labels in legend
  CN_at <- intersect(0:max_CN, sort(unique(as.numeric(CN_mtx))))
  CN_labels <- paste0(CN_at, " ")
  #Annotation
  ha_column = HeatmapAnnotation(odd = anno_empty(border = F),
                                even = anno_empty(border = F)) #Chromosome annotation at bottom
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(matrix = CN_mtx, cluster_rows = hclust, row_km = km, row_split = row_split, cluster_row_slices = FALSE, show_row_dend = show_raw_dend, 
                                               row_title_rot = 0, cluster_columns = F, show_row_names = F, show_column_names = F, row_title_gp = grid::gpar(fontsize = 24),
                                               heatmap_legend_param = list(at = CN_at, labels = CN_labels, title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 20), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), border = "black", nrow=1),
                                               use_raster = FALSE,
                                               bottom_annotation = ha_column, left_annotation = row_ann, col = colour_scheme, column_title = column_title, name = title),
                       heatmap_legend_side="bottom", annotation_legend_side="bottom")
  #Add chr lines
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
### Part 3: Plot inferCNV CN heatmap per sample (including clusters)
####################################################################################################################################


for (sample in pairs[1:7]) {
  #sample = pairs[1]

  segmented_mtx <- readRDS(paste0(inferCNV.dir,"MPNST_",sample,"_cluster_segmented_mtx",suffix,".rds"))
  phylos <- readRDS(paste0(inferCNV.dir,"MPNST_",sample,"_cluster_trees",suffix,".rds"))
  merge_cluster_ids <- readRDS(paste0(inferCNV.dir,"MPNST_", sample, "_clusters",suffix,".rds"))

  
  #Generate integer CN matrix
  #Convert segmented mtx to medicc input 
  medicc_mtx <- disc_mtx(segmented_mtx, breaks = discretize_breaks, labels = discretize_labels)
  
  #convert CN profiles
  medicc_totCN_mtx <- do.call(rbind, lapply(1:nrow(medicc_mtx), function(c) {
    return(R1_ref_CN %>% mutate(cn_a = ifelse(is.na(totCN + medicc_mtx[c,]), totCN, totCN + medicc_mtx[c,])) %>% pull(cn_a))
  }))
  rownames(medicc_totCN_mtx) <- rownames(medicc_mtx)
  
  medicc_totCN_plot_mtx <- do.call(rbind, rep(lapply(1:nrow(medicc_mtx), function(c) {
    CN <- R1_ref_CN %>% mutate(cn_a = ifelse(is.na(totCN + medicc_mtx[c,]), totCN, totCN + medicc_mtx[c,])) %>% pull(cn_a)
    return(rep(CN, as.list(R1_ref_CN$n.probes)))
  }), lapply(merge_cluster_ids, function(c) length(c)/10)))
  
  #Plot
  png(filename = paste0("figure_4B_",sample,".png"), width = 4000, height = 2000, res = 200)
  sc_totCN_heatmap2(medicc_totCN_plot_mtx, hclust = F, row_split = factor(rep(names(merge_cluster_ids), unlist(lapply(merge_cluster_ids, length))/10), levels = names(merge_cluster_ids)),
                    probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), column_title = NA,
                    colour_scheme = scDNA_CN_colour)
  dev.off()
  
}

