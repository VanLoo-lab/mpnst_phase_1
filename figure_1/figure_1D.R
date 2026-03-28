### title: Analysis and plotting of scDNA UMAPS

# Comment: Part 1 creates the UMAP but can be skipped by directly loading
# the scDNA_umap object from the input repository if desired. Alternatively,
# create scDNA_umap object from the CN profile.

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
#library(umap)

### Read in the CN profile for the scDNA-seq data (with gamma = 5 used for clustering)
gamma = 5
scDNA_CN_mtx    <- readRDS(paste0(input.dir,"CN_profiles/MPNST_all_CN_mtx_mpcf_", gamma,".rds"))

### Read in the UMAP coordinates for the scDNA-seq data
scDNA_umap <- readRDS(paste0(input.dir, "UMAP/10X_DLP_totCN_umap.rds"))

### Read in the cluster IDs for the scDNA-seq data generated after running kmeans
K=22
scDNA_CN_cluster_ids <- readRDS(paste0(input.dir, "CN_profiles/MPNST_all_k_means_K",K,"_clusters.rds"))

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

totCN_mtx_clean <- scDNA_CN_mtx[c(unlist(scDNA_CN_cluster_ids),scDNA_normal_ids),]


####################################################################################################################################
### Part 2: Plot UMAPs for scDNA-seq data colored by region and data type
####################################################################################################################################

### Define a data frame with cluster and data type information for each barcode
cluster_ids <- data.frame(barcode = c(unlist(scDNA_CN_cluster_ids), scDNA_normal_ids),
                          Cluster = factor(c(unlist(lapply(1:length(scDNA_CN_cluster_ids), function(c) return(rep(names(scDNA_CN_cluster_ids)[c], length(scDNA_CN_cluster_ids[[c]]))))),
                                             rep("Diploid", length(scDNA_normal_ids))), levels = c(names(scDNA_CN_cluster_ids), "Diploid")),
                          Data = gsub(".*_(10X|DLP)_.*", "\\1", c(unlist(scDNA_CN_cluster_ids), scDNA_normal_ids)))

### Create a data frame for plotting the UMAP coordinates with ggplot2
umap_plot <- scDNA_umap[["layout"]] %>% as.data.frame() %>% 
  `colnames<-`(c("UMAP1", "UMAP2")) %>% 
  #rename("UMAP1" = "V1", "UMAP2" = "V2") %>% 
  rownames_to_column("barcode") %>%
  mutate(Region = gsub("_.*", "", rownames(totCN_mtx_clean))) %>% 
  left_join(cluster_ids, by = "barcode")

### Plot UMAP colored by region
pdf(file = paste0("figure_1D1.pdf"), width = 7, height = 7)
print(
  umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = Region)) + 
    geom_point(shape = 16, size = 1) + 
    theme_classic(base_size = 24) + 
    theme(legend.position = c(0.9, 0.3),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    guides(color = guide_legend(override.aes = list(size=6)))
)
dev.off()


### Plot UMAP colored by data type
pdf(file = paste0("figure_1D2.pdf"), width = 7, height = 7)
print(
  umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = Data)) + 
    geom_point(shape = 16, size = 1) + 
    theme_classic(base_size = 24) + 
    theme(legend.position = c(0.9, 0.3),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + 
    scale_color_manual(values = c("mediumpurple1", "olivedrab3")) + guides(color = guide_legend(override.aes = list(size=6)))
)
dev.off()

