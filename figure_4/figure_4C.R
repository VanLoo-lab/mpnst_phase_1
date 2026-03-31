### title: UMAP of 10X scDNA total CN profiles and projection of Visium total CN profiles

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
library(pheatmap)
library(pbmcapply)
library(dendextend)
library(umap)

### Load scDNA data
gamma = 5
scDNA_CN_mtx    <- readRDS(paste0(input.dir,"scDNA/CN_profiles/MPNST_all_CN_mtx_mpcf_", gamma,".rds"))
all_probes      <- readRDS(paste0(input.dir,"scDNA/CN_profiles/MPNST_all_probes_mpcf_", gamma,".rds"))
all_chr_probes  <- readRDS(paste0(input.dir,"scDNA/CN_profiles/MPNST_all_chr_probes_mpcf_", gamma,".rds"))

### Read in the cluster IDs for the scDNA-seq data generated after running kmeans
K=22
scDNA_CN_cluster_ids <- readRDS(paste0(input.dir, "scDNA/CN_profiles/MPNST_all_k_means_K",K,"_clusters.rds"))

### Load CN profiles from inferCNV
visium.infercnv.dir <- paste0(input.dir, "spRNA/infercnv/")
# loaded into object segmented_mtx by sample in the for loop below

### Load stored UMAP
scDNA_umap <- readRDS(paste0(input.dir, "scDNA/UMAP/10X_DLP_totCN_umap.rds"))

### Load spRNA CN reference
R1_ref_CN <- readRDS(paste0(input.dir, "spRNA/MPNST_R1_ref_CN.rds"))

####################################################################################################################################
### Part 1: Prepare metadata
####################################################################################################################################


### Define sample names
visium_samples = c("Pa", "R1", "R3", "Pb", "R2", "R4", "R5a", "R5b")
names(visium_samples) = c("VER683A1", "VER683A2", "VER683A3", "VER683A4", "VER683A5", "VER683A6", "VER683A7", "VER683A8")
# sample = "C"
color_scale <- colorRampPalette(c("darkblue", "white", "darkred"))(100)

combos <- combn(visium_samples[c(2:8,1)], 2)
pairs <- paste0(combos[1,], "_", combos[2,])


### Define sample name and colors
samples = c("R1", "R2", "R3", "R4", "R5", "P")
names(samples) = c("FIT208A3", "FIT208A4", "FIT208A5", "FIT208A6", "FIT208A7", "FIT208A8")
# SUBSETS = c("P", "R1", "R2", "R3", "R4", "R5", "G2MS", "Ribo")

#Colours
region_colours <- c("#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3", "#F8766D")[c(6,1:5)]
names(region_colours) <- samples[c(6,1:5)]

### Define clustering parameters
rm_small_cluster <- 0 #Set to 50 or 0
adj_cluster <- F
high_res <- F
suffix <- paste0(ifelse(rm_small_cluster > 0, paste0("_rm_sub", rm_small_cluster), ""), ifelse(adj_cluster, "_adj", ""), ifelse(high_res, "_res", ""))

### Function to discretize
disc_mtx <- function(mtx, breaks, labels) {
  d_mtx <- matrix(as.matrix(cut(mtx, breaks = breaks, labels = labels)) %>% as.numeric(), nrow(mtx))
  rownames(d_mtx) <- rownames(mtx)
  colnames(d_mtx) <- colnames(mtx)
  return(d_mtx)
}

discretize_breaks <- c(-Inf,-0.2,-0.05,0.05,0.2,Inf)
discretize_labels <- c(-2,-1,0,1,2)

### Convert Visium CN profiles to total CN profiles and prepare for projection onto scDNA UMAP
visium_totCN_profiles <- do.call(rbind, lapply(visium_samples[-2], function(slide) {
  ##Load discritized matrix
  segmented_mtx <- readRDS(paste0(visium.infercnv.dir, "MPNST_R1_",slide,"_cluster_segmented_mtx",suffix,".rds"))
  
  #Convert segmented mtx to medicc input 
  medicc_mtx <- disc_mtx(segmented_mtx, breaks = discretize_breaks, labels = discretize_labels)
  
  #convert CN profiles
  medicc_totCN_mtx <- do.call(rbind, lapply(1:nrow(medicc_mtx), function(c) {
    return(R1_ref_CN %>% mutate(cn_a = ifelse(is.na(totCN + medicc_mtx[c,]), totCN, totCN + medicc_mtx[c,])) %>% pull(cn_a))
  }))
  rownames(medicc_totCN_mtx) <- rownames(medicc_mtx)
  
  medicc_totCN_plot_mtx <- do.call(rbind, lapply(1:nrow(medicc_mtx), function(c) {
    CN <- R1_ref_CN %>% mutate(cn_a = ifelse(is.na(totCN + medicc_mtx[c,]), totCN, totCN + medicc_mtx[c,])) %>% pull(cn_a)
    return(rep(CN, as.list(R1_ref_CN$n.probes)))
  }))
  rownames(medicc_totCN_plot_mtx) <- rownames(segmented_mtx)
  return(medicc_totCN_plot_mtx)
}))

### scDNA metadata
K = 22
run <- "22_kmeans_10X_DLP"
remove_kmeans2 = T
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

cluster_ids <- data.frame(barcode = c(unlist(scDNA_CN_cluster_ids), scDNA_normal_ids),
                          Cluster = factor(c(unlist(lapply(1:length(scDNA_CN_cluster_ids), function(c) return(rep(c, length(scDNA_CN_cluster_ids[[c]]))))),
                                             rep("Diploid", length(scDNA_normal_ids))), levels = c(1:length(scDNA_CN_cluster_ids), "Diploid")),
                          Data = gsub(".*_(10X|DLP)_.*", "\\1", c(unlist(scDNA_CN_cluster_ids), scDNA_normal_ids)))

names(scDNA_CN_cluster_ids) <- ave(scDNA_CN_clusters_origins, scDNA_CN_clusters_origins, FUN = function(i) paste0(i, "_", seq_along(i)))

totCN_mtx_clean <- scDNA_CN_mtx[c(unlist(scDNA_CN_cluster_ids),scDNA_normal_ids),]



####################################################################################################################################
### Make the UMAPs
####################################################################################################################################

### Define UMAP axes through scDNA_umap
umap_plot <- scDNA_umap[["layout"]] %>% as.data.frame() %>% 
  rename("UMAP1" = "V1", "UMAP2" = "V2") %>% 
  rownames_to_column("barcode") %>%
  mutate(Region = gsub("_.*", "", rownames(totCN_mtx_clean)), Section = " ") %>% 
  left_join(cluster_ids, by = "barcode")

### Project onto spRNA UMAP
scDNA_10X_visium_umap <- predict(scDNA_umap, visium_totCN_profiles)


scDNA_10X_visium_umap_plot <- rbind(umap_plot, 
                                    data.frame(barcode = rownames(scDNA_10X_visium_umap),
                                               UMAP1 = scDNA_10X_visium_umap[,1],
                                               UMAP2 = scDNA_10X_visium_umap[,2],
                                               Region = gsub("a|b|_.*", "", rownames(visium_totCN_profiles)),
                                               Section = gsub("_.*", "", rownames(visium_totCN_profiles)),
                                               Cluster = rownames(visium_totCN_profiles),
                                               Data = "Visium") %>% remove_rownames())



### UMAP with colors by region
pdf(file = paste0("figure_4C_region.pdf"), width = 7, height = 7)
print(
  scDNA_10X_visium_umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = Region, size = Data, alpha = Data)) + 
    geom_point() + 
    scale_size_manual(values = c(0.5, 0.5, 8)) + 
    scale_alpha_manual(values = c(1, 1, 0.2)) + 
    theme_classic(base_size = 24) + 
    theme(legend.position = c(0.77, 0.3), legend.box = "horizontal",
          axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
    guides(color = guide_legend(override.aes = list(size=6)))
)
dev.off()

### UMAP with colors by section
scDNA_10X_visium_umap_plot <- scDNA_10X_visium_umap_plot %>% mutate(Section = ifelse(Data == "Visium", Section, ifelse(Region %in% c("P", "R5"), paste0(Region,"a"), Region)))

pdf(file = paste0("figure_4C_section.pdf"), width = 7, height = 7)
print(
  scDNA_10X_visium_umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = Section, size = Data, alpha = Data)) + 
    geom_point() + 
    scale_color_manual(values = c(
      "R1"  = unname(region_colours["R1"]),
      "R2"  = unname(region_colours["R2"]),
      "R3"  = unname(region_colours["R3"]),
      "R4"  = unname(region_colours["R4"]),
      "R5a" = unname(region_colours["R5"]),
      "R5b" = unname("darkorchid4"),
      "Pa"  = unname(region_colours["P"]),
      "Pb"  = unname("darkgoldenrod1")
    )) + 
    scale_size_manual(values = c(0.5, 0.5, 8)) + 
    scale_alpha_manual(values = c(1, 1, 0.2)) + 
    theme_classic(base_size = 24) + 
    theme(
      legend.position = c(0.77, 0.3),
      legend.box = "horizontal",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) + 
    guides(color = guide_legend(override.aes = list(size = 6)))
)
dev.off()

