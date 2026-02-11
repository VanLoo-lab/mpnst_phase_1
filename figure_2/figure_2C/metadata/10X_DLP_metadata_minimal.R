####################################################################################################################################
### Some directory names
####################################################################################################################################

medicc.dir = paste0(INPUTDIR, "10x_DLP/MEDICC2/")

samples = c("R1", "R2", "R3", "R4", "R5", "P")
names(samples) = c("FIT208A3", "FIT208A4", "FIT208A5", "FIT208A6", "FIT208A7", "FIT208A8")
# SUBSETS = c("P", "R1", "R2", "R3", "R4", "R5", "G2MS", "Ribo")

bulk_samples = c("VER236A1", "VER236A2", "VER236A3", "VER236A4", "VER236A5", "VER236A6")
names(bulk_samples) <- samples

#Colours
region_colours <- c("#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3", "#F8766D")[c(6,1:5)]
names(region_colours) <- samples[c(6,1:5)]

####################################################################################################################################
### Cluster info generated after running kmeans in ASCAT.sc.R then cleaned up
####################################################################################################################################

INPUTDIR
tenX_DLP_ASCAT.dir <- paste0(INPUTDIR,"10x_DLP/ASCAT.sc/")

K = 22
run <- "22_kmeans_10X_DLP"
remove_kmeans2 = T
scDNA_CN_cluster_ids <- readRDS(paste0(tenX_DLP_ASCAT.dir, "MPNST_all_k_means_K",K,"_clusters.rds"))
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


