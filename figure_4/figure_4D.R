### title: Subclones and phylogenies in space

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
library(ape)
library(ggtree)
library(Seurat)

### Load scDNA data
all_probes <- readRDS(paste0(input.dir, "scDNA/CN_profiles/MPNST_all_probes_mpcf_5.rds"))
all_chr_probes <- readRDS(paste0(input.dir, "scDNA/CN_profiles/MPNST_all_chr_probes_mpcf_5.rds"))

### Load spRNA data
MPNST_sp_markers <- readRDS(paste0(input.dir, "spRNA/MPNST_sp_markers.rds"))

### Load spRNA reference
R1_ref_CN <- readRDS(paste0(input.dir, "spRNA/MPNST_R1_ref_CN.rds"))

### Load inferCNV results
inferCNV.dir <- paste0(input.dir,"spRNA/inferCNV/")
# loaded per sample in for loop below

### Load MEDICC2 trees from LCM spots
medicc2.dir <- paste0(input.dir, "LCM/MEDICC2/")

####################################################################################################################################
### Part 1: Prepare metadata
####################################################################################################################################

### Define color pallette
color_scale <- colorRampPalette(c("darkblue", "grey95", "darkred"))(100)
scRNA_mtx_breaks = seq(0.8,1.2, by = 0.004)
scRNA_adj_mtx_breaks = seq(-0.2,0.2, by = 0.004)

# Variables for plotting CN heatmap
scDNA_CN_colour <- c("#00008B", "#7B7BC3", "#FFFFFF", "#FFCCCC", "#FF8080", "#FF3333", "#E60000", "#B31000", "#8B0000", "#660000", "#330000")
scDNA_CN_breaks = seq(-0.5,10.5, by = 1)

### Function to discretize
disc_mtx <- function(mtx, breaks, labels) {
  d_mtx <- matrix(as.matrix(cut(mtx, breaks = breaks, labels = labels)) %>% as.numeric(), nrow(mtx))
  rownames(d_mtx) <- rownames(mtx)
  colnames(d_mtx) <- colnames(mtx)
  return(d_mtx)
}

discretize_breaks <- c(-Inf,-0.2,-0.05,0.05,0.2,Inf)
discretize_labels <- c(-2,-1,0,1,2)
discretize_colors = c("darkblue", "royalblue1", "grey95", "firebrick3", "darkred")

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

####################################################################################################################################
### Part 2: Plot phylogeny and subclones on image
####################################################################################################################################


for (sample in pairs[1:7]) {
  ### Load spRNA data
  segmented_mtx <- readRDS(paste0(inferCNV.dir,"MPNST_",sample,"_cluster_segmented_mtx",suffix,".rds"))
  phylos <- readRDS(paste0(inferCNV.dir,"MPNST_",sample,"_cluster_trees",suffix,".rds"))
  merge_cluster_ids <- readRDS(paste0(inferCNV.dir,"MPNST_", sample, "_clusters",suffix,".rds"))
  
  ### Plot medicc tree on spots
  if (T & !single_spot) {
    medicc_input_name <- paste0("MPNST_visium_clusters_", gsub(".*_", "", sample),suffix)
    medicc_visium_tree <- read.tree(paste0(medicc2.dir, medicc_input_name, "_final_tree.new"))
    
    #Calculate centroid of clusters
    if (!kd_centroid) {
      #Work out average location of clone
      if (!file.exists(paste0("MPNST_sp_",gsub(".*_", "", sample),"_subclone_locations",suffix,".rds"))) {
        max_imagerow = max(GetTissueCoordinates(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]])$imagerow)
        min_imagerow = min(GetTissueCoordinates(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]])$imagerow)
        subclone_locations <- do.call(rbind, lapply(merge_cluster_ids, function (c) {
          locations <- GetTissueCoordinates(subset(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]], cells = c)) %>% rownames_to_column(var = "spot")
          return(data.frame(imagecol = mean(locations$imagecol), imagerow = max_imagerow+min_imagerow-mean(locations$imagerow))) 
        })) %>% add_column(cluster = as.character(1:length(merge_cluster_ids)), cluster_name = names(merge_cluster_ids), .before = "imagecol")
      } else {
        subclone_locations <- readRDS(paste0("MPNST_sp_",gsub(".*_", "", sample),"_subclone_locations",suffix,".rds"))
      }
    } else {
      #Attempt to use kernel density for subclone location
      if (!file.exists(paste0("MPNST_sp_",gsub(".*_", "", sample),"_subclone_kdlocations",suffix,".rds"))) {
        max_imagerow = max(GetTissueCoordinates(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]])$imagerow)
        min_imagerow = min(GetTissueCoordinates(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]])$imagerow)
        
        locations_all <- GetTissueCoordinates(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]]) %>% rownames_to_column(var = "spot") %>%
          mutate(imagerow = max_imagerow+min_imagerow-imagerow) %>%
          left_join(unnest(enframe(merge_cluster_ids)), by = c("spot" = "value")) %>% dplyr::rename("Cluster" = "name")
        
        subclone_locations <- do.call(rbind, lapply(names(merge_cluster_ids), function (c) {
          locations <- locations_all %>% filter(Cluster == c)
          de2d <- MASS::kde2d(x = locations$imagerow, y = locations$imagecol, n = 100)
          max_density_grid <- which(de2d$z == max(de2d$z), arr.ind = T)
          density_row <- de2d$x[max_density_grid[1,1]]
          density_col <- de2d$y[max_density_grid[1,2]]
          return(data.frame(cluster_name = c, imagecol = density_col, imagerow = density_row))
        })) %>% add_column(cluster = as.character(1:length(merge_cluster_ids)), .before = "cluster_name") %>% mutate(cluster_name = as.character(cluster_name))
        rownames(subclone_locations) <- subclone_locations$cluster_name
      } else {
        subclone_locations <- readRDS(paste0("MPNST_sp_",gsub(".*_", "", sample),"_subclone_kdlocations",suffix,".rds"))
      }
    }
    
    #Add hclust subclone to metadata
    MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]]$subclone <- gsub("Spatial_", "", MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]]@active.ident)
    MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]]$subclone <- unlist(lapply(1:length(merge_cluster_ids), function(k) {
      clusters <- rep(names(merge_cluster_ids)[k], length(merge_cluster_ids[[k]]))
      names(clusters) <- merge_cluster_ids[[k]]
      return(clusters)
    }))
    
    #Generate location of intermediate nodes
    all <- c(medicc_visium_tree[["tip.label"]], "root", medicc_visium_tree[["node.label"]][-1])
    nodes <- medicc_visium_tree[["node.label"]][-1]
    edges <- data.frame(parent = all[medicc_visium_tree[["edge"]][,1]],
                        child = all[medicc_visium_tree[["edge"]][,2]],
                        length = medicc_visium_tree[["edge.length"]]) %>% left_join(subclone_locations %>% select(cluster_name, imagecol, imagerow), by = c("child" = "cluster_name")) %>%
      filter(str_detect(child, "diploid", negate = T)) %>% mutate(edge = 1:nrow(.)) %>% rename(coord.x = imagecol, coord.y = imagerow)
    
    #Iterate through to add locations
    while(any(is.na(edges))) {
      for(i in 1:nrow(print(edges))) {
        node = edges[i,"parent"]
        children = edges %>% filter(parent == node) %>% pull(child)
        edges[edges$child == node,"coord.x"] <- ifelse(any(is.na(edges %>% filter(child %in% children))), NA, 
                                                       (edges %>% filter(child == children[1]) %>% pull(coord.x) + edges %>% filter(child == children[2]) %>% pull(coord.x))/2)
        edges[edges$child == node,"coord.y"] <- ifelse(any(is.na(edges %>% filter(child %in% children))), NA, 
                                                       (edges %>% filter(child == children[1]) %>% pull(coord.y) + edges %>% filter(child == children[2]) %>% pull(coord.y))/2)
      }
      print(edges)
    }
    
    #Adjust locations
    pct_adj = 0.3
    edges_adj <- edges
    for (i in nrow(edges):1) {
      current_node <- edges_adj[i,"child"]
      if(grepl("internal", current_node)) {
        print(current_node)
        child_x <- edges_adj[edges_adj[,"child"] == current_node, "coord.x"]
        child_y <- edges_adj[edges_adj[,"child"] == current_node, "coord.y"]
        #Get parent node location
        parent_node <- edges_adj[edges_adj[,"child"] == current_node, "parent"]
        if (parent_node != "root") {
          parent_x <- edges_adj[edges_adj[,"child"] == parent_node, "coord.x"]
          parent_y <- edges_adj[edges_adj[,"child"] == parent_node, "coord.y"]
          #Adjust and update current node
          adj_x <- child_x + (parent_x-child_x)*pct_adj
          adj_y <- child_y + (parent_y-child_y)*pct_adj
          edges_adj[edges_adj[,"child"] == current_node, "coord.x"] <- adj_x
          edges_adj[edges_adj[,"child"] == current_node, "coord.y"] <- adj_y
        }
      }
    }
    
    #Rename to highest node (first row) to MRCA
    edges_adj <- edges_adj %>% mutate(parent = str_replace(parent, paste0(edges_adj[1,"child"],"$"), "MRCA"))
    edges_adj[1,"child"] <- "MRCA"
    
    edges_ggplot <- rbind(edges_adj[,c(2,4,5,6)] %>% mutate(direction = "end"), 
                          edges_adj[,c(1,6)] %>% left_join(edges_adj[,c(2,4,5)], by = c("parent" = "child")) %>% select(1,3,4,2) %>% rename(child = parent) %>% mutate(direction = "start")) %>% 
      mutate(type = ifelse(str_detect(child, gsub(".*_", "", sample)), child, ifelse(str_detect(child, "MRCA"), "MRCA", "Node"))) %>%filter(str_detect(child, "root", negate = T))
    spot_colors <- c("navy", NA, hcl(h = seq(15, 375, length = length(unique(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]]$subclone)) + 1), l = 65, c = 100) %>% head(-1))
    names(spot_colors) <- c("MRCA", "Node", sort(unique(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]]$subclone)))

    MPNST_sp_markers <- lapply(MPNST_sp_markers, UpdateSeuratObject)
    
    pdf(file = paste0("figure_4D_", sample,".pdf"), width = 7, height = 7)
    print(SpatialPlot(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]], group.by = "subclone", pt.size.factor = 45, stroke = NA) + labs(color = "Subclone") + scale_fill_manual(values = spot_colors[-c(1:2)], guide = "none") +
            geom_point(data = edges_ggplot, mapping = aes(x = coord.x, y = coord.y, color = type), shape = 16, size = 8, inherit.aes = F) + scale_color_manual(values = spot_colors) +
            geom_path(data = edges_ggplot %>% arrange(edge, desc(direction)), mapping = aes(x = coord.x, y = coord.y, group = edge), arrow = arrow(length=unit(12,"pt")), colour = "navy", size = 1, inherit.aes = F) +
            theme(legend.position = "right"))
    dev.off()
  }
}


