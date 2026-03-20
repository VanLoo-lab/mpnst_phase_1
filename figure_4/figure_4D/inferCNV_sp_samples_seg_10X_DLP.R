
####################################################################################################################################
### Part 0: Code preperations
####################################################################################################################################

# options(bitmapType='cairo') #to solve plotting issue on CAMP

library(tidyverse)
library(pheatmap)
library(dendextend)
library(pbmcapply)
library(grid)
library(gridExtra)
library(ape)
library(phangorn)
library(ggtree)
library(Seurat)
library(ComplexHeatmap)

### Load metadata
#source(file = "/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/code/10X_DLP_metadata.R")
#source(file = "/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DNA/code/MPNST_common_functions.R")
#source(file = "/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_spatial/code/Visium_metadata.R")

### Define input and output directories
INPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_1/data/"
CODEDIR <- "~/Documents/GitHub/mpnst_phase_1/figure_4/figure_4D/"
OUTPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_4/results/figure_4D/"

### Load metadata
source( paste0(CODEDIR, "metadata/MPNST_common_functions.R") )
source( paste0(CODEDIR, "metadata/10X_DLP_metadata_minimal.R") )
source( paste0(CODEDIR, "metadata/Visium_metadata_minimal.R") )

### Define directories
output.dir = OUTPUTDIR
setwd(output.dir)

####################################################################################################################################
### Part 1: Load data 
####################################################################################################################################

#Load segments from 10X DNA
all_10X_DLP.dir = paste0(INPUTDIR, "10x_DLP/ASCAT.sc/")
# pcf.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DNA/results/pcf/")
# if (F) {
#   load(paste0(pcf.dir, "mpcf/multi_pcf_all_gammas_filtered_cells.Rda"))
#   totCN_chr_probes <- do.call(rbind, lapply(allsegs_pf[["gamma_17"]], function(c) {
#     return(c[,c(1, 3:5)])
#   })) %>% mutate(chrom = paste0("chr", ifelse(chrom == 23, "X", chrom))) %>% rename("chr" = "chrom", "startpos" = "start.pos", "endpos" = "end.pos")
#   saveRDS(totCN_chr_probes, paste0(pcf.dir, "MPNST_scDNA_totCN_probes.rds"))
#   
# }
all_probes <- readRDS(paste0(all_10X_DLP.dir, "MPNST_all_probes_mpcf_5.rds"))
all_chr_probes <- readRDS(paste0(all_10X_DLP.dir, "MPNST_all_chr_probes_mpcf_5.rds"))

color_scale <- colorRampPalette(c("darkblue", "grey95", "darkred"))(100)
scRNA_mtx_breaks = seq(0.8,1.2, by = 0.004)
scRNA_adj_mtx_breaks = seq(-0.2,0.2, by = 0.004)

#Load in variables for plotting CN heatmap
scDNA_CN_colour <- c("#00008B", "#7B7BC3", "#FFFFFF", "#FFCCCC", "#FF8080", "#FF3333", "#E60000", "#B31000", "#8B0000", "#660000", "#330000")
scDNA_CN_breaks = seq(-0.5,10.5, by = 1)

#Function to discretize
disc_mtx <- function(mtx, breaks, labels) {
  d_mtx <- matrix(as.matrix(cut(mtx, breaks = breaks, labels = labels)) %>% as.numeric(), nrow(mtx))
  rownames(d_mtx) <- rownames(mtx)
  colnames(d_mtx) <- colnames(mtx)
  return(d_mtx)
}

#Load spatial
input.dir <- paste0(INPUTDIR, "10X_spatial/")
if (T) {
  MPNST_sp_markers <- readRDS(paste0(input.dir, "MPNST_sp_markers.rds"))
  for (s in visium_samples) {
    MPNST_sp_markers[[which(visium_samples == s)]] <- RenameCells(MPNST_sp_markers[[which(visium_samples == s)]],
                                                                  new.names = gsub(paste0("(................)-1"), paste0(s, "_\\1","SP"),colnames(MPNST_sp_markers[[which(visium_samples == s)]])))
    sp.cluster.ids <- paste0(s, "_",0:(length(levels(MPNST_sp_markers[[which(visium_samples == s)]]))-1))
    names(sp.cluster.ids) <- levels(MPNST_sp_markers[[which(visium_samples == s)]])
    MPNST_sp_markers[[which(visium_samples == s)]] <- RenameIdents(MPNST_sp_markers[[which(visium_samples == s)]], sp.cluster.ids)
  }
}

####################################################################################################################################
### Part 2: 2nd attempt (per cluster)
####################################################################################################################################
#Generate ref R1 CN profile
if (!file.exists("MPNST_R1_ref_CN.rds")) {
  #Load R1 subclone and create median CN
  # scDNA_CN_mtx <- readRDS(paste0(pcf.dir, "fittedCN/MPNST_scDNA_CN_mtx.rds"))
  scDNA_CN_mtx <- readRDS(paste0(all_10X_DLP.dir, "MPNST_all_CN_mtx_mpcf_5.rds"))
  
  #Plot R1_1
  # png(filename = paste0("MPNST_scDNA_CN_K",run,"_R1_2.png"), width = 4000, height = 1000, res = 200)
  # 
  # ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[scDNA_CN_cluster_ids[["R1_2"]],])),
  #                                        Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx[scDNA_CN_cluster_ids[["R1_2"]],]))),
  #                        col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
  #                                   Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")), show_annotation_name = F)
  # print(sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[scDNA_CN_cluster_ids[["R1_2"]],], hclust = F, row_ann = ha_row, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), title = "MPNST All \nTotal CN "))
  # dev.off()
  
  R1_2_split <- cutree(hclust(as.dist(1-cor(t(scDNA_CN_mtx[scDNA_CN_cluster_ids[["R1_2"]],]))), method = "ward.D2"), k=2, order_clusters_as_data = F)
  R1_2_sub <- names(R1_2_split[R1_2_split == 1])
  
  #Plot R1 subclone
  # png(filename = paste0("MPNST_scDNA_CN_R1_subclone.png"), width = 4000, height = 1000, res = 200)
  # sc_totCN_heatmap(scDNA_CN_mtx[R1_2_sub,], hclust = F, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")),
  #                  colour_scheme = scDNA_CN_colour)
  # dev.off()
  
  #Consensus CN mtx
  R1_mode_mtx <- matrix(apply(scDNA_CN_mtx[R1_2_sub,], 2, function(m) {
    names(sort(table(m), decreasing = T))[1] %>% as.numeric()
  }), nrow = 1)
  
  R1_ref_CN <- all_probes %>% mutate(probes.start = lag(cumsum(n.probes)+1, default = 1), probes.end = cumsum(n.probes))
  R1_ref_CN <- R1_ref_CN %>% mutate(totCN = unlist(lapply(1:nrow(R1_ref_CN), function(s) {
    return(R1_mode_mtx[,R1_ref_CN[s,"probes.start"]:R1_ref_CN[s,"probes.end"]] %>% table %>% sort(decreasing = T) %>% names %>% .[1] %>% as.numeric())
  })))
  saveRDS(R1_ref_CN, "MPNST_R1_ref_CN.rds")
  
  #Plot
  # R1_ref_mtx <- matrix(rep(R1_ref_CN$totCN, as.list(R1_ref_CN$n.probes)), nrow = 1)
  # png(filename = paste0("MPNST_scDNA_CN_K",run,"_R1_2_sub_mode.png"), width = 4000, height = 1000, res = 200)
  # sc_totCN_heatmap(R1_ref_mtx, hclust = F, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")),
  #                  colour_scheme = scDNA_CN_colour)
  # print(pheatmap(mat = R1_ref_mtx, cluster_rows = F, cluster_cols = F,
  #                show_rownames = F, show_colnames = F, color = scDNA_CN_colour, breaks = scDNA_CN_breaks, fontsize = 14, main = "MPNST scDNA Fitted CN Heatmap"))
  # grid.lines(x=chr_probes$cum.probes[1]/ncol(scDNA_CN_mtx)*0.961+0.004, y=c(0.004,0.984), gp=gpar(col="black", lwd=2))
  # for (c in chr_probes$chrom[-1]) {
  #   grid.lines(x=chr_probes$cum.probes[c+1]/ncol(scDNA_CN_mtx)*0.961+0.004, y=c(0.004,0.984), gp=gpar(col="black", lwd=2))
  #   grid.text(chr_probes$chrom[c+1], x=chr_probes$cum.probes[c]/ncol(scDNA_CN_mtx)*0.961+0.014, y=0.01, gp=gpar(cex = 2))
  # }
  # dev.off()
} else {
  R1_ref_CN <- readRDS("MPNST_R1_ref_CN.rds")
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


####################################################################################################################################
### Part 3: Figure 4D
####################################################################################################################################


# #Run medicc (using infered total CN) only for clusters (not single spots)
# if (F & !single_spot) {
#   no_wgd <- F 
#   visium_clusters_input <- do.call(rbind, lapply(1:nrow(medicc_totCN_mtx), function(s) {
#     all_probes %>% add_column(sample_id = rownames(medicc_totCN_mtx)[s], .before = "chrom") %>% mutate(cn_a = replace_na(medicc_totCN_mtx[s,], 2), cn_b = 0) %>% select(-n.probes)
#   }))
#   
#   #Make medicc input
#   medicc_input_name <- paste0("MPNST_visium_clusters_", gsub(".*_", "", sample),suffix)
#   options(scipen=10)
#   write.table(visium_clusters_input, file = paste0(medicc_input_name, ".tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
#   options(scipen=0)
#   
#   #Run medicc manually
#   min_seg = 1.5
#   cat("ml purge\nml Anaconda3\nconda activate medicc_env_2108\n")
#   cat(paste0("python /camp/home/yanh/lab_share_working/yanh/medicc2_2108/medicc2/medicc2.py --input-type t ", ifelse(no_wgd, "--no-wgd ", ""),
#              "--filter-segment-length ", min_seg*1000000, " --no-plot ",
#              output.dir, medicc_input_name, ".tsv ",
#              output.dir, medicc_input_name, "\n"))
# }

# invisible(readline(prompt="Press [enter] after medicc run"))

input.dir2 <- paste0(INPUTDIR,"10X_spatial/inferCNV_samples_seg_10X_DLP/")

for (sample in pairs[1:7]) {
  
  segmented_mtx <- readRDS(paste0(input.dir2,"MPNST_",sample,"_cluster_segmented_mtx",suffix,".rds"))
  phylos <- readRDS(paste0(input.dir2,"MPNST_",sample,"_cluster_trees",suffix,".rds"))
  merge_cluster_ids <- readRDS(paste0(input.dir2,"MPNST_", sample, "_clusters",suffix,".rds"))
  
  ### Plot medicc tree on spots
  
  if (T & !single_spot) {
    medicc_input_name <- paste0("MPNST_visium_clusters_", gsub(".*_", "", sample),suffix)
    medicc_visium_tree <- read.tree(paste0(output.dir, "final_tree/", medicc_input_name, "_final_tree.new"))
    
    #Calculate centroid of clusters
    if (!kd_centroid) {
      #Work out average location of clone
      if (!file.exists(paste0("MPNST_sp_",gsub(".*_", "", sample),"_subclone_locations",suffix,".rds"))) {
        #Get locations
        # test <- lapply(merge_cluster_ids, function (c) {
        #      locations <- GetTissueCoordinates(subset(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]], cells = c)) %>% rownames_to_column(var = "spot")
        #      return(data.frame(imagecol = (locations$imagecol), imagerow = (locations$imagerow))) 
        #  })
        # max_imagerow = max(GetTissueCoordinates(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]])$imagerow)*1.25 #Silly Seurat bug reverse the y axis means have to flip axis manually with largest imagerow point
        max_imagerow = max(GetTissueCoordinates(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]])$imagerow)
        min_imagerow = min(GetTissueCoordinates(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]])$imagerow)
        subclone_locations <- do.call(rbind, lapply(merge_cluster_ids, function (c) {
          locations <- GetTissueCoordinates(subset(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]], cells = c)) %>% rownames_to_column(var = "spot")
          return(data.frame(imagecol = mean(locations$imagecol), imagerow = max_imagerow+min_imagerow-mean(locations$imagerow))) 
        })) %>% add_column(cluster = as.character(1:length(merge_cluster_ids)), cluster_name = names(merge_cluster_ids), .before = "imagecol")
        saveRDS(subclone_locations, paste0("MPNST_sp_",gsub(".*_", "", sample),"_subclone_locations",suffix,".rds"))
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
          left_join(unnest(enframe(merge_cluster_ids)), by = c("spot" = "value")) %>% rename("Cluster" = "name")
        
        # ggplot(locations_all, aes(x = imagerow, y = imagecol, colour = Cluster)) + geom_point()
        # ggplot(locations_all, aes(x = imagerow, y = imagecol)) + stat_density_2d(aes(alpha = ..level.., fill = Cluster), geom = "polygon", colour="white", bins = 100)
        
        subclone_locations <- do.call(rbind, lapply(names(merge_cluster_ids), function (c) {
          locations <- locations_all %>% filter(Cluster == c)
          de2d <- MASS::kde2d(x = locations$imagerow, y = locations$imagecol, n = 100)
          max_density_grid <- which(de2d$z == max(de2d$z), arr.ind = T)
          density_row <- de2d$x[max_density_grid[1,1]]
          density_col <- de2d$y[max_density_grid[1,2]]
          return(data.frame(cluster_name = c, imagecol = density_col, imagerow = density_row))
        })) %>% add_column(cluster = as.character(1:length(merge_cluster_ids)), .before = "cluster_name") %>% mutate(cluster_name = as.character(cluster_name))
        rownames(subclone_locations) <- subclone_locations$cluster_name
        # ggplot(locations_all, aes(x = imagerow, y = imagecol)) + geom_point(aes(colour = Cluster)) +
        #   geom_point(data = subclone_locations, aes(colour = cluster_name, size = 5))
        saveRDS(subclone_locations, paste0("MPNST_sp_",gsub(".*_", "", sample),"_subclone_kdlocations",suffix,".rds"))
      } else {
        subclone_locations <- readRDS(paste0("MPNST_sp_",gsub(".*_", "", sample),"_subclone_kdlocations",suffix,".rds"))
      }
    }
    
    #Add hclust subclone to metadata
    # infercnv_subclone <- cutree(MPNST_hclust_cor_ward, k=K, order_clusters_as_data = F)
    MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]]$subclone <- gsub("Spatial_", "", MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]]@active.ident)
    MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]]$subclone <- unlist(lapply(1:length(merge_cluster_ids), function(k) {
      clusters <- rep(names(merge_cluster_ids)[k], length(merge_cluster_ids[[k]]))
      names(clusters) <- merge_cluster_ids[[k]]
      return(clusters)
    }))
    #Plot
    pdf(file = paste0("MPNST_", sample, "_spatialplot_subclones",suffix,".pdf"), width = 7, height = 7)
    print(SpatialPlot(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]], group.by = "subclone", pt.size.factor = 1.2, stroke = NA) +
            geom_point(data = subclone_locations, mapping = aes(x = imagecol, y = imagerow, color = cluster_name), shape = 16, size = 8, inherit.aes = F, show.legend = F) +
            theme(legend.position = "right"))
    dev.off()
    
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
    
    pdf(file = paste0("MPNST_", sample, "_spatialplot_subclones_tree",suffix,".pdf"), width = 7, height = 7)
    print(SpatialPlot(MPNST_sp_markers[[which(visium_samples == gsub(".*_", "", sample))]], group.by = "subclone", pt.size.factor = 45, stroke = NA) + labs(color = "Subclone") + scale_fill_manual(values = spot_colors[-c(1:2)], guide = "none") +
            geom_point(data = edges_ggplot, mapping = aes(x = coord.x, y = coord.y, color = type), shape = 16, size = 8, inherit.aes = F) + scale_color_manual(values = spot_colors) +
            geom_path(data = edges_ggplot %>% arrange(edge, desc(direction)), mapping = aes(x = coord.x, y = coord.y, group = edge), arrow = arrow(length=unit(12,"pt")), colour = "navy", size = 1, inherit.aes = F) +
            theme(legend.position = "right"))
    dev.off()
  }
}


