### title: analysis and plotting of scDNA UMAPS
# the code generates Figure 1D in the manuscript (as well as other outputs)

### To run the code:
# (1) Make sure to load SessionInfo() for correct package version
# (2) Set both INPUTDIR, CODEDIR and OUTPUTDIR to correct path
# (3) Only Part 0 and 1 are needed for Figure 1D

# Computational cost is not too high such that the code should be feasible to 
#   run on local machines.

####################################################################################################################################
### Part 0: Code preperation 
####################################################################################################################################

# In case plots look weird, the following setting may fix it
#options(bitmapType='cairo')

### Load libraries
library(tidyverse)
library(pheatmap)
library(pbmcapply)
library(dendextend)
library(umap)

### Define input and output directories
INPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_1/data/"
CODEDIR <- "~/Documents/GitHub/mpnst_phase_1/figure_4/figure_4C/"
OUTPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_4/results/figure_4C/"

### Load metadata
source( paste0(CODEDIR, "metadata/MPNST_common_functions.R") )
source( paste0(CODEDIR, "metadata/10X_DLP_metadata_minimal.R") )
source( paste0(CODEDIR, "metadata/Visium_metadata_minimal.R") )


#source("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DNA/code/MPNST_common_functions.R")
#source(file = "/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/code/10X_DLP_metadata.R")
#source(file = "/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_spatial/code/Visium_metadata.R")

output.dir = OUTPUTDIR
setwd(output.dir)

visium.infercnv.dir <- paste0(INPUTDIR, "10X_spatial/inferCNV_samples_seg_10X_DLP/")

# run <- runs[5,"name"]
# scDNA_CN_cluster_ids <- run_cluster_ids[[run]][["cluster_ids"]]
# scDNA_CN_clusters_coverage <- run_cluster_ids[[run]][["cluster_coverage"]]
# scDNA_CN_cluster_orig_names <- run_cluster_ids[[run]][["cluster_orig_names"]]
# scDNA_normal_ids <- run_cluster_ids[[run]][["normal_ids"]]


gamma = 5
# MPNST_CN_hclust_man_ward <- hclust_save_load(scDNA_CN_mtx, sample = paste0("../ASCAT.sc/MPNST_all_CN_all_mpcf_", gamma,"_raw"), distance = "manhattan", method = "ward.D2" )
# nonwgd_cells <- names(cutree(MPNST_CN_hclust_man_ward, k = 10)[cutree(MPNST_CN_hclust_man_ward, k = 10) == 1]) #Cluster 2 and 3 is WGD cells
scDNA_CN_mtx    <- readRDS(paste0(INPUTDIR,"10X_DLP/ASCAT.sc/MPNST_all_CN_mtx_mpcf_", gamma,".rds"))
all_probes      <- readRDS(paste0(INPUTDIR,"10X_DLP/ASCAT.sc/MPNST_all_probes_mpcf_", gamma,".rds"))
all_chr_probes  <- readRDS(paste0(INPUTDIR,"10X_DLP/ASCAT.sc/MPNST_all_chr_probes_mpcf_", gamma,".rds"))

totCN_mtx_clean <- scDNA_CN_mtx[c(unlist(scDNA_CN_cluster_ids),scDNA_normal_ids),]

####################################################################################################################################
### Part 1: 06/06/22 Generate UMAPs for totalCN from 10X scDNA
####################################################################################################################################
if (F) {
  umap.settings <- umap.defaults
  umap.settings$verbose = T
  umap.settings$metric = "manhattan"
  umap.settings$random_state = 101
  
  scDNA_umap <- readRDS(paste0(INPUTDIR, "10x_DLP/10X_DLP_totCN_umap.rds"))
  
  cluster_ids <- data.frame(barcode = c(unlist(scDNA_CN_cluster_ids), scDNA_normal_ids),
                            Cluster = factor(c(unlist(lapply(1:length(scDNA_CN_cluster_ids), function(c) return(rep(names(scDNA_CN_cluster_ids)[c], length(scDNA_CN_cluster_ids[[c]]))))),
                                                  rep("Diploid", length(scDNA_normal_ids))), levels = c(names(scDNA_CN_cluster_ids), "Diploid")),
                            Data = gsub(".*_(10X|DLP)_.*", "\\1", c(unlist(scDNA_CN_cluster_ids), scDNA_normal_ids)))

  umap_plot <- scDNA_umap[["layout"]] %>% as.data.frame() %>% 
    `colnames<-`(c("UMAP1", "UMAP2")) %>% 
    #rename("UMAP1" = "V1", "UMAP2" = "V2") %>% 
    rownames_to_column("barcode") %>%
    mutate(Region = gsub("_.*", "", rownames(totCN_mtx_clean))) %>% 
    left_join(cluster_ids, by = "barcode")

  #png(filename = paste0("MPNST_10X_DLP_totCN_umap_region.png"), width = 2000, height = 2000, res = 200)
  pdf(file = paste0("MPNST_10X_DLP_totCN_umap_region.pdf"), width = 7, height = 7)
  print(
    umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = Region)) + geom_point(shape = 16, size = 1) + 
      theme_classic(base_size = 24) + 
      theme(legend.position = c(0.9, 0.3),
            axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(color = guide_legend(override.aes = list(size=6)))
  )
  dev.off()
  
  #png(filename = paste0("MPNST_10X_DLP_totCN_umap_cluster.png"), width = 2000, height = 2000, res = 200)
  pdf(file = paste0("MPNST_10X_DLP_totCN_umap_cluster.pdf"), width = 7, height = 7)
  col_len <- length(unique(umap_plot$Cluster))
  print(
    umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = Cluster)) + 
      geom_point(shape = 16, size = 1) + 
      theme_classic(base_size = 24) + 
      theme(legend.position = c(0.8, 0.3),
            axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
      guides(color = guide_legend(ncol=2, override.aes = list(size=6))) +
      scale_color_manual(values = c(hcl(h = seq(15, 375, length = col_len + 1), l = 65, c = 100)[c(seq(1, col_len, 3)[-1], seq(2, col_len, 3), seq(3, col_len, 3))], "grey"))
  )
  dev.off()
  
  # png(filename = paste0("MPNST_10X_DLP_totCN_umap_cluster_R1_1.png"), width = 2000, height = 2000, res = 200)
  # umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = cluster)) + geom_point(shape = 16, size = 1) + theme_classic() + theme(legend.position = c(0.9, 0.8)) + guides(color=guide_legend(ncol=2)) +
  #   scale_color_manual(values = c("coral", rep("grey", col_len-1)))
  # dev.off()
  
  #png(filename = paste0("MPNST_10X_DLP_totCN_umap_tech.png"), width = 2000, height = 2000, res = 200)
  pdf(file = paste0("MPNST_10X_DLP_totCN_umap_tech.pdf"), width = 7, height = 7)
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
  
  umap_plot_labels <- umap_plot %>% group_by(Cluster) %>% summarise(UMAP1=median(UMAP1), UMAP2=median(UMAP2))
  
  #png(filename = paste0("MPNST_10X_DLP_totCN_umap_cluster_lab.png"), width = 2000, height = 2000, res = 200)
  pdf(file = paste0("MPNST_10X_DLP_totCN_umap_cluster_lab.pdf"), width = 7, height = 7)
  print(
    umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = Cluster)) + 
      geom_point(shape = 16, size = 1) + 
      theme_classic(base_size = 24) + 
      theme(legend.position = c(0.8, 0.3),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) + 
      guides(color = guide_legend(ncol=2, override.aes = list(size=6))) +
      scale_color_manual(values = c(hcl(h = seq(15, 375, length = col_len + 1), l = 65, c = 100)[c(seq(1, col_len, 3)[-1], seq(2, col_len, 3), seq(3, col_len, 3))], "grey")) + geom_label(data = umap_plot_labels, aes(label = Cluster))
  )
  dev.off()
  
  #Load bulk profiles
  totCN_chr_probes <- all_probes
  totCN_chr_probes_pos <- do.call(rbind, lapply(1:nrow(totCN_chr_probes), function(seg) {
    return(data.frame(chr = gsub("chr", "", totCN_chr_probes$chr[seg]),
                      startpos = seq(from = totCN_chr_probes$startpos[seg], by = 500000, length.out = totCN_chr_probes$n.probes[seg]),
                      endpos = seq(from = totCN_chr_probes$startpos[seg] + 499999, by = 500000, length.out = totCN_chr_probes$n.probes[seg]), stringsAsFactors = F))
  }))
  
  if (!file.exists("bb_profiles.rds")) {
    bb_profiles <- do.call(rbind, lapply(bb_subclones, function(CN_profile) {
      totCN <- unlist(lapply(1:nrow(totCN_chr_probes_pos), function(p) {
        CN_seg <- CN_profile %>% filter(chr == totCN_chr_probes_pos[p,"chr"], startpos <=  totCN_chr_probes_pos[p,"startpos"], endpos >= totCN_chr_probes_pos[p,"endpos"])
        return(ifelse(nrow(CN_seg) == 0,
                      2, #Set segments missing in BB to CN2
                      ifelse(CN_seg$frac1_A >= 0.5, CN_seg$nMaj1_A + CN_seg$nMin1_A, CN_seg$nMaj2_A + CN_seg$nMin2_A))) 
      }))
      return(totCN)
    }))
    saveRDS(bb_profiles, "bb_profiles.rds")
  } else {
    bb_profiles <- readRDS("bb_profiles.rds")
  }

  #Project onto UMAP
  if (!file.exists(paste0(output.dir, "10X_DLP_bulk_umap.rds"))) {
    scDNA_10X_bulk_umap <- predict(scDNA_umap, bb_profiles)
    saveRDS(scDNA_10X_bulk_umap, paste0(output.dir, "10X_DLP_bulk_umap.rds"))
  } else {
    scDNA_10X_bulk_umap <- readRDS(paste0(output.dir, "10X_DLP_bulk_umap.rds"))
  }
  
  scDNA_10X_bulk_umap_plot <- rbind(umap_plot, 
                                      data.frame(barcode = rownames(scDNA_10X_bulk_umap),
                                                 UMAP1 = scDNA_10X_bulk_umap[,1],
                                                 UMAP2 = scDNA_10X_bulk_umap[,2],
                                                 Region = names(bb_subclones),
                                                 Cluster = paste0("Bulk_", names(bb_subclones)),
                                                 Data = "Bulk") %>% remove_rownames())
  
  png(filename = paste0("MPNST_10X_DLP_bulk_totCN_umap_region.png"), width = 2000, height = 2000, res = 200)
  print(
    scDNA_10X_bulk_umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = Region, size = Data, alpha = Data)) + geom_point(shape = 16) + 
      scale_size_manual(values = c(1, 1, 20)) + scale_alpha_manual(values = c(1, 1, 0.3)) + theme_classic(base_size = 16) + theme(legend.position = c(0.8, 0.3), legend.box = "horizontal")
  )
  dev.off()
  
  pdf(file = paste0("MPNST_10X_DLP_bulk_totCN_umap_region.pdf"), width = 7, height = 7)
  print(
    scDNA_10X_bulk_umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = Region, size = Data, alpha = Data)) + 
      geom_point(shape = 16) + 
      scale_size_manual(values = c(1, 1, 12)) + 
      scale_alpha_manual(values = c(1, 1, 0.3)) + 
      theme_classic(base_size = 24) + 
      theme(legend.position = c(0.75, 0.3), legend.box = "horizontal",
            axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
      guides(color = guide_legend(override.aes = list(size=6)))
  )
  dev.off()
}

####################################################################################################################################
### Part 2: 06/06/22 Project Visium totalCN onto UMAP
####################################################################################################################################

rm_small_cluster <- 0 #Set to 50 or 0
adj_cluster <- F
high_res <- F
suffix <- paste0(ifelse(rm_small_cluster > 0, paste0("_rm_sub", rm_small_cluster), ""), ifelse(adj_cluster, "_adj", ""), ifelse(high_res, "_res", ""))
R1_ref_CN <- readRDS(paste0(visium.infercnv.dir, "MPNST_R1_ref_CN.rds"))

#Function to discretize
disc_mtx <- function(mtx, breaks, labels) {
  d_mtx <- matrix(as.matrix(cut(mtx, breaks = breaks, labels = labels)) %>% as.numeric(), nrow(mtx))
  rownames(d_mtx) <- rownames(mtx)
  colnames(d_mtx) <- colnames(mtx)
  return(d_mtx)
}

discretize_breaks <- c(-Inf,-0.2,-0.05,0.05,0.2,Inf)
discretize_labels <- c(-2,-1,0,1,2)

##Import Visium 
if (!file.exists(paste0(output.dir, "visium_totCN_profiles.rds"))) {
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
  saveRDS(visium_totCN_profiles, paste0("visium_totCN_profiles.rds"))
} else {
  visium_totCN_profiles <- readRDS(paste0(output.dir, "visium_totCN_profiles.rds"))
}

#Load stored UMAP
scDNA_umap <- readRDS(paste0(output.dir, "10X_DLP_totCN_umap.rds"))
cluster_ids <- data.frame(barcode = c(unlist(scDNA_CN_cluster_ids), scDNA_normal_ids),
                          Cluster = factor(c(unlist(lapply(1:length(scDNA_CN_cluster_ids), function(c) return(rep(c, length(scDNA_CN_cluster_ids[[c]]))))),
                                             rep("Diploid", length(scDNA_normal_ids))), levels = c(1:length(scDNA_CN_cluster_ids), "Diploid")),
                          Data = gsub(".*_(10X|DLP)_.*", "\\1", c(unlist(scDNA_CN_cluster_ids), scDNA_normal_ids)))

# umap_plot <- scDNA_umap[["layout"]] %>% as.data.frame() %>% rename("UMAP1" = "V1", "UMAP2" = "V2") %>% rownames_to_column("barcode") %>%
#   mutate(Region = gsub("_.*", "", rownames(totCN_mtx_clean)), Section = " ") %>% left_join(cluster_ids, by = "barcode")
umap_plot <- scDNA_umap[["layout"]] %>% as.data.frame() %>% 
  rename("UMAP1" = "V1", "UMAP2" = "V2") %>% 
  rownames_to_column("barcode") %>%
  mutate(Region = gsub("_.*", "", rownames(totCN_mtx_clean)), Section = " ") %>% 
  left_join(cluster_ids, by = "barcode")

#Project onto UMAP
if (!file.exists(paste0(output.dir, "10X_DLP_visium_umap.rds"))) {
  scDNA_10X_visium_umap <- predict(scDNA_umap, visium_totCN_profiles)
  saveRDS(scDNA_10X_visium_umap, paste0(output.dir, "10X_DLP_visium_umap.rds"))
} else {
  scDNA_10X_visium_umap <- readRDS(paste0(output.dir, "10X_DLP_visium_umap.rds"))
}

scDNA_10X_visium_umap_plot <- rbind(umap_plot, 
                                    data.frame(barcode = rownames(scDNA_10X_visium_umap),
                                               UMAP1 = scDNA_10X_visium_umap[,1],
                                               UMAP2 = scDNA_10X_visium_umap[,2],
                                               Region = gsub("a|b|_.*", "", rownames(visium_totCN_profiles)),
                                               Section = gsub("_.*", "", rownames(visium_totCN_profiles)),
                                               Cluster = rownames(visium_totCN_profiles),
                                               Data = "Visium") %>% remove_rownames())

# png(filename = paste0("MPNST_10X_Visium_totCN_umap_region.png"), width = 2000, height = 2000, res = 200)
# scDNA_10X_visium_umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = Region, size = Data, alpha = Data)) + geom_point() + 
#   scale_size_manual(values = c(0.5, 0.5, 8)) + scale_alpha_manual(values = c(1, 1, 0.2)) + 
#   theme_classic(base_size = 16) + theme(legend.position = c(0.8, 0.3), legend.box = "horizontal")
# dev.off()

pdf(file = paste0("MPNST_10X_Visium_totCN_umap_region.pdf"), width = 7, height = 7)
print(
  scDNA_10X_visium_umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = Region, size = Data, alpha = Data)) + geom_point() + 
    scale_size_manual(values = c(0.5, 0.5, 8)) + scale_alpha_manual(values = c(1, 1, 0.2)) + 
    theme_classic(base_size = 24) + theme(legend.position = c(0.77, 0.3), legend.box = "horizontal",
                                          axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
    guides(color = guide_legend(override.aes = list(size=6)))
)
dev.off()

scDNA_10X_visium_umap_plot <- scDNA_10X_visium_umap_plot %>% mutate(Section = ifelse(Data == "Visium", Section, ifelse(Region %in% c("P", "R5"), paste0(Region,"a"), Region)))

# png(filename = paste0("MPNST_10X_Visium_totCN_umap_section.png"), width = 2000, height = 2000, res = 200)
# print(
#   scDNA_10X_visium_umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = Section, size = Data, alpha = Data)) + geom_point() +
#   scale_color_manual(values = c(region_colours[1], "darkgoldenrod1", region_colours[2:6], "darkorchid4")) +
#   scale_size_manual(values = c(0.5, 0.5, 8)) + scale_alpha_manual(values = c(1, 1, 0.2)) + theme_classic(base_size = 16) + theme(legend.position = c(0.8, 0.3), legend.box = "horizontal")
# )
# dev.off()


# pdf(file = paste0("MPNST_10X_Visium_totCN_umap_section.pdf"), width = 7, height = 7)
# print(
#   scDNA_10X_visium_umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = Section, size = Data, alpha = Data)) + geom_point() + 
#     scale_color_manual(values = c(region_colours[1], "darkgoldenrod1", region_colours[2:6], "darkorchid4")) + 
#     scale_size_manual(values = c(0.5, 0.5, 8)) + scale_alpha_manual(values = c(1, 1, 0.2)) + 
#     theme_classic(base_size = 24) + theme(legend.position = c(0.77, 0.3), legend.box = "horizontal", axis.text.x=element_blank(),
#                                           axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
#     guides(color = guide_legend(override.aes = list(size=6)))
# )
# dev.off()

pdf(file = paste0("MPNST_10X_Visium_totCN_umap_section.pdf"), width = 7, height = 7)
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

