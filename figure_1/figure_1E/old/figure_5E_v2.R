### title: Plot the copy number profiles of single cells. 
# The following code is used to generate Figure 1E (and other figures)

### To run the code:
# (1) Make sure to load SessionInfo() for correct package version 
# (2) Set both INPUTDIR and OUTPUTDIR to correct path
# (3) Part 4 contains essential code to recreate Figure 1B.
#     Parts can be ignored by setting False in the "if" commands, i.e. "if(F)".

# Computational cost is not too high such that the code should be feasible to 
#   run on local machines.

####################################################################################################################################
### Part 0: Preperation 
####################################################################################################################################

# In case plots look weird, the following setting may fix it
#options(bitmapType='cairo')

### Load libraries
library(tidyverse)
library(ASCAT.sc)
library(parallel)
library(pbmcapply)
library(fastcluster)
library(ComplexHeatmap)
library(grid)
library(gridExtra)
library(magrittr)


INPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_1/data/"
OUTPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_1/results/figure_1e/"
CODEDIR <- "~/Documents/GitHub/MPNST-Phase-1/figure_1/figure_1E/"

### Load metadata
source( paste0(CODEDIR, "metadata/MPNST_common_functions.R") )
source( paste0(CODEDIR, "metadata/10X_DLP_metadata_minimal.R") )

### Prepare data input
#FASTA <- "/camp/project/proj-vanloo/reference_files/human/references/alignment/hs38DH/hs38DH.fa"

scDNA_CN_colour <- c("#00008B", "#7B7BC3", "#FFFFFF", "#FFCCCC", "#FF8080", "#FF3333", "#E60000", "#B31000", "#8B0000", "#660000", "#330000")
names(scDNA_CN_colour) <- 0:10
scDNA_CN_breaks = seq(-0.5,10.5, by = 1)

####################################################################################################################################
### Part 1: Get probes for plotting non-mpcf run
####################################################################################################################################

if (T) {
  data("lSe_filtered_30000.hg38", package = "ASCAT.sc")
  binsize = 500000
  nlSe <- treatlSe(lSe.hg38.filtered[1:23], window = ceiling(binsize/30000))
  
  chr_probe_positions <- do.call(rbind, lapply(1:length(nlSe), function(c){
    return(data.frame(chr = names(nlSe)[c], probe.pos = nlSe[[c]][["starts"]]))
  }))
  chr_probes <- rbind(data.frame(chr = 0, probes = 0), chr_probe_positions %>% 
                        mutate(chr = as.numeric(gsub("X", 23, gsub("chr", "", chr)))) %>% 
                        group_by(chr) %>% 
                        summarise(probes = n())) %>% mutate(cum.probes = cumsum(probes))
}

# ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx))),
#                        col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D")), show_annotation_name = F)

output.dir = OUTPUTDIR
setwd(output.dir)

####################################################################################################################################
### Part 5.1 : Load data
####################################################################################################################################


gamma = 5

### Read in data
as_res <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/all_ASCAT.sc_mpcf_", gamma,".rds"))
asCN_probes <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/MPNST_all_asCN_probes_mpcf_", gamma,".rds"))
asCN_chr_probes <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/MPNST_all_asCN_chr_probes_mpcf_", gamma,".rds"))

#as_res <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/all_ASCAT.sc_ascn_mpcf_", gamma,".rds"))
#asCN_probes <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/MPNST_all_asCN_probes_mpcf_", gamma,".rds"))
#asCN_chr_probes <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/MPNST_all_asCN_chr_probes_mpcf_", gamma,".rds"))

#Note refit results were mostly identical to without refit so used without as a refit put a few cells into WGD state

scDNA_CN_mtx <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/MPNST_all_CN_mtx_mpcf_", gamma,".rds"))

####################################################################################################################################
### Part 5.2 Create heatmap
####################################################################################################################################



#Named kmeans of scDNA_CN_cluster_ids (cleaned kmeans)
png(filename = paste0("MPNST_all_scDNA_CN_mpcf_", gamma,"_K", K, "_named.png"), width = 4000, height = 4000, res = 200)
named_clusters <- unlist(lapply(1:length(scDNA_CN_cluster_ids), function(i) return(rep(names(scDNA_CN_cluster_ids)[i], length(scDNA_CN_cluster_ids[[i]]))))) %>% set_names(unlist(scDNA_CN_cluster_ids))
ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[names(named_clusters),])),
                                       Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx[names(named_clusters),]))),
                       col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
                                  Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")),
                       annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1),
                                                      Tech = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
sc_totCN_heatmap2(CN_mtx = scDNA_CN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, probes = asCN_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), row_ann = ha_row)
dev.off()

####################################################################################################################################
### Part 5.3 Unnecessary other plots
####################################################################################################################################

#Generate raw heatmap
MPNST_CN_hclust_man_ward <- hclust_save_load(scDNA_CN_mtx, sample = paste0("MPNST_all_CN_all_mpcf_", gamma,"_raw"), distance = "manhattan", method = "ward.D2" )
ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx)),
                                       Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx))),
                       col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
                                  Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")), 
                       annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1),
                                                      Tech = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
png(filename = paste0("MPNST_all_scDNA_CN_mpcf_", gamma,"_raw.png"), width = 4000, height = 4000, res = 200)
print(sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx, hclust = MPNST_CN_hclust_man_ward, row_ann = ha_row, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X"))))
dev.off()

K = 70
named_clusters <- cutree(MPNST_CN_hclust_man_ward, k = K) %>% sort()
ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[names(named_clusters),])),
                                       Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx[names(named_clusters),]))),
                       col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
                                  Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")), show_annotation_name = F)
png(filename = paste0("MPNST_all_scDNA_CN_mpcf_", gamma,"_raw_hclust_K", K, ".png"), width = 4000, height = 4000, res = 200)
print(sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), row_ann = ha_row, title = "MPNST All \nTotal CN "))
dev.off()

#Generate pass heatmap
MPNST_CN_hclust_man_ward <- hclust_save_load(scDNA_CN_mtx[all_res$filters,], sample = paste0("MPNST_all_CN_all_mpcf_", gamma,"_pass"), distance = "manhattan", method = "ward.D2" )
ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[all_res$filters,])),
                                       Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx[all_res$filters,]))),
                       col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
                                  Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")), show_annotation_name = F)
png(filename = paste0("MPNST_all_scDNA_CN_mpcf_", gamma,"_pass.png"), width = 4000, height = 4000, res = 200)
print(sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[all_res$filters,], hclust = MPNST_CN_hclust_man_ward, row_ann = ha_row, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), title = "MPNST All \nTotal CN "))
dev.off()

#Filter actually doesn't work well given it is two distributions (10X and DLP)
MPNST_CN_hclust_man_ward <- hclust_save_load(scDNA_CN_mtx, sample = paste0("MPNST_all_CN_all_mpcf_", gamma,"_raw"), distance = "manhattan", method = "ward.D2" )
nonwgd_cells <- names(cutree(MPNST_CN_hclust_man_ward, k = 10)[cutree(MPNST_CN_hclust_man_ward, k = 10) == 1]) #Cluster 2 and 3 is WGD cells

#Check quality of cells fitted with extra WGD
wgd_cells <- names(cutree(MPNST_CN_hclust_man_ward, k = 10)[cutree(MPNST_CN_hclust_man_ward, k = 10) %in% c(2:3)]) #Cluster 2 and 3 is WGD cells

#Plot on nreads/noise plot
if (F) {
  allT <- all_res$allTracks.processed
  allS <- all_res$allSolutions
  getloess <- function(qu, nr) {
    nms <- paste0("n", 1:length(qu))
    names(qu) <- nms
    quo <- qu[order(nr, decreasing = F)]
    fitted <- stats::runmed(quo, k = 31, endrule = "keep")
    names(fitted) <- names(quo)
    list(fitted = fitted[nms], residuals = quo[nms] - 
           fitted[nms])
  }
  getQuality.SD <- function(allT) {
    sapply(allT, function(x) {
      median(abs(diff(unlist(lapply(x$lCTS, function(y) y$smoothed)))))
    })
  }
  nrecords <- sapply(allT, function(x) sum(unlist(lapply(x$lCTS, 
                                                         function(y) y$records))))
  # thresholdNrec <- quantile(nrecords, probs = probs)
  # ambiguous <- sapply(allS, function(x) x$ambiguous)
  # doublet <- sapply(allS, function(x) if (!is.null(x$bestfit)) 
  #   !x$bestfit$ambiguous
  #   else F)
  qualities <- getQuality.SD(allT)
  # thresholdQual <- quantile(qualities, probs = 1 - probs)
  # keep <- qualities <= thresholdQual & nrecords >= thresholdNrec & 
  #   !ambiguous & !doublet
  # keep2 <- !(qualities < thresholdQual & nrecords < thresholdNrec)
  # keep2 <- keep2 & !ambiguous
  # ll <- getloess(qualities[keep2], log2(nrecords)[keep2])
  # ord <- order(log2(nrecords)[keep2], decreasing = F)
  # filters <- (1:length(nrecords)) %in% (which(keep2)[abs(ll$residuals) <= 0.02]) & keep
  
  #Plot
  pdf(paste0("MPNST_WGD_log2nrecord_vs_qualities.pdf"))
  plot(qualities, log2(nrecords), xlab = "Noise logr", ylab = "Total number of reads", pch = 19, 
       cex = ifelse(names(allT) %in% wgd_cells, 0.3, 0.1),
       col = ifelse(names(allT) %in% wgd_cells, "red", "black"))
  dev.off()
}

#Get kmeans clusters
K = 22

png(filename = paste0("MPNST_all_scDNA_CN_mpcf_", gamma,"_K", K, ".png"), width = 4000, height = 4000, res = 200)
if (!file.exists(paste0("MPNST_all_k_means_K",K,"_clusters.rds"))) {
  set.seed(123)
  hm <- ComplexHeatmap::Heatmap(matrix = scDNA_CN_mtx[nonwgd_cells,], cluster_rows = T, cluster_row_slices = FALSE, row_km = K, row_km_repeats = 1000)
  clusters <- row_order(hm)
  k_means_clusters <- lapply(1:length(clusters), function(j) {
    return(rownames(scDNA_CN_mtx[nonwgd_cells,])[clusters[[j]]])
  })
  k_means_clusters_hclust <- lapply(k_means_clusters, function(k) {
    if(length(k) > 1) {
      hclust <- fastcluster::hclust(dist(scDNA_CN_mtx[k,], method = "manhattan"), method = "ward.D2")
      return(hclust[["labels"]][hclust[["order"]]])
    } else {
      return(k)
    }
  }) #perform hclust on each kmeans cluster
  saveRDS(k_means_clusters, paste0("MPNST_all_k_means_K",K,"_clusters.rds"))
  saveRDS(k_means_clusters_hclust, paste0("MPNST_all_k_means_K",K,"_clusters_hclust.rds"))
} else {
  k_means_clusters <- readRDS(paste0("MPNST_all_k_means_K",K,"_clusters_hclust.rds")) #using hclust on each kmeans
}
named_clusters <- unlist(lapply(1:length(k_means_clusters), function(i) return(rep(i, length(k_means_clusters[[i]]))))) %>% set_names(unlist(k_means_clusters))
ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[names(named_clusters),])),
                                       Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx[names(named_clusters),]))),
                       col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
                                  Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")), 
                       annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1),
                                                      Tech = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), row_ann = ha_row)
dev.off()

#Plot profiles of wgd cells and mixed cluster
cluster_name <- "WGD"
idx_to_print <- which(names(all_res$allProfiles) %in% wgd_cells)

cluster_name <- "kmeans_2"
idx_to_print <- which(names(all_res$allProfiles) %in% k_means_clusters[[2]])

pdf(paste0(cluster_name, "_profiles.pdf"), width = 15, height = 4)
for (i in idx_to_print) {
  try({
    plotSolution(all_res$allTracks.processed[[i]], purity = all_res$allSolutions[[i]]$purity,
                 ploidy = all_res$allSolutions[[i]]$ploidy, gamma = 1,
                 sol = all_res$allSolutions[[i]])
    title(names(all_res$allTracks)[i])
  })
}
dev.off()

#Named kmeans of scDNA_CN_cluster_ids (cleaned kmeans)
png(filename = paste0("MPNST_all_scDNA_CN_mpcf_", gamma,"_K", K, "_named.png"), width = 4000, height = 4000, res = 200)
named_clusters <- unlist(lapply(1:length(scDNA_CN_cluster_ids), function(i) return(rep(names(scDNA_CN_cluster_ids)[i], length(scDNA_CN_cluster_ids[[i]]))))) %>% set_names(unlist(scDNA_CN_cluster_ids))
ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[names(named_clusters),])),
                                       Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx[names(named_clusters),]))),
                       col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
                                  Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")),
                       annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1),
                                                      Tech = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), row_ann = ha_row)
dev.off()

# #Plot only 10X or DLP
# subset = "DLP"
# png(filename = paste0("MPNST_", subset, "_scDNA_CN_mpcf_", gamma,"_K", K, ".png"), width = 4000, height = 4000, res = 200)
# k_means_clusters <- readRDS(paste0("MPNST_all_k_means_K",K,"_clusters.rds"))
# named_clusters <- unlist(lapply(1:length(scDNA_CN_cluster_ids), function(i) return(rep(names(scDNA_CN_cluster_ids)[i], length(scDNA_CN_cluster_ids[[i]]))))) %>% set_names(unlist(scDNA_CN_cluster_ids))
# named_clusters <- named_clusters[grep(subset, names(named_clusters))]
# ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[names(named_clusters),])),
#                                        Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx[names(named_clusters),]))),
#                        col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
#                                   Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")),
#                        annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1),
#                                                       Tech = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
# sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), row_ann = ha_row)
# dev.off()