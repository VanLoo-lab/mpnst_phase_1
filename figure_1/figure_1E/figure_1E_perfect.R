### title: Plot the copy number profiles of single cells. 
# The following code is used to generate Figure 1E (and other figures)

### To run the code:
# (1) Make sure to load SessionInfo() for correct package version 
# (2) Set both INPUTDIR and OUTPUTDIR to correct path
# (3) Part 4 contains essential code to recreate Figure 1B.
#     Parts can be ignored by setting False in the "if" commands, i.e. "if(F)".

# Computational cost is not too high such that the code should be feasible to 
#   run on local machines with sufficient RAM

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
### Part 6
####################################################################################################################################

library(ComplexHeatmap)
library(circlize)
library(tibble)
library(dplyr)
library(grid)

gamma = 5
as_res <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/all_ASCAT.sc_ascn_mpcf_", gamma,".rds"))
asCN_probes <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/MPNST_all_asCN_probes_mpcf_", gamma,".rds"))
asCN_chr_probes <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/MPNST_all_asCN_chr_probes_mpcf_", gamma,".rds"))

#Generate asCN heatmap
scDNA_asCN_mtx <- exist_load_save(file = "MPNST_DLP_asfree_asCN_mtx.rds",
                                  do.call(rbind, lapply(which(lengths(as_res[["allProfiles_AS"]])!=0), function(b){
                                    prof <- as_res$allProfiles_AS[[b]]$nprof.fixed[!is.na(as_res$allProfiles_AS[[b]]$nprof.fixed$total_copy_number),] #Using fixed
                                    prof$nA[is.na(prof$nA)] <- 0
                                    prof$nB[is.na(prof$nB)] <- prof$total_copy_number[which(is.na(prof$nB))]
                                    asCN_df <- data.frame(allele1 = prof$nA,
                                                          allele2 = prof$nB) %>% #flip alleles so A is larger
                                      mutate(alleleA = ifelse(allele1>=allele2, allele1, allele2),
                                             alleleB = ifelse(allele1>=allele2, allele2, allele1)) %>%
                                      mutate(alleleA = ifelse(alleleA<0, 0, alleleA), #code to fix CN of -1
                                             alleleB = ifelse(alleleB<0, 0, alleleB)) %>%
                                      mutate(alleleT = alleleA + alleleB) %>% #get total CN
                                      mutate(plotT = alleleT*10+alleleB) %>% select(-c(1:2))#Convert to matrix format
                                    rep(na.omit(as.numeric(asCN_df$plotT)), asCN_probes$n.probes) #NAs are omitted...
                                  })) %>% set_rownames(names(as_res[["allProfiles_AS"]])))


# MPNST_asCN_hclust_ward <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/MPNST_all_asCN_hclust_ward.rds"))

####################################################################################################################################
### Part 7
####################################################################################################################################

# --- 1) build asCN matrix (rows = cells, cols = genome bins) ---
# scDNA_asCN_mtx <- do.call(rbind, lapply(which(lengths(as_res[["allProfiles_AS"]]) != 0), function(b){
#   
#   prof <- as_res$allProfiles_AS[[b]]$nprof.fixed
#   prof <- prof[!is.na(prof$total_copy_number), ]
#   
#   prof$nA[is.na(prof$nA)] <- 0
#   prof$nB[is.na(prof$nB)] <- prof$total_copy_number[is.na(prof$nB)]
#   
#   asCN_df <- data.frame(allele1 = prof$nA, allele2 = prof$nB) %>%
#     mutate(
#       alleleA = ifelse(allele1 >= allele2, allele1, allele2),
#       alleleB = ifelse(allele1 >= allele2, allele2, allele1),
#       alleleA = pmax(alleleA, 0),
#       alleleB = pmax(alleleB, 0),
#       alleleT = alleleA + alleleB,
#       plotT   = alleleT * 10 + alleleB
#     ) %>%
#     select(plotT)
#   
#   rep(na.omit(as.numeric(asCN_df$plotT)), asCN_probes$n.probes)
# })) %>% set_rownames(names(as_res[["allProfiles_AS"]]))

# --- 2) clustering (ward) ---
#MPNST_asCN_hclust_ward <- fastcluster::hclust(dist(scDNA_asCN_mtx), method = "ward.D2")
MPNST_asCN_hclust_ward <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/MPNST_all_asCN_hclust_ward.rds"))


# --- 3) row annotation (Region + Tech) ---
Region <- factor(gsub("_.*", "", rownames(scDNA_asCN_mtx)), levels = c("R1","R2","R3","R4","R5","P"))
Tech   <- factor(gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_asCN_mtx)), levels = c("10X","DLP"))

ha_row <- rowAnnotation(
  df = data.frame(Region = Region, Tech = Tech),
  col = list(
    Region = c("R1"="#B79F00","R2"="#00BA38","R3"="#00BFC4","R4"="#619CFF","R5"="#F564E3","P"="#F8766D"),
    Tech   = c("10X"="mediumpurple1","DLP"="olivedrab3")
  ),
  show_annotation_name = FALSE
)

# --- 4) probes for chromosome boundaries (named 1:22,X) ---
probes <- asCN_chr_probes$cum.probes[-1]
names(probes) <- c(as.character(1:22), "X")

# --- 5) draw heatmap (writes PNG) ---
png("MPNST_all_asfixed_asCN_heatmap.png", width = 4000, height = 4000, res = 200)
sc_asCN_heatmap2(
  CN_mtx = scDNA_asCN_mtx,
  hclust = MPNST_asCN_hclust_ward,
  row_ann = ha_row,
  probes = probes
)
dev.off()

# --- 6) draw heatmap with clusters

#k_means_clusters <- readRDS(paste0("MPNST_all_k_means_K",K,"_clusters_hclust.rds")) #using hclust on each kmeans
k_means_clusters <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/MPNST_all_k_means_K", K,"_clusters_hclust.rds"))

ha_row = rowAnnotation(
  df = data.frame(Region = gsub("_.*", "", rownames(scDNA_asCN_mtx[names(named_clusters),])), 
                  Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_asCN_mtx[names(named_clusters),]))),
  col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"), 
             Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")),
  annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 22), 
                                               grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1), 
                                 Tech = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 22), grid_height = unit(0.8, "cm"), 
                                             grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
png(filename = paste0("MPNST_all_asfixed_asCN_K", K, "_named.png"), width = 4000, height = 4000, res = 200)
named_clusters <- unlist(lapply(1:length(scDNA_CN_cluster_ids), function(i) return(rep(names(scDNA_CN_cluster_ids)[i], length(scDNA_CN_cluster_ids[[i]]))))) %>% set_names(unlist(scDNA_CN_cluster_ids))
print(sc_asCN_heatmap3(CN_mtx = scDNA_asCN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, row_ann = ha_row, column_title = NA, probes = asCN_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X"))))
dev.off()


