### title: Plot the copy number profiles of single cells. 
# The following code is used to generate Figure 1E (and other figures)

### To run the code:
# (1) Make sure to load SessionInfo() for correct package version 
# (2) Set both INPUTDIR and OUTPUTDIR to correct path
# (3) ...

# Computational cost is not too high such that the code should be feasible to 
#   run on local machines with sufficient RAM

####################################################################################################################################
### Part 0: Preperation 
####################################################################################################################################

# In case plots look weird, the following setting may fix it
#options(bitmapType='cairo')

### Load libraries
library(tidyverse)
library(ape)
library(phangorn)
library(treeio)
library(ggtree)
library(gridExtra)

### Define input, code and output dir
INPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_1/data/"
OUTPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_2/results/figure_2c/"
CODEDIR <- "~/Documents/GitHub/MPNST-Phase-1/figure_2/figure_2C/"

#Load metadata
#source(file = "/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/code/10X_DLP_metadata.R")
#source(file = "/camp/project/proj-vanloo/analyses/hyan/mpnst/LCM/code/LCM_metadata.R")

source(paste0(CODEDIR,"metadata/10X_DLP_metadata_minimal.R"))
source(paste0(CODEDIR,"metadata/10X_DLP_metadata_minimal.R"))



output.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/results/one_tree_LCM_GnT/")
setwd(output.dir)

####################################################################################################################################
#Choose which clustering run to use
####################################################################################################################################
filter_n_cells <- 2 #SNVs with less than or equal to filter_n_cells cell removed 
sub_folder <- paste0(run, "_filter_n_", filter_n_cells)
system(paste0("mkdir -p ", output.dir, sub_folder))

#Generate trees (only need to run once)
if (F) {
  scaling_ratio = 10
  SNV_medicc_tree <- read.tree(file = "../one_tree/22_kmeans_10X_DLP_filter_n_2/SNV_medicc_tree_scale10.new")
  
  ####################################################################################################################################
  #Add LCM tips
  ####################################################################################################################################
  medicc_pdm <- read.delim(paste0("../MEDICC2/MPNST_all_single_cells_LCM_GnT_2.5/MPNST_all_single_cells_LCM_GnT_2.5_ds_6310_wgd_seed100/MPNST_all_single_cells_LCM_GnT_2.5_ds_6310_seed100_pdm_total.tsv"), row.names = 1)
  medicc_sc_pdm <- medicc_pdm[,gsub("-1", ".1", SNV_medicc_tree[["tip.label"]][-1])]
  LCM_medicc_spots <- LCM_coordinates_tumour %>% mutate(sample = gsub("_.*", "", sample)) %>% pull(sample)
  LCM_medicc_spots <- LCM_medicc_spots[LCM_medicc_spots %in% gsub("_.*", "", rownames(medicc_sc_pdm))]
  
  # SC_medicc_cells <- paste0(SNV_medicc_tree[["tip.label"]][-1], ".1")
  GnT_medicc_cells <- grep("EGA", rownames(medicc_pdm), value = T)
  GnT_medicc_tumour_cells <- data.frame(label = GnT_medicc_cells, Barcode = gsub(".*_", "", GnT_medicc_cells)) %>%
    left_join (readRDS("/camp/project/proj-vanloo/analyses/hyan/mpnst/GnT/results/Genotype_SNPs/both_haplo_count_combined_by_cell.rds")) %>% 
    mutate(Haplotype_Ratio = Haplo_1_Count/Haplo_2_Count) %>% filter(Haplotype_Ratio > 2) %>% pull(label)
  # medicc_sc_pdm <- medicc_pdm[,!colnames(medicc_pdm) %in% c(LCM_medicc_spots, GnT_medicc_cells)]
  
  #Edit tree as some lengths are 0 and can't add tips so changed to 0.1
  SNV_medicc_tree_add <- as_tibble(SNV_medicc_tree) %>% 
    mutate(branch.length = ifelse(branch.length == 0 & !str_detect(label, "internal"), 0.1, branch.length)) %>% as.phylo()
  
  for (b in GnT_medicc_tumour_cells) {
    closest_cell <- sort(medicc_sc_pdm[b,])[1] %>% names()
    closest_cell_distance <- sort(medicc_sc_pdm[b,])[1] %>% as.numeric()
    closest_cell_node <- SNV_medicc_tree_add %>% as_tibble() %>% filter(label == gsub("\\.1", "-1", closest_cell)) %>% pull(node)
    SNV_medicc_tree_add <- bind.tree(SNV_medicc_tree_add, read.tree(text = paste0("(",b,":",closest_cell_distance,");")), where = closest_cell_node, position = 0.1)
  }
  write.tree(SNV_medicc_tree_add, paste0(sub_folder, "/SNV_medicc_10X_GnT_tree_scale",scaling_ratio,".new"))
  
  SNV_medicc_tree_add <- as_tibble(SNV_medicc_tree) %>% 
    mutate(branch.length = ifelse(branch.length == 0 & !str_detect(label, "internal"), 0.1, branch.length)) %>% as.phylo()
  
  for (b in LCM_medicc_spots) {
    closest_cell <- sort(medicc_sc_pdm[b,])[1] %>% names()
    closest_cell_distance <- sort(medicc_sc_pdm[b,])[1] %>% as.numeric()
    closest_cell_node <- SNV_medicc_tree_add %>% as_tibble() %>% filter(label == gsub("\\.1", "-1", closest_cell)) %>% pull(node)
    SNV_medicc_tree_add <- bind.tree(SNV_medicc_tree_add, read.tree(text = paste0("(",b,":",closest_cell_distance,");")), where = closest_cell_node, position = 0.1)
  }
  write.tree(SNV_medicc_tree_add, paste0(sub_folder, "/SNV_medicc_10X_LCM_tree_scale",scaling_ratio,".new"))
  
  for (b in GnT_medicc_tumour_cells) {
    closest_cell <- sort(medicc_sc_pdm[b,])[1] %>% names()
    closest_cell_distance <- sort(medicc_sc_pdm[b,])[1] %>% as.numeric()
    closest_cell_node <- SNV_medicc_tree_add %>% as_tibble() %>% filter(label == gsub("\\.1", "-1", closest_cell)) %>% pull(node)
    SNV_medicc_tree_add <- bind.tree(SNV_medicc_tree_add, read.tree(text = paste0("(",b,":",closest_cell_distance,");")), where = closest_cell_node, position = 0.1)
  }
  
  write.tree(SNV_medicc_tree_add, paste0(sub_folder, "/SNV_medicc_10X_GnT_LCM_tree_scale",scaling_ratio,".new"))
}

####################################################################################################################################
#Make prettier plot
####################################################################################################################################
scaling_ratio = 10
tree = "GnT_LCM" #GnT, LCM or GnT_LCM

if (T) {
  GnT_samples = c("R1", "R2", "R3", "R4", "R5", "P")
  names(GnT_samples) = c("26787_8", "30095_2", "30363_3", "30177_3", "30177_4", "30177_2")
  
  download.dir <- "/camp/project/proj-vanloo/analyses/mtarabichi/sarcoma/mpnst/GnT/dl_from_ega/"
  GnT.input.dir <- paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/GnT/data/hg38_bwa_BAM/")
  
  GnT_IDs <- list.files(download.dir)[-c(1:2)]
  GnT_file_names <- unlist(lapply(GnT_IDs, function(i) {
    return(gsub("\\..*", "", list.files(paste0(download.dir, i))) %>% unique())
  }))
  
  all_odd_bams <- dir(GnT.input.dir, pattern = "[13579]_markdup.bam$")
  names(all_odd_bams) <- GnT_file_names[seq(2,length(GnT_file_names),2)]
  
  region_odd_IDs <- lapply(names(GnT_samples), function(s) {
    return(gsub("_markdup.bam", "", all_odd_bams[grep(s, names(all_odd_bams))]))
  })
  names(region_odd_IDs) <- GnT_samples
  
  region_bams <- lapply(region_odd_IDs, function(b) {
    return(all_odd_bams[all_odd_bams %in% paste0(b,"_markdup.bam")])
  })
}
raw_tree <- read_file(paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_scale",scaling_ratio,".new"))
for (b in unlist(region_odd_IDs)) {
  print(b)
  raw_tree <- gsub(b, paste0(samples[grep(b, region_odd_IDs)], "_", b), raw_tree)
}
SNV_medicc_tree_add <- read.tree(text = raw_tree)

SNV_medicc_tree_metadata <- data.frame(label=SNV_medicc_tree_add$tip.label, region = gsub("_.*", "", SNV_medicc_tree_add$tip.label), stringsAsFactors = F) %>%
  left_join(LCM_coordinates_tumour %>% select(sample, Cell_Type, side) %>% mutate(sample = gsub("_.*", "", sample), side = gsub(".*-", "", side)), by = c("label" = "sample")) %>%
  mutate(region = ifelse(is.na(Cell_Type), region, Cell_Type)) %>% mutate(Cell_Type = ifelse(str_detect(label, "EGAF"), "G&T", ifelse(is.na(Cell_Type), "10X", "LCM"))) %>% rename(Type = Cell_Type, Side = side) %>%
  mutate(LCM_region = ifelse(Type == "LCM", region , Type)) %>% mutate(Side = ifelse(Type == "LCM", Side, "")) %>% mutate(GnT_region = ifelse(Type == "G&T", region , Type)) #Add LCM/GnT region info 
SNV_medicc_tree_groups <- full_join(SNV_medicc_tree_add, SNV_medicc_tree_metadata, by = "label") 

SNV_medicc_tree_colours <- c("grey", region_colours[names(region_colours) %in% SNV_medicc_tree_groups@data[["region"]]])

SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups) + geom_tippoint(aes(shape = Type, colour = region), size = .5) + #geom_text(aes(label=node)) + #geom_tiplab(size = 8) +
  scale_color_manual(values = SNV_medicc_tree_colours, labels = c("diploid", samples)) +
  theme_tree2()

#Exploratory Plots
if (F) {
  pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,".pdf"), width = 10, height = 14)
  print(SNV_medicc_ggtree)
  dev.off()
  
  pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_lab_reorder.pdf"), width = 10, height = 56)
  (SNV_medicc_ggtree + geom_text(aes(label = node)))
  dev.off()
  
  pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_reorder.pdf"), width = 10, height = 14)
  SNV_medicc_ggtree %>% flip(12041, 6314)
  dev.off()
}

#Annotate 
cluster_colours_reorder <- c("black", "darkslategray4", "darkgreen", "#00BA38", "#619CFF", "#00BFC4", "coral", "#B79F00", "#F564E3", "#F8766D")
cluster_names <- c("XXXXXX", "0XXX00", "0X0X00", "0X0000", "000X00", "00X000", "X000X0", "X00000", "0000X0", "00000X")
node_names <- c("MRCA", "P", "P_1", "P_2", "P_3",
                "R1", "R1_1", "R1_2", "R1_3", "R1_4", "R1_R5",
                "R2", "R2_1", "R2_2", "R2_3", "R2_R3_R4", "R2_R4",
                "R3_1", "R4", "R4_1", "R4_2", "R4_3",
                "R5", "R5_1", "R5_2", "R5_3", "R5_4")

node_colours <- c("white", rep(cluster_colours_reorder[10],4),
                  rep(cluster_colours_reorder[8],5), cluster_colours_reorder[7],
                  rep(cluster_colours_reorder[4],4), cluster_colours_reorder[2], cluster_colours_reorder[3],
                  cluster_colours_reorder[6], rep(cluster_colours_reorder[5],4), rep(cluster_colours_reorder[9],5))
names(node_colours) <- node_names

SNV_medicc_info <- data.frame(node = unlist(lapply(node_names, function(n) {which(as_tibble(SNV_medicc_tree_groups)$label == n)})),
                              name = node_names,
                              fill_col = node_colours,
                              SNV_branch = "SNV_T", stringsAsFactors = F) %>% arrange(name) #Anndata for trunk branches

SNV_medicc_info <- rbind(SNV_medicc_info, anti_join(data.frame(node = SNV_medicc_tree_groups@data$node,
                                                               name = NA,
                                                               fill_col = NA,
                                                               SNV_branch = "SNV_F") ,SNV_medicc_info, by = "node")) #Add non trunk branches

SNV_medicc_info[SNV_medicc_info$node %in% c(1),"SNV_branch"] <- "SNV_T" #Change root branch to SNV too

# SNV_medicc_info <- data.frame(node = c(3262, 3263, 3456, 3264, 3370, 3431,
#                                        3651),
#                               name = c("MRCA", "P", "P_4", "P_1", "P_2", "P_3",
#                                        "R1_R5"),
#                               fill_col = c("white", rep(cluster_colours_reorder[10],5),
#                                            cluster_colours_reorder[7]), stringsAsFactors = F)
SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups)
#Edit so Fertilised Egg appears
# SNV_medicc_ggtree[["data"]][which(SNV_medicc_ggtree[["data"]]$label == "Fertilised_egg"), "isTip"] <- FALSE

# pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"tree_ggtree_scale",scaling_ratio,"_nodelab_reorder.pdf"), width = 10, height = 14)
# SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = region), size = .5) + 
#   scale_color_manual(values = SNV_medicc_tree_colours, labels = c("diploid", samples)) +
#     geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_info$fill_col,  labels=SNV_medicc_info$name) +
#   theme_tree2()
# dev.off()

if (tree == "LCM") {
  pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder.pdf"), width = 10, height = 14)
  SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(color = SNV_branch, size = SNV_branch)) %>% flip(11683, 6135)
  SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = region), size = .5) + 
    scale_color_manual(values = c(SNV_medicc_tree_colours,"black", "black"), labels = c("diploid", samples, "yes")) + scale_size_manual(values = c(0.3,2), labels = c("CNA", "SNV")) +
    geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_info$fill_col,  labels=SNV_medicc_info$name) +
    guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2()
  dev.off()
  
  #Hightlight LCM spots
  pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder_LCM.pdf"), width = 10, height = 14)
  SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(color = SNV_branch, size = SNV_branch)) %>% flip(11683, 6135)
  SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = LCM_region, shape = Type), size = .5) + 
    scale_color_manual(values = c(SNV_medicc_tree_colours[c(1,2,3,6)],"black", "black"), labels = c("diploid", samples, "yes")) + 
    scale_size_manual(values = c(0.3,2), labels = c("CNA", "SNV")) +
    scale_shape_manual(values = c(19,15)) +
    geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_info$fill_col,  labels=SNV_medicc_info$name) +
    guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2()
  dev.off()
  
  pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder_side.pdf"), width = 10, height = 14)
  SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(color = SNV_branch, size = SNV_branch)) %>% flip(11683, 6135)
  SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = Side, shape = Type), size = .5, show.legend = T) + 
    scale_color_manual(values = c("grey", hcl(h = seq(15, 375, length = 3 + 1), l = 65, c = 100)[-4], "black", "black"), labels = c("", "Back", "Front", "Side")) +
    scale_size_manual(values = c(0.3,2), labels = c("CNA", "SNV")) +
    scale_shape_manual(values = c(19,15)) +
    geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_info$fill_col,  labels=SNV_medicc_info$name) +
    guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2()
  dev.off()
  
  #Prune
  for (sample in c("P", "R1", "R4")) {
    plot_height = round(length(unlist(Descendants(SNV_medicc_tree_add,  SNV_medicc_info[which(SNV_medicc_info$name == sample), "node"], type = "tips")))/length(SNV_medicc_tree_add[["tip.label"]])*35)
    
    # Old code using viewCLade, not very clean
    # pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder_side_", sample, ".pdf"), width = 10, height = plot_height)
    # SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(color = SNV_branch, size = SNV_branch)) %>% flip(11683, 6135)
    # print(viewClade(SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = Side, shape = Type), size = 1) + 
    #                   scale_color_manual(values = c("grey", hcl(h = seq(15, 375, length = 3 + 1), l = 65, c = 100)[-4], "black", "black"), labels = c("Front", "Side", "Back", "")) +
    #                   scale_size_manual(values = c(0.3,2), labels = c("CNA", "SNV")) +
    #                   scale_shape_manual(values = c(19,15)) +
    #                   geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_info$fill_col,  labels=SNV_medicc_info$name) +
    #                   guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2(), 
    #                 SNV_medicc_info[which(SNV_medicc_info$name == sample), "node"]))
    # dev.off()
    # 
    # pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder_LCM_", sample, ".pdf"), width = 10, height = plot_height)
    # SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(color = SNV_branch, size = SNV_branch)) %>% flip(11683, 6135)
    # print(viewClade(SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = LCM_region, shape = Type), size = 1) + 
    #                   scale_color_manual(values = c(SNV_medicc_tree_colours[c(1,2,3,6)],"black", "black"), labels = c("diploid", samples, "yes")) + 
    #                   scale_size_manual(values = c(0.3,2), labels = c("CNA", "SNV")) +
    #                   scale_shape_manual(values = c(19,15)) +
    #                   geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_info$fill_col,  labels=SNV_medicc_info$name) +
    #                   guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2(), 
    #                 SNV_medicc_info[which(SNV_medicc_info$name == sample), "node"]))
    # dev.off()
    
    SNV_medicc_ggtree_sub <- tree_subset(tree = SNV_medicc_tree_groups, node = sample, levels_back = 0)
    if (sample == "P") {
      sub_node_names <- c("P", "P_1", "P_2", "P_3")
    } else {
      if (sample == "R1") {
        sub_node_names <- c("R1", "R1_2", "R1_3", "R1_4")
      } else{
        sub_node_names <- c("R4", "R4_1", "R4_2", "R4_3")
      }
    }
    
    SNV_medicc_sub_info <- data.frame(node = unlist(lapply(sub_node_names, function(n) {which(as_tibble(SNV_medicc_ggtree_sub)$label == n)})),
                                      name = sub_node_names,
                                      fill_col = node_colours[sub_node_names],
                                      SNV_branch = "SNV_T", stringsAsFactors = F) %>% arrange(name) #Anndata for trunk branches
    SNV_medicc_sub_info <- rbind(SNV_medicc_sub_info, anti_join(data.frame(node = SNV_medicc_ggtree_sub@data$node,
                                                                           name = NA,
                                                                           fill_col = NA,
                                                                           SNV_branch = "SNV_F") ,SNV_medicc_sub_info, by = "node")) #Add non trunk branches
    # SNV_medicc_sub_info[SNV_medicc_sub_info$node %in% c(1),"SNV_branch"] <- "SNV_T" #Change root branch to SNV too #don't need this when pruned
    
    pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder_side_", sample, ".pdf"), width = 10, height = plot_height)
    SNV_medicc_ggtree <- ggtree(SNV_medicc_ggtree_sub, aes(color = SNV_branch, size = SNV_branch))
    print(SNV_medicc_ggtree %<+% SNV_medicc_sub_info + geom_tippoint(aes(colour = Side, shape = Type, size = Type)) + xlim(0,150) + 
            scale_color_manual(values = c("grey", hcl(h = seq(15, 375, length = 3 + 1), l = 65, c = 100)[-4], "black", "black"), labels = c("", "Back", "Front", "Side")) +
            scale_size_manual(values = c(1,3,0.3,2), labels = c("10X", "LCM", "SNV", "CNA")) +
            scale_shape_manual(values = c(19,15)) +
            geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_sub_info$fill_col,  labels=SNV_medicc_sub_info$name) +
            guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2(), 
          SNV_medicc_sub_info[which(SNV_medicc_sub_info$name == sample), "node"]) 
    dev.off()
    
    pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder_LCM_", sample, ".pdf"), width = 10, height = plot_height)
    SNV_medicc_ggtree <- ggtree(SNV_medicc_ggtree_sub, aes(color = SNV_branch, size = SNV_branch))
    print(SNV_medicc_ggtree %<+% SNV_medicc_sub_info + geom_tippoint(aes(colour = LCM_region, shape = Type, size = Type)) + xlim(0,150) + 
            scale_color_manual(values = c(SNV_medicc_tree_colours[c(1,which(names(SNV_medicc_tree_colours) == sample))],"black", "black"), labels = c("diploid", sample, "yes")) +
            scale_size_manual(values = c(1,3,0.3,2), labels = c("10X", "LCM", "SNV", "CNA")) +
            scale_shape_manual(values = c(19,15)) +
            geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_sub_info$fill_col,  labels=SNV_medicc_sub_info$name) +
            guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2(),
          SNV_medicc_info[which(SNV_medicc_sub_info$name == sample), "node"])
    dev.off()
  }
}
  
if (tree == "GnT_LCM") {
  pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder.pdf"), width = 10, height = 14)
  SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(color = SNV_branch, size = SNV_branch)) %>% flip(12041, 6314)
  SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = region), size = .5) + 
    scale_color_manual(values = c(SNV_medicc_tree_colours,"black", "black"), labels = c("diploid", samples, "yes")) + scale_size_manual(values = c(0.3,2), labels = c("CNA", "SNV")) +
    geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_info$fill_col,  labels=SNV_medicc_info$name) +
    guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2()
  dev.off()
  
  #Hightlight LCM spots
  pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder_LCM.pdf"), width = 10, height = 14)
  SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(color = SNV_branch, size = SNV_branch)) %>% flip(12041, 6314)
  SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = LCM_region, shape = Type), size = .5) + 
    scale_color_manual(values = c(SNV_medicc_tree_colours[c(1,1,2,3,6)],"black", "black"), labels = c("diploid", samples, "yes")) + 
    scale_size_manual(values = c(0.3,2), labels = c("CNA", "SNV")) +
    scale_shape_manual(values = c(19,17,15)) +
    geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_info$fill_col,  labels=SNV_medicc_info$name) +
    guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2()
  dev.off()
  
  pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder_side.pdf"), width = 10, height = 14)
  SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(color = SNV_branch, size = SNV_branch)) %>% flip(12041, 6314)
  SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = Side, shape = Type), size = .5) + 
    scale_color_manual(values = c("grey", hcl(h = seq(15, 375, length = 3 + 1), l = 65, c = 100)[-4], "black", "black"), labels = c("", "Back", "Front", "Side")) +
    scale_size_manual(values = c(0.3,2), labels = c("CNA", "SNV")) +
    scale_shape_manual(values = c(19,17,15)) +
    geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_info$fill_col,  labels=SNV_medicc_info$name) +
    guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2()
  dev.off()
  
  #Hightlight G&T spots
  pdf(file = paste0(sub_folder, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder_G&T.pdf"), width = 10, height = 14)
  SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(color = SNV_branch, size = SNV_branch)) %>% flip(12041, 6314)
  SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = GnT_region, shape = Type), size = .5) + 
    scale_color_manual(values = c(SNV_medicc_tree_colours[c(1,1,2,4,5,6,7,1)],"black", "black"), labels = c("diploid", samples, "yes")) + 
    scale_size_manual(values = c(0.3,2), labels = c("CNA", "SNV")) +
    scale_shape_manual(values = c(19,17,15)) +
    geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_info$fill_col,  labels=SNV_medicc_info$name) +
    guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2()
  dev.off()
}




