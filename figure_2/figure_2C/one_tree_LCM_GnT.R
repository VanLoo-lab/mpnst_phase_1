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
library(tidytree)

### Define input, code and output dir
INPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_1/data/"
OUTPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_2/results/figure_2c/"
CODEDIR <- "~/Documents/GitHub/mpnst_phase_1/figure_2/figure_2C/"

#Load metadata
#source(file = "/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/code/10X_DLP_metadata.R")
#source(file = "/camp/project/proj-vanloo/analyses/hyan/mpnst/LCM/code/LCM_metadata.R")

source(paste0(CODEDIR, "metadata/10X_DLP_metadata_minimal.R"))
source(paste0(CODEDIR, "metadata/LCM_metadata.R"))

output.dir = OUTPUTDIR
setwd(output.dir)

####################################################################################################################################
### Part 1: Create tree group and some plots
####################################################################################################################################
scaling_ratio = 10
tree = "LCM" #GnT, LCM or GnT_LCM

raw_tree <- read_file(paste0(INPUTDIR, "LCM/22_kmeans_10X_DLP_filter_n_2/SNV_medicc_10X_",tree,"_tree_scale",scaling_ratio,".new"))
# for (b in unlist(region_odd_IDs)) {
#   print(b)
#   raw_tree <- gsub(b, paste0(samples[grep(b, region_odd_IDs)], "_", b), raw_tree)
# }
SNV_medicc_tree_add <- read.tree(text = raw_tree)

### Alex Stein
### Warning: "LCM_coordinates_tumour" from the LCM_metadata.R file is missing
###           I tried to build this file myself but not sure it is exactly the same

SNV_medicc_tree_metadata <- data.frame(label=SNV_medicc_tree_add$tip.label, region = gsub("_.*", "", SNV_medicc_tree_add$tip.label), stringsAsFactors = F) %>%
  left_join(LCM_coordinates_tumour %>% select(sample, Cell_Type, side) %>% mutate(sample = gsub("_.*", "", sample), side = gsub(".*-", "", side)), by = c("label" = "sample")) %>%
  mutate(region = ifelse(is.na(Cell_Type), region, Cell_Type)) %>% mutate(Cell_Type = ifelse(str_detect(label, "EGAF"), "G&T", ifelse(is.na(Cell_Type), "10X", "LCM"))) %>% rename(Type = Cell_Type, Side = side) %>%
  mutate(LCM_region = ifelse(Type == "LCM", region , Type)) %>% mutate(Side = ifelse(Type == "LCM", Side, ""))  #Add LCM region info 
SNV_medicc_tree_groups <- full_join(SNV_medicc_tree_add, SNV_medicc_tree_metadata, by = "label") 

SNV_medicc_tree_colours <- c("grey", region_colours[names(region_colours) %in% SNV_medicc_tree_groups@data[["region"]]])

SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups) + geom_tippoint(aes(shape = Type, colour = region), size = .5) + #geom_text(aes(label=node)) + #geom_tiplab(size = 8) +
  scale_color_manual(values = SNV_medicc_tree_colours, labels = c("diploid", samples)) +
  theme_tree2()

#Exploratory Plots
if (F) {
  pdf(file = paste0(OUTPUTDIR, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,".pdf"), width = 10, height = 14)
  print(SNV_medicc_ggtree)
  dev.off()
  
  pdf(file = paste0(OUTPUTDIR, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_lab_reorder.pdf"), width = 10, height = 56)
  (SNV_medicc_ggtree + geom_text(aes(label = node)))
  dev.off()
  
  pdf(file = paste0(OUTPUTDIR, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_reorder.pdf"), width = 10, height = 14)
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


####################################################################################################################################
### Some error fix

# SNV_medicc_info <- rbind(SNV_medicc_info, anti_join(data.frame(node = SNV_medicc_tree_groups@data$node,
#                                                                name = NA,
#                                                                fill_col = NA,
#                                                                SNV_branch = "SNV_F") ,SNV_medicc_info, by = "node")) #Add non trunk branches

# nodes from the ggtree object (recommended)
all_nodes <- SNV_medicc_ggtree$data$node   # or: tree_df <- SNV_medicc_ggtree$data

missing_df <- tibble::tibble(
  node = all_nodes,
  name = rep(NA_character_, length(all_nodes)),
  fill_col = rep(NA_character_, length(all_nodes)),
  SNV_branch = rep("SNV_F", length(all_nodes))
)

SNV_medicc_info <- dplyr::bind_rows(
  SNV_medicc_info,
  dplyr::anti_join(missing_df, SNV_medicc_info, by = "node")
)
####################################################################################################################################

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

# pdf(file = paste0(OUTPUTDIR, "/SNV_medicc_10X_",tree,"tree_ggtree_scale",scaling_ratio,"_nodelab_reorder.pdf"), width = 10, height = 14)
# SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = region), size = .5) + 
#   scale_color_manual(values = SNV_medicc_tree_colours, labels = c("diploid", samples)) +
#     geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_info$fill_col,  labels=SNV_medicc_info$name) +
#   theme_tree2()
# dev.off()


pdf(file = paste0(OUTPUTDIR, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder.pdf"), width = 10, height = 14)
SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(color = SNV_branch, size = SNV_branch)) %>% flip(11683, 6135)
SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = region), size = .5) + 
  scale_color_manual(values = c(SNV_medicc_tree_colours,"black", "black"), labels = c("diploid", samples, "yes")) + scale_size_manual(values = c(0.3,2), labels = c("CNA", "SNV")) +
  geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_info$fill_col,  labels=SNV_medicc_info$name) +
  guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2()
dev.off()

#Hightlight LCM spots
pdf(file = paste0(OUTPUTDIR, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder_LCM.pdf"), width = 10, height = 14)
SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(color = SNV_branch, size = SNV_branch)) %>% flip(11683, 6135)
SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = LCM_region, shape = Type), size = .5) + 
  scale_color_manual(values = c(SNV_medicc_tree_colours[c(1,2,3,6)],"black", "black"), labels = c("diploid", samples, "yes")) + 
  scale_size_manual(values = c(0.3,2), labels = c("CNA", "SNV")) +
  scale_shape_manual(values = c(19,15)) +
  geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_info$fill_col,  labels=SNV_medicc_info$name) +
  guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2()
dev.off()

pdf(file = paste0(OUTPUTDIR, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder_side.pdf"), width = 10, height = 14)
SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(color = SNV_branch, size = SNV_branch)) %>% flip(11683, 6135)
SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = Side, shape = Type), size = .5, show.legend = T) + 
  scale_color_manual(values = c("grey", hcl(h = seq(15, 375, length = 3 + 1), l = 65, c = 100)[-4], "black", "black"), labels = c("", "Back", "Front", "Side")) +
  scale_size_manual(values = c(0.3,2), labels = c("CNA", "SNV")) +
  scale_shape_manual(values = c(19,15)) +
  geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_info$fill_col,  labels=SNV_medicc_info$name) +
  guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2()
dev.off()

####################################################################################################################################
### Part 3: Create plots separate for P, R1 and R4
####################################################################################################################################

#Prune
# for (sample in c("P", "R1", "R4")) {
#   plot_height = round(length(unlist(Descendants(SNV_medicc_tree_add,  SNV_medicc_info[which(SNV_medicc_info$name == sample), "node"], type = "tips")))/length(SNV_medicc_tree_add[["tip.label"]])*35)
#   
#   SNV_medicc_ggtree_sub <- tree_subset(tree = SNV_medicc_tree_groups, node = sample, levels_back = 0)
#   if (sample == "P") {
#     sub_node_names <- c("P", "P_1", "P_2", "P_3")
#   } else {
#     if (sample == "R1") {
#       sub_node_names <- c("R1", "R1_2", "R1_3", "R1_4")
#     } else{
#       sub_node_names <- c("R4", "R4_1", "R4_2", "R4_3")
#     }
#   }
#   
#   SNV_medicc_sub_info <- data.frame(node = unlist(lapply(sub_node_names, function(n) {which(as_tibble(SNV_medicc_ggtree_sub)$label == n)})),
#                                     name = sub_node_names,
#                                     fill_col = node_colours[sub_node_names],
#                                     SNV_branch = "SNV_T", stringsAsFactors = F) %>% arrange(name) #Anndata for trunk branches
#   
#   ####################################################################################################################################
#   
#   ### Some error fix
#   
#   # SNV_medicc_sub_info <- rbind(SNV_medicc_sub_info, anti_join(data.frame(node = SNV_medicc_ggtree_sub@data$node,
#   #                                                                        name = NA,
#   #                                                                        fill_col = NA,
#   #                                                                        SNV_branch = "SNV_F") ,SNV_medicc_sub_info, by = "node")) #Add non trunk branches
#   # SNV_medicc_sub_info[SNV_medicc_sub_info$node %in% c(1),"SNV_branch"] <- "SNV_T" #Change root branch to SNV too #don't need this when pruned
#   
#   #get the node vector from the subset object safely
#   sub_nodes <- NULL
#   if (!is.null(SNV_medicc_ggtree_sub$data)) {                 # if it's a ggtree plot
#     sub_nodes <- SNV_medicc_ggtree_sub$data$node
#   } else if (!is.null(SNV_medicc_ggtree_sub@data$node)) {     # if it's treedata with @data
#     sub_nodes <- SNV_medicc_ggtree_sub@data$node
#   } else {
#     # fallback: fortify
#     sub_nodes <- ggtree::fortify(SNV_medicc_ggtree_sub)$node
#   }
#   #sub_nodes <- SNV_medicc_ggtree_sub@data$node
#   
#   missing_df <- tibble::tibble(
#     node = sub_nodes,
#     name = rep(NA_character_, length(sub_nodes)),
#     fill_col = rep(NA_character_, length(sub_nodes)),
#     SNV_branch = rep("SNV_F", length(sub_nodes))
#   )
#   
#   # SNV_medicc_sub_info <- dplyr::bind_rows(
#   #   SNV_medicc_sub_info,
#   #   dplyr::anti_join(missing_df, SNV_medicc_sub_info, by = "node")
#   # )
#   
#   sub_df <- as_tibble(SNV_medicc_ggtree_sub)
#   
#   SNV_medicc_sub_info <- tibble::tibble(
#     node = sapply(sub_node_names, function(n) sub_df$node[sub_df$label == n][1]),
#     name = sub_node_names,
#     fill_col = unname(node_colours[sub_node_names]),
#     SNV_branch = "SNV_T"
#   ) %>% dplyr::arrange(name)
# 
#   
#   
#   ####################################################################################################################################
#   
#   
#   pdf(file = paste0(OUTPUTDIR, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder_side_", sample, ".pdf"), width = 10, height = plot_height)
#   SNV_medicc_ggtree <- ggtree(SNV_medicc_ggtree_sub, aes(color = SNV_branch, size = SNV_branch))
#   print(SNV_medicc_ggtree %<+% SNV_medicc_sub_info + geom_tippoint(aes(colour = Side, shape = Type, size = Type)) + xlim(0,150) + 
#           scale_color_manual(values = c("grey", hcl(h = seq(15, 375, length = 3 + 1), l = 65, c = 100)[-4], "black", "black"), labels = c("", "Back", "Front", "Side")) +
#           scale_size_manual(values = c(1,3,0.3,2), labels = c("10X", "LCM", "SNV", "CNA")) +
#           scale_shape_manual(values = c(19,15)) +
#           geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_sub_info$fill_col,  labels=SNV_medicc_sub_info$name) +
#           guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2(), 
#         SNV_medicc_sub_info[which(SNV_medicc_sub_info$name == sample), "node"]) 
#   dev.off()
#   
#   pdf(file = paste0(OUTPUTDIR, "/SNV_medicc_10X_",tree,"_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder_LCM_", sample, ".pdf"), width = 10, height = plot_height)
#   SNV_medicc_ggtree <- ggtree(SNV_medicc_ggtree_sub, aes(color = SNV_branch, size = SNV_branch))
#   print(SNV_medicc_ggtree %<+% SNV_medicc_sub_info + geom_tippoint(aes(colour = LCM_region, shape = Type, size = Type)) + xlim(0,150) + 
#           scale_color_manual(values = c(SNV_medicc_tree_colours[c(1,which(names(SNV_medicc_tree_colours) == sample))],"black", "black"), labels = c("diploid", sample, "yes")) +
#           scale_size_manual(values = c(1,3,0.3,2), labels = c("10X", "LCM", "SNV", "CNA")) +
#           scale_shape_manual(values = c(19,15)) +
#           geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=SNV_medicc_sub_info$fill_col,  labels=SNV_medicc_sub_info$name) +
#           guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2(),
#         SNV_medicc_info[which(SNV_medicc_sub_info$name == sample), "node"])
#   dev.off()
# }



####################################################################################################################################
### Part 3: Create plots separate for P, R1 and R4 -- Version 2
####################################################################################################################################

# Prune
for (sample in c("P", "R1", "R4")) {
  
  plot_height <- round(
    length(unlist(Descendants(SNV_medicc_tree_add,
                              SNV_medicc_info[SNV_medicc_info$name == sample, "node"],
                              type = "tips"))) /
      length(SNV_medicc_tree_add$tip.label) * 35
  )
  
  sub_tree <- tidytree::tree_subset(tree = SNV_medicc_tree_groups, node = sample, levels_back = 0)
  #sub_df   <- tibble::as_tibble(ggtree::ggtree(sub_tree))  # ensures we have node/label
  
  p_sub  <- ggtree::ggtree(sub_tree)
  sub_df <- p_sub$data

  sub_node_names <- if (sample == "P") {
    c("P","P_1","P_2","P_3")
  } else if (sample == "R1") {
    c("R1","R1_2","R1_3","R1_4")
  } else {
    c("R4","R4_1","R4_2","R4_3")
  }
  
  # trunk nodes (SNV_T)
  SNV_medicc_sub_info <- tibble::tibble(
    node = sapply(sub_node_names, function(n) sub_df$node[sub_df$label == n][1]),
    name = sub_node_names,
    fill_col = unname(node_colours[sub_node_names]),
    SNV_branch = "SNV_T"
  ) %>% dplyr::arrange(name)
  
  # add all other nodes as SNV_F (so SNV_branch exists everywhere)
  all_nodes <- sub_df$node
  missing_df <- tibble::tibble(
    node = all_nodes,
    name = NA_character_,
    fill_col = NA_character_,
    SNV_branch = "SNV_F"
  )
  SNV_medicc_sub_info <- dplyr::bind_rows(
    SNV_medicc_sub_info,
    dplyr::anti_join(missing_df, SNV_medicc_sub_info, by = "node")
  )
  
  # ---- Plot 1 (tips coloured by Side) ----
  pdf(paste0(OUTPUTDIR, "/SNV_medicc_10X_", tree, "_tree_ggtree_scale", scaling_ratio,
             "_sized_nodelab_reorder_side_", sample, ".pdf"),
      width = 10, height = plot_height)
  
  p <- ggtree::ggtree(sub_tree, aes(size = SNV_branch), color = "black") #%>%
    #flip(SNV_medicc_sub_info[SNV_medicc_sub_info$name == sample, "node"])
  
  print(
    (p %<+% SNV_medicc_sub_info) +
      geom_tippoint(aes(colour = Side, shape = Type, size = Type)) +
      xlim(0,150) +
      scale_color_manual(values = c("Back"="#1b9e77","Front"="#d95f02","Side"="#7570b3"), na.value="grey50") +
      scale_size_manual(values = c("10X"=1, "LCM"=3, "SNV_F"=0.3, "SNV_T"=2)) +
      scale_shape_manual(values = c("10X"=19, "LCM"=15)) +
      geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) +
      scale_fill_manual(values = setNames(SNV_medicc_sub_info$fill_col, SNV_medicc_sub_info$name)) +
      guides(colour="none", size = guide_legend(title="Branch Type"), fill="none") +
      theme_tree2()
  )
  dev.off()
  
  # ---- Plot 2 (tips coloured by LCM_region) ----
  pdf(paste0(OUTPUTDIR, "/SNV_medicc_10X_", tree, "_tree_ggtree_scale", scaling_ratio,
             "_sized_nodelab_reorder_LCM_", sample, ".pdf"),
      width = 10, height = plot_height)
  
  p2 <- ggtree::ggtree(sub_tree, aes(size = SNV_branch), color = "black") #%>%
    #flip(SNV_medicc_sub_info[SNV_medicc_sub_info$name == sample, "node"])
  
  print(
    (p2 %<+% SNV_medicc_sub_info) +
      geom_tippoint(aes(colour = LCM_region, shape = Type, size = Type)) +
      xlim(0,150) +
      scale_color_manual(values = SNV_medicc_tree_colours, na.value="grey50") +
      scale_size_manual(values = c("10X"=1, "LCM"=3, "SNV_F"=0.3, "SNV_T"=2)) +
      scale_shape_manual(values = c("10X"=19, "LCM"=15)) +
      geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) +
      scale_fill_manual(values = setNames(SNV_medicc_sub_info$fill_col, SNV_medicc_sub_info$name)) +
      guides(colour="none", size = guide_legend(title="Branch Type"), fill="none") +
      theme_tree2()
  )
  dev.off()
}


