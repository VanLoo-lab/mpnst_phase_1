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
OUTPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_1/results/figure_1f/"
CODEDIR <- "~/Documents/GitHub/MPNST-Phase-1/figure_1/figure_1F/"

### Load metadata
source(paste0(CODEDIR,"metadata/10X_DLP_metadata_minimal.R"))

### Set output directory
output.dir <- OUTPUTDIR
setwd(output.dir)

####################################################################################################################################
### Part 1: Choose which clustering run to use
####################################################################################################################################
filter_n_cells <- 2 #SNVs with less than or equal to filter_n_cells cell removed 
sub_folder <- paste0(run, "_filter_n_", filter_n_cells)
system(paste0("mkdir -p ", output.dir, sub_folder))

####################################################################################################################################
### Part 2: Load medicc CN cluster trees
####################################################################################################################################
medicc.newick.trees <- paste0(medicc.dir, "scDNA_", run, "_2.5", "/scDNA_", run, "_2.5_", names(scDNA_CN_cluster_ids), "/scDNA_", run, "_2.5_", names(scDNA_CN_cluster_ids), "_final_tree.new")
medicc.trees <- lapply(medicc.newick.trees, function(t) {
  read.tree(t)
})
names(medicc.trees) <- names(scDNA_CN_cluster_ids)
#Prune tree (check how to prune below)
# plot(medicc.trees[["R4_2"]])
# plot(extract.clade(medicc.trees[["R4_2"]], medicc.trees[["R4_2"]][["edge"]][1,2]))
medicc.trees.pruned <- lapply(1:length(medicc.trees), function(t) {
  node <- medicc.trees[[t]][["edge"]][1,2] #First two is root and second column is int node with all descendants except diploid (last row)
  pruned <- extract.clade(medicc.trees[[t]], node)
  pruned$node.label[1] <- names(medicc.trees)[t] #Rename root node as subclone
  return(pruned)
})
names(medicc.trees.pruned) <- names(medicc.trees)

#Decide scaling ratio (2 SNVs = 1 CNA)
scaling_ratio = 10

#Load number of mutations for branches in tree
if (T) {
  #Load no. SNVs from bulk tree
  bulk_SNV_consensus_clusters <- readRDS(paste0(INPUTDIR,"bulk/DPClust/CICC_ex1_ex2/MPNST_1_consensus_clusters_40_mod_clean.rds"))
  bulk_SNVs <- bulk_SNV_consensus_clusters %>% group_by(cluster) %>% summarise(SNVs = n())
  bulk_SNVs$SNVs <- bulk_SNVs$SNVs/scaling_ratio
  
  #Load no. SNVs from subclone tree
  upper_limit = 2
  J = 32
  CN_based_SNV_clusters_ids <- readRDS(paste0(INPUTDIR,"10x_DLP/snv_mnv_indel_tumouronly/SNV_analysis/",run,"_filter_n_", filter_n_cells, "/MPNST_all_sub_clonal_SNV_G1000filtered_clusters_ex_WGD_N_ids_",J,".rds"))
  CN_based_SNV_clusters_sizes <- unlist(lapply(CN_based_SNV_clusters_ids, length))
  
  #Reorder cluster
  if (run == "22_kmeans_10X_DLP") {
    SNV_clusters <- 1:J
    names(SNV_clusters) <- c("R4_2", "R4_3", "R4_2+R4_3", "R4", "R5_4",
                             "R5_2", "Clonal_A", "Clonal_B", "Clonal_C", "R5_2+R5_3",
                             "R5_3", "R1_4", "R1+R5_A", "R1+R5_B", "R1+R5_C",
                             "R5_1+R5_2+R5_3", "R5", "R5_2+R5_3+R5_4", "R5_2+R5_4", "R2_2",
                             "R2_3", "R2", "R2+R3+R4", "R2_1+R2_2", "R2_1",
                             "R3", "R4_1", "R1_1", "P", "P_2",
                             "R1_2", "R1_3")
    names(CN_based_SNV_clusters_sizes) <- names(SNV_clusters)
    
    SNV_cluster_reorder <- c(7,8,9,#Clonal
                             29,30,#P
                             28,#R1extra
                             13,14,15,#R1/R5
                             32,31,12,#R1
                             17,18,16,19,10,6,5,11,#R5
                             23,#R2,3,4
                             22,24,20,25,21,#R2
                             26,#R3
                             4,3,1,2,27#R4
    )
  }
  CN_based_SNV_clusters_sizes_reorder <- CN_based_SNV_clusters_sizes[SNV_cluster_reorder]
}

####################################################################################################################################
### Part 3: Make trunk and stem
####################################################################################################################################
SNV.trunk.newick.tree <- paste0("(Fertilised_egg:1,
(P:",bulk_SNVs[2,2],",
(R1:",bulk_SNVs[14,2],",R5:",bulk_SNVs[4,2],"):",bulk_SNVs[15,2],",
((R2:",bulk_SNVs[10,2],",R4:",bulk_SNVs[6,2],"):",bulk_SNVs[11,2],",R3:",bulk_SNVs[8,2],"):",bulk_SNVs[12,2],"):",bulk_SNVs[16,2],");")

SNV.trunk.tree <- read.tree(text = SNV.trunk.newick.tree)
# plot(SNV.trunk.tree)

#Use subclone SNV sizes
if (T) {
  if (F) {
    #To get SNVs from cluster that's too small
    scDNA_CN_clusters_SNVs_resize <- readRDS(paste0("../snv_mnv_indel_tumouronly/SNV_analysis/", run, "_filter_n_", filter_n_cells,"/MPNST_all_sub_clonal_SNV_G1000filtered_clusters_ex_WGD_N_mtx.rds"))
    # scDNA_CN_clusters_SNVs_resize[scDNA_CN_clusters_SNVs_resize[,"P_1"] != 0,"P_1"]
    #Find SNVs which aren't present in any cluster except cluster of interest
    unique_SNVs <- unlist(lapply(1:nrow(scDNA_CN_clusters_SNVs_resize), function(v) {
      all(scDNA_CN_clusters_SNVs_resize[v,-which(colnames(scDNA_CN_clusters_SNVs_resize) == "R5_1")] == 0)
    }))
    sum(unique_SNVs)
  } #check manually
  
  subclone_sizes <- CN_based_SNV_clusters_sizes/scaling_ratio
  # missing_sizes <- c(21, 4, 6, 7, 12)/scaling_ratio
  # names(missing_sizes) <- c("P_1", "P_2", "R1_3", "R5_1", "R1_5") #4.1.2
  missing_sizes <- c(3, 1, 6)/scaling_ratio
  names(missing_sizes) <- c("P_1", "P_3", "R5_1")
  
  SNV.newick.tree <- paste0("(Fertilised_egg:1,
((P_1:",missing_sizes["P_1"],",P_2:", subclone_sizes["P_2"],",P_3:",missing_sizes["P_3"],")P:",bulk_SNVs[bulk_SNVs$cluster=="00000X",2],",
((R1_2:",subclone_sizes["R1_2"],",R1_3:",subclone_sizes["R1_3"],",R1_4:",subclone_sizes["R1_4"],")R1:",bulk_SNVs[bulk_SNVs$cluster=="X00000",2],",
(R5_1:",missing_sizes["R5_1"],",R5_2:",subclone_sizes["R5_2"],",R5_3:",subclone_sizes["R5_3"],",R5_4:",subclone_sizes["R5_4"],")R5:",bulk_SNVs[bulk_SNVs$cluster=="0000X0",2],")R1_R5:",bulk_SNVs[bulk_SNVs$cluster=="X000X0",2],",
(((R2_1:",subclone_sizes["R2_1"],",R2_2:",subclone_sizes["R2_2"],",R2_3:",subclone_sizes["R2_3"],")R2:",bulk_SNVs[bulk_SNVs$cluster=="0X0000",2],",(R4_1:",subclone_sizes["R4_1"],",R4_2:",subclone_sizes["R4_2"],",R4_3:",subclone_sizes["R4_3"],")R4:",bulk_SNVs[bulk_SNVs$cluster=="000X00",2],")R2_R4:",bulk_SNVs[bulk_SNVs$cluster=="0X0X00",2],",R3_1:",bulk_SNVs[bulk_SNVs$cluster=="00X000",2],")R2_R3_R4:",bulk_SNVs[bulk_SNVs$cluster=="0XXX00",2],
",R1_1:",subclone_sizes["R1_1"],")MRCA:",bulk_SNVs[bulk_SNVs$cluster=="XXXXXX",2],");")
} else {
  SNV.newick.tree <- paste0("(Fertilised_egg:1,
((P_1:10,P_2:10,P_3:10,P_4:10)P:",bulk_SNVs[13,2],",
((R1_1:10,R1_2:10,R1_3:10,R1_6:10)R1:",bulk_SNVs[5,2],",
((R5_2:10,R5_3:10)R5_2_R5_3:10,R5_1:7)R5:",bulk_SNVs[8,2],")R1_R5:",bulk_SNVs[12,2],",
(((R2_1:10,R2_2:10,R2_3:10)R2:",bulk_SNVs[4,2],",(R4_1:10,R4_2:10)R4:",bulk_SNVs[3,2],")R2_R4:",bulk_SNVs[15,2],",R3_1:",bulk_SNVs[9,2],")R2_R3_R4:",bulk_SNVs[14,2],",
(R1_4:10,R1_5:10)R1e:10)MRCA:",bulk_SNVs[2,2],");") #old ignore
}

SNV.tree <- read.tree(text = SNV.newick.tree)
# ggtree(SNV.tree) + geom_tiplab() + geom_nodelab()

####################################################################################################################################
### Part 4: Add medicc trees
####################################################################################################################################
SNV_medicc_tree <- SNV.tree
for (i in 1:length(medicc.trees.pruned)) {
  where = which(SNV_medicc_tree[["tip.label"]] == names(medicc.trees.pruned)[i])
  SNV_medicc_tree <- bind.tree(SNV_medicc_tree, medicc.trees.pruned[[i]], where = where, position = 0, interactive = FALSE)
}
#Manually rename
# plot(SNV_medicc_tree, show.tip.label = F)
write.tree(SNV_medicc_tree, paste0(sub_folder, "/SNV_medicc_tree_scale",scaling_ratio,".new"))

####################################################################################################################################
### Part 5: Make prettier plot
####################################################################################################################################
SNV_medicc_tree_groups <- full_join(SNV_medicc_tree, data.frame(label=SNV_medicc_tree$tip.label, region = gsub("_.*", "", SNV_medicc_tree$tip.label), stringsAsFactors = F), by = "label") 

SNV_medicc_tree_colours <- c("grey", region_colours[names(region_colours) %in% SNV_medicc_tree_groups@data[["region"]]])

SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups) + geom_tippoint(aes(colour = region), size = .5) + #geom_text(aes(label=node)) + #geom_tiplab(size = 8) +
  scale_color_manual(values = SNV_medicc_tree_colours, labels = c("diploid", samples)) +
  theme_tree2()

# Make a ggtree plot and attach tip metadata properly
tip_df <- tibble::tibble(
  label  = SNV_medicc_tree$tip.label,
  region = gsub("_.*", "", SNV_medicc_tree$tip.label)
)

SNV_medicc_tree_colours <- c("grey", region_colours[names(region_colours) %in% tip_df$region])

SNV_medicc_ggtree <- ggtree(SNV_medicc_tree) %<+% tip_df +
  geom_tippoint(aes(colour = region), size = .5) +
  scale_color_manual(values = SNV_medicc_tree_colours, labels = c("diploid", samples)) +
  theme_tree2()



#Plots
if (T) {
  pdf(file = paste0(sub_folder, "/SNV_medicc_tree_ggtree_scale",scaling_ratio,".pdf"), width = 10, height = 14)
  print(SNV_medicc_ggtree)
  dev.off()
  
  pdf(file = paste0(sub_folder, "/SNV_medicc_tree_ggtree_scale",scaling_ratio,"_lab_reorder.pdf"), width = 10, height = 56)
  (SNV_medicc_ggtree + geom_text(aes(label = node))) %>% flip(11499, 6043) 
  # SNV_medicc_ggtree %>% flip(6042, 5000) %>% flip(6042, 3651)
  dev.off()
  
  pdf(file = paste0(sub_folder, "/SNV_medicc_tree_ggtree_scale",scaling_ratio,"_reorder.pdf"), width = 10, height = 14)
  SNV_medicc_ggtree %>% flip(11499, 6043) 
  # SNV_medicc_ggtree %>% flip(6042, 5000) %>% flip(6042, 3651)
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
                              SNV_branch = "SNV_F"),SNV_medicc_info, by = "node")) #Add non trunk branches

tree_df <- SNV_medicc_ggtree$data
nodes <- tree_df$node

all_nodes <- SNV_medicc_ggtree$data$node

missing_df <- tibble::tibble(
  node = all_nodes,
  name = NA_character_,
  fill_col = NA_character_,
  SNV_branch = "SNV_F"
)

SNV_medicc_info <- dplyr::bind_rows(
  SNV_medicc_info,
  dplyr::anti_join(missing_df, SNV_medicc_info, by = "node")
)


# SNV_medicc_info[SNV_medicc_info$node %in% c(1,4208),"SNV_branch"] <- "SNV_T" #Change root branch to SNV too
SNV_medicc_info[SNV_medicc_info$node %in% c(1),"SNV_branch"] <- "SNV_T" #Change root branch to SNV too

# SNV_medicc_info <- data.frame(node = c(3262, 3263, 3456, 3264, 3370, 3431,
#                                        3651),
#                               name = c("MRCA", "P", "P_4", "P_1", "P_2", "P_3",
#                                        "R1_R5"),
#                               fill_col = c("white", rep(cluster_colours_reorder[10],5),
#                                            cluster_colours_reorder[7]), stringsAsFactors = F)
SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups)

pdf(file = paste0(sub_folder, "/SNV_medicc_tree_ggtree_scale",scaling_ratio,"_nodelab_reorder.pdf"), width = 10, height = 14)
SNV_medicc_ggtree %>% flip(11499, 6043) %<+% SNV_medicc_info + geom_tippoint(aes(colour = region), size = .5) + 
  scale_color_manual(values = SNV_medicc_tree_colours, labels = c("diploid", samples)) +
    geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=node_colours, labels=names(node_colours)) +
  theme_tree2()
dev.off()

pdf(file = paste0(sub_folder, "/SNV_medicc_tree_ggtree_scale",scaling_ratio,"_sized_nodelab_reorder.pdf"), width = 10, height = 14)
#SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(color = SNV_branch, size = SNV_branch)) %>% flip(11499, 6043)
SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(size = SNV_branch), color="black") %>% flip(11499, 6043)
SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = region), size = .5, show.legend = F) + 
  scale_color_manual(values = c(SNV_medicc_tree_colours,"black", "black"), labels = c("diploid", samples, "yes")) + scale_size_manual(values = c(0.3,2), labels = c("CNA", "SNV")) +
  geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=node_colours,  labels=names(node_colours)) +
  guides(color = guide_legend(title = "Region"), size = guide_legend(title = "Branch Type"), fill = "none") + theme_tree2()
  # guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2()
dev.off()


