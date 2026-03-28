### title: Plot singl-cell resolved phylogeny combining SNV and CN data

# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"

# Define input directory with data and output directory to save the figure
input.dir <- paste0(zenodo.dir, "data/")
output.dir <- paste0(zenodo.dir, "results/figure_1/")
setwd(output.dir)

####################################################################################################################################
### Part 0: Load libraries and data
####################################################################################################################################

### Load libraries
library(tidyverse)
library(ape)
library(phangorn)
library(treeio)
library(ggtree)
library(gridExtra)

### Load number of SNVs from bulk tree
bulk_SNV_consensus_clusters <- readRDS(paste0(input.dir,"bulk/DPClust/MPNST_1_consensus_clusters_40_mod_clean.rds"))

### Load number SNVs from subclone tree
filter_n_cells <- 2 #SNVs with less than or equal to filter_n_cells cell removed
CN_based_SNV_clusters_ids <- readRDS(paste0(input.dir,"bulk/SNV_analysis/MPNST_all_sub_clonal_SNV_G1000filtered_clusters_ex_WGD_N_ids_32.rds"))
CN_based_SNV_clusters_sizes <- unlist(lapply(CN_based_SNV_clusters_ids, length))

### Load k-means cluster info for scDNA-seq
K <- 22
scDNA_CN_cluster_ids <- readRDS(paste0(input.dir, "scDNA/CN_profiles/MPNST_all_k_means_K",K,"_clusters.rds"))

### Load MEDICC2 trees after loading metadata
# --> see medicc.trees below

### Define scaling ratio
scaling_ratio = 10 # Decide scaling ratio (10 SNVs = 1 CNA)


####################################################################################################################################
### Part 1: Sort metadata from scDNA-seq
####################################################################################################################################

### Some sample names and color definitions

samples = c("R1", "R2", "R3", "R4", "R5", "P")
names(samples) = c("FIT208A3", "FIT208A4", "FIT208A5", "FIT208A6", "FIT208A7", "FIT208A8")
# SUBSETS = c("P", "R1", "R2", "R3", "R4", "R5", "G2MS", "Ribo")

bulk_samples = c("VER236A1", "VER236A2", "VER236A3", "VER236A4", "VER236A5", "VER236A6")
names(bulk_samples) <- samples

#Colours
region_colours <- c("#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3", "#F8766D")[c(6,1:5)]
names(region_colours) <- samples[c(6,1:5)]

### Cluster info generated after running kmeans in ASCAT.sc.R then cleaned up

K = 22
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

names(scDNA_CN_cluster_ids) <- ave(scDNA_CN_clusters_origins, scDNA_CN_clusters_origins, FUN = function(i) paste0(i, "_", seq_along(i)))

####################################################################################################################################
### Part 2: Load medicc CN cluster trees
####################################################################################################################################

### Load MEDICC2 trees
medicc.dir = paste0(input.dir, "scDNA/MEDICC2/")
run <- "22_kmeans_10X_DLP"
medicc.newick.trees <- paste0(medicc.dir, "scDNA_", run, "_2.5_", names(scDNA_CN_cluster_ids), "/scDNA_", run, "_2.5_", names(scDNA_CN_cluster_ids), "_final_tree.new")
medicc.trees <- lapply(medicc.newick.trees, function(t) {
  read.tree(t)
})
names(medicc.trees) <- names(scDNA_CN_cluster_ids)

### Prune tree (check how to prune below)
medicc.trees.pruned <- lapply(1:length(medicc.trees), function(t) {
  node <- medicc.trees[[t]][["edge"]][1,2] #First two is root and second column is int node with all descendants except diploid (last row)
  pruned <- extract.clade(medicc.trees[[t]], node)
  pruned$node.label[1] <- names(medicc.trees)[t] #Rename root node as subclone
  return(pruned)
})
names(medicc.trees.pruned) <- names(medicc.trees)


#Load number of mutations for branches in tree
bulk_SNVs <- bulk_SNV_consensus_clusters %>% group_by(cluster) %>% summarise(SNVs = n())
bulk_SNVs$SNVs <- bulk_SNVs$SNVs/scaling_ratio


#Reorder cluster
if (run == "22_kmeans_10X_DLP") {
  J = 32
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


####################################################################################################################################
### Part 3: Make trunk and stem
####################################################################################################################################

### Make trunk tree with bulk SNV numbers as branch lengths
SNV.trunk.newick.tree <- paste0("(Fertilised_egg:1,
(P:",bulk_SNVs[2,2],",
(R1:",bulk_SNVs[14,2],",R5:",bulk_SNVs[4,2],"):",bulk_SNVs[15,2],",
((R2:",bulk_SNVs[10,2],",R4:",bulk_SNVs[6,2],"):",bulk_SNVs[11,2],",R3:",bulk_SNVs[8,2],"):",bulk_SNVs[12,2],"):",bulk_SNVs[16,2],");")

SNV.trunk.tree <- read.tree(text = SNV.trunk.newick.tree)

### Use subclone SNV sizes
subclone_sizes <- CN_based_SNV_clusters_sizes/scaling_ratio
missing_sizes <- c(3, 1, 6)/scaling_ratio
names(missing_sizes) <- c("P_1", "P_3", "R5_1")

SNV.newick.tree <- paste0("(Fertilised_egg:1, 
((P_1:",missing_sizes["P_1"],",P_2:", subclone_sizes["P_2"],",P_3:",missing_sizes["P_3"],")P:",bulk_SNVs[bulk_SNVs$cluster=="00000X",2],",
((R1_2:",subclone_sizes["R1_2"],",R1_3:",subclone_sizes["R1_3"],",R1_4:",subclone_sizes["R1_4"],")R1:",bulk_SNVs[bulk_SNVs$cluster=="X00000",2],",
(R5_1:",missing_sizes["R5_1"],",R5_2:",subclone_sizes["R5_2"],",R5_3:",subclone_sizes["R5_3"],",R5_4:",subclone_sizes["R5_4"],")R5:",bulk_SNVs[bulk_SNVs$cluster=="0000X0",2],")R1_R5:",bulk_SNVs[bulk_SNVs$cluster=="X000X0",2],",
(((R2_1:",subclone_sizes["R2_1"],",R2_2:",subclone_sizes["R2_2"],",R2_3:",subclone_sizes["R2_3"],")R2:",bulk_SNVs[bulk_SNVs$cluster=="0X0000",2],",(R4_1:",subclone_sizes["R4_1"],",R4_2:",subclone_sizes["R4_2"],",R4_3:",subclone_sizes["R4_3"],")R4:",bulk_SNVs[bulk_SNVs$cluster=="000X00",2],")R2_R4:",bulk_SNVs[bulk_SNVs$cluster=="0X0X00",2],",R3_1:",bulk_SNVs[bulk_SNVs$cluster=="00X000",2],")R2_R3_R4:",bulk_SNVs[bulk_SNVs$cluster=="0XXX00",2],
                          ",R1_1:",subclone_sizes["R1_1"],")MRCA:",bulk_SNVs[bulk_SNVs$cluster=="XXXXXX",2],");")

SNV.tree <- read.tree(text = SNV.newick.tree)

####################################################################################################################################
### Part 4: Add medicc trees
####################################################################################################################################

SNV_medicc_tree <- SNV.tree

for (i in 1:length(medicc.trees.pruned)) {
  where = which(SNV_medicc_tree[["tip.label"]] == names(medicc.trees.pruned)[i])
  SNV_medicc_tree <- bind.tree(SNV_medicc_tree, medicc.trees.pruned[[i]], where = where, position = 0, interactive = FALSE)
}

### Save the tree
write.tree(SNV_medicc_tree, paste0("SNV_medicc_tree_scale",scaling_ratio,".new"))

####################################################################################################################################
### Part 5: Make a pretty plot
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

SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups)

pdf(file = "figure_1F.pdf", width = 10, height = 14)
SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups, aes(size = SNV_branch), color="black") %>% flip(11499, 6043)
SNV_medicc_ggtree %<+% SNV_medicc_info + geom_tippoint(aes(colour = region), size = .5, show.legend = F) + 
  scale_color_manual(values = c(SNV_medicc_tree_colours,"black", "black"), labels = c("diploid", samples, "yes")) + scale_size_manual(values = c(0.3,2), labels = c("CNA", "SNV")) +
  geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) + scale_fill_manual(values=node_colours,  labels=names(node_colours)) +
  guides(color = guide_legend(title = "Region"), size = guide_legend(title = "Branch Type"), fill = "none") + theme_tree2()
  # guides(color = "none", size = guide_legend(title = "Branch Type"), fill = guide_legend(title = "Subclone", override.aes = aes(label = ""))) + theme_tree2()
dev.off()


