### title: Plot scDNA MEDICC2 tree with mapped LCM spots

# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"

# Define input and output directories requested by the user
input.dir <- paste0(zenodo.dir, "data/")
output.dir <- paste0(zenodo.dir, "results/figure_2/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
setwd(output.dir)

####################################################################################################################################
### Part 0: Load libraries
####################################################################################################################################

library(tidyverse)
library(ape)
library(phangorn)
library(treeio)
library(ggtree)
library(gridExtra)
library(tidytree)

####################################################################################################################################
### Part 1: Inline metadata (minimal subset copied from metadata scripts)
####################################################################################################################################

### Metadata from 10X_DLP_metadata_minimal.R (only required objects)
samples <- c("R1", "R2", "R3", "R4", "R5", "P")
names(samples) <- c("FIT208A3", "FIT208A4", "FIT208A5", "FIT208A6", "FIT208A7", "FIT208A8")

region_colours <- c("#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3", "#F8766D")[c(6, 1:5)]
names(region_colours) <- samples[c(6, 1:5)]

### Metadata from LCM_metadata.R (only required objects)
LCM.dir <- paste0(input.dir, "LCM/")

LCM_locations <- read.delim(paste0(LCM.dir, "LCM_data_final_manual.csv"), sep = ",")

LCM_barcodes <- read.delim(paste0(LCM.dir, "barcodes.tsv")) %>%
	mutate(region = ifelse(str_detect(sample, "006"), "Primary", ifelse(str_detect(sample, "005"), "T1.1", "T4.1")))

both_haplo_count_combined_by_cell <- readRDS(paste0(LCM.dir, "Genotype_SNP/both_haplo_count_combined_by_cell.rds")) %>%
	mutate(Haplotype_Ratio = Haplo_1_Count / Haplo_2_Count)

LCM_asCN_profile_samples <- data.frame(
	sample = unique(both_haplo_count_combined_by_cell$Barcode),
	CN_profile = "Yes",
	stringsAsFactors = FALSE
)

LCM_coordinates <- data.frame(
	barcode = gsub("\\)", "", gsub(".*\\(", "", LCM_locations$UDF.Index)),
	coord.x = as.numeric(gsub(",.*", "", LCM_locations$Coordinates..X.Y..mm)),
	coord.y = as.numeric(gsub(".*,", "", LCM_locations$Coordinates..X.Y..mm)),
	region = LCM_locations$Region,
	side = LCM_locations$Side,
	sample_number = LCM_locations$Sample.Name.LCM.area
) %>%
	left_join(LCM_barcodes, by = c("barcode", "region")) %>%
	left_join(LCM_asCN_profile_samples, by = "sample") %>%
	mutate(CN_profile = ifelse(is.na(CN_profile), "No", "Yes"))

LCM_coordinates_tumour <- LCM_coordinates %>%
	left_join(
		both_haplo_count_combined_by_cell,
		by = c("barcode", "coord.x", "coord.y", "region", "side", "sample_number", "CN_profile", "sample" = "Barcode")
	) %>%
	filter(Haplotype_Ratio > 2) %>%
	filter(sample != "VER315A36_S270_L006") %>%
	filter(sample != "VER315A110_S319_L007")

####################################################################################################################################
### Part 2: Build the LCM tree and annotate nodes
####################################################################################################################################

scaling_ratio <- 10
tree <- "LCM"

tree_file <- paste0(input.dir, "LCM/22_kmeans_10X_DLP_filter_n_2/SNV_medicc_10X_", tree, "_tree_scale", scaling_ratio, ".new")
SNV_medicc_tree_add <- read.tree(file = tree_file)

SNV_medicc_tree_metadata <- data.frame(
	label = SNV_medicc_tree_add$tip.label,
	region = gsub("_.*", "", SNV_medicc_tree_add$tip.label),
	stringsAsFactors = FALSE
) %>%
	left_join(
		LCM_coordinates_tumour %>%
			select(sample, Cell_Type, side) %>%
			mutate(sample = gsub("_.*", "", sample), side = gsub(".*-", "", side)),
		by = c("label" = "sample")
	) %>%
	mutate(region = ifelse(is.na(Cell_Type), region, Cell_Type)) %>%
	mutate(Cell_Type = ifelse(str_detect(label, "EGAF"), "G&T", ifelse(is.na(Cell_Type), "10X", "LCM"))) %>%
	rename(Type = Cell_Type, Side = side) %>%
	mutate(LCM_region = ifelse(Type == "LCM", region, Type)) %>%
	mutate(Side = ifelse(Type == "LCM", Side, ""))

SNV_medicc_tree_groups <- full_join(treeio::as.treedata(SNV_medicc_tree_add), SNV_medicc_tree_metadata, by = "label")

cluster_colours_reorder <- c("black", "darkslategray4", "darkgreen", "#00BA38", "#619CFF", "#00BFC4", "coral", "#B79F00", "#F564E3", "#F8766D")
node_names <- c("MRCA", "P", "P_1", "P_2", "P_3",
								"R1", "R1_1", "R1_2", "R1_3", "R1_4", "R1_R5",
								"R2", "R2_1", "R2_2", "R2_3", "R2_R3_R4", "R2_R4",
								"R3_1", "R4", "R4_1", "R4_2", "R4_3",
								"R5", "R5_1", "R5_2", "R5_3", "R5_4")

node_colours <- c("white", rep(cluster_colours_reorder[10], 4),
									rep(cluster_colours_reorder[8], 5), cluster_colours_reorder[7],
									rep(cluster_colours_reorder[4], 4), cluster_colours_reorder[2], cluster_colours_reorder[3],
									cluster_colours_reorder[6], rep(cluster_colours_reorder[5], 4), rep(cluster_colours_reorder[9], 5))
names(node_colours) <- node_names

SNV_medicc_info <- data.frame(
	node = unlist(lapply(node_names, function(n) which(as_tibble(SNV_medicc_tree_groups)$label == n))),
	name = node_names,
	fill_col = node_colours,
	SNV_branch = "SNV_T",
	stringsAsFactors = FALSE
) %>%
	arrange(name)

SNV_medicc_ggtree <- ggtree(SNV_medicc_tree_groups)
missing_df <- tibble(
	node = SNV_medicc_ggtree$data$node,
	name = NA_character_,
	fill_col = NA_character_,
	SNV_branch = "SNV_F"
)

SNV_medicc_info <- bind_rows(
	SNV_medicc_info,
	anti_join(missing_df, SNV_medicc_info, by = "node")
)
SNV_medicc_info[SNV_medicc_info$node %in% 1, "SNV_branch"] <- "SNV_T"

tip_colours_lcm <- c(
	"10X" = "grey70",
	"P" = unname(region_colours["P"]),
	"R1" = unname(region_colours["R1"]),
	"R4" = unname(region_colours["R4"]),
	"LCM" = "black",
	"G&T" = "black"
)

####################################################################################################################################
### Part 3: Generate only the three requested LCM subtree figures
####################################################################################################################################

for (sample in c("R4", "P", "R1")) {
	plot_height <- round(
		length(unlist(Descendants(
			SNV_medicc_tree_add,
			SNV_medicc_info[SNV_medicc_info$name == sample, "node"],
			type = "tips"
		))) / length(SNV_medicc_tree_add$tip.label) * 35
	)

	sub_tree <- tidytree::tree_subset(tree = SNV_medicc_tree_groups, node = sample, levels_back = 0)
	sub_df <- ggtree::ggtree(sub_tree)$data

	sub_node_names <- if (sample == "P") {
		c("P", "P_1", "P_2", "P_3")
	} else if (sample == "R1") {
		c("R1", "R1_2", "R1_3", "R1_4")
	} else {
		c("R4", "R4_1", "R4_2", "R4_3")
	}

	SNV_medicc_sub_info <- tibble(
		node = sapply(sub_node_names, function(n) sub_df$node[sub_df$label == n][1]),
		name = sub_node_names,
		fill_col = unname(node_colours[sub_node_names]),
		SNV_branch = "SNV_T"
	) %>%
		arrange(name)

	missing_sub_df <- tibble(
		node = sub_df$node,
		name = NA_character_,
		fill_col = NA_character_,
		SNV_branch = "SNV_F"
	)

	SNV_medicc_sub_info <- bind_rows(
		SNV_medicc_sub_info,
		anti_join(missing_sub_df, SNV_medicc_sub_info, by = "node")
	)

	# Add extra space below the tree so the molecular-distance scale bar does not overlap the phylogeny.
	y_lower <- min(sub_df$y, na.rm = TRUE) - 100

	pdf(
		file = paste0(output.dir, "SNV_medicc_10X_", tree, "_tree_ggtree_scale", scaling_ratio, "_sized_nodelab_reorder_LCM_", sample, ".pdf"),
		width = 10,
		height = plot_height
	)

	p <- ggtree::ggtree(sub_tree, aes(size = SNV_branch), color = "black")

	print(
		(p %<+% SNV_medicc_sub_info) +
			geom_tippoint(aes(colour = LCM_region, shape = Type, size = Type)) +
			xlim(0, 150) +
			expand_limits(y = y_lower) +
			scale_color_manual(values = tip_colours_lcm, na.value = "grey50") +
			scale_size_manual(values = c("10X" = 1, "LCM" = 3, "SNV_F" = 0.3, "SNV_T" = 2)) +
			scale_shape_manual(values = c("10X" = 19, "LCM" = 15)) +
			geom_nodelab(aes(label = name, fill = name), geom = "label", size = 3) +
			scale_fill_manual(values = setNames(SNV_medicc_sub_info$fill_col, SNV_medicc_sub_info$name)) +
			guides(colour = "none", size = guide_legend(title = "Branch Type"), fill = "none") +
			theme_tree2()
	)

	dev.off()
}

