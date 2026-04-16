### title: Plot LCM MEDICC2 trees with mapped spots for T4.1 region sides

# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"

# Define input and output directories
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
library(jpeg)

####################################################################################################################################
### Part 1: Inline metadata (minimal subset for T4.1 LCM analysis)
####################################################################################################################################

### Region definitions
LCM.dir <- paste0(input.dir, "LCM/")

# Sides for T4.1 region
sides <- c("Front", "Back", "Side")
region_name <- "T4.1"
min_seg <- 2.5
tumour_name <- paste0("LCM_sides_", min_seg)

# User-provided locations for required inputs
lcm_image_dir <- paste0(LCM.dir, "images/")
medicc_tree_dir <- paste0(LCM.dir, "MEDICC2_slides/")

### Load LCM coordinate data
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
	filter(sample != "VER315A110_S319_L007") %>%
	filter(region == region_name)

####################################################################################################################################
### Part 2: Load LCM images
####################################################################################################################################

LCM_images_list <- list(
	"T4.1_Front" = "T4.1_Front_OVERVIEW_CUT_6.3x_positions.jpeg",
	"T4.1_Back" = "T4.1_Back_OV_20x_cut.jpg",
	"T4.1_Side" = "T4.1_side_OV_cut_20x.jpg"
)

side_map <- c(
	"Front" = "A -Front",
	"Back" = "B -Back",
	"Side" = "C -Side"
)

LCM_image_jpg <- list()
image_paths_found <- c()

for (side_name in sides) {
	sample_label <- paste0(region_name, "_", side_name)
	side_suffix <- tolower(side_name)
	img_name <- LCM_images_list[[sample_label]]
	
	# Try multiple possible paths for the images
	possible_paths <- c(
		paste0(lcm_image_dir, img_name)
	)
	
	img_found <- FALSE
	for (path in possible_paths) {
		if (file.exists(path)) {
			LCM_image_jpg[[sample_label]] <- readJPEG(path)
			image_paths_found <- c(image_paths_found, path)
			img_found <- TRUE
			break
		}
	}
	
	if (!img_found) {
		warning(paste("Image not found for", sample_label, "in paths:"))
		for (path in possible_paths) warning(paste("  ", path))
	}
}

####################################################################################################################################
### Part 3: Generate figures for each T4.1 side
####################################################################################################################################

run_name <- "wgd"

for (side_name in sides) {
	sample_label <- paste0(region_name, "_", side_name)
	side_suffix <- tolower(side_name)
	
	# Construct the tree file path
	alt_tree_paths <- c(
		paste0(medicc_tree_dir, tumour_name, "_", sample_label, "_final_tree.new")
	)
	
	tree_file_exists <- FALSE
	actual_tree_file <- NULL
	
	for (path in alt_tree_paths) {
		if (file.exists(path)) {
			tree_file_exists <- TRUE
			actual_tree_file <- path
			break
		}
	}
	
	if (!tree_file_exists) {
		warning(paste("MEDICC2 tree file not found for", sample_label))
		warning(paste("Expected paths:"))
		for (path in alt_tree_paths) warning(paste("  ", path))
		next
	}
	
	# Load the tree
	medicc_LCM_tree <- read.tree(file = actual_tree_file)
	
	####################################################################################################################################
	### Part 3a: Generate location edges for tree visualization
	####################################################################################################################################
	
	# Generate location of intermediate nodes based on LCM coordinates
	all <- c(medicc_LCM_tree[["tip.label"]], "root", medicc_LCM_tree[["node.label"]][-1])
	nodes <- medicc_LCM_tree[["node.label"]][-1]
	edges <- data.frame(
		parent = all[medicc_LCM_tree[["edge"]][,1]],
		child = all[medicc_LCM_tree[["edge"]][,2]],
		length = medicc_LCM_tree[["edge.length"]]
	) %>%
		left_join(
			LCM_coordinates_tumour %>%
				filter(side == side_map[[side_name]]) %>%
				select(sample, coord.x, coord.y) %>%
				distinct(),
			by = c("child" = "sample")
		) %>%
		filter(str_detect(child, "diploid", negate = TRUE)) %>%
		mutate(edge = 1:nrow(.))
	
	# Iterate to fill in coordinates for internal nodes
	max_iterations <- 100
	iteration <- 0
	while(any(is.na(edges$coord.x)) && iteration < max_iterations) {
		for (i in 1:nrow(edges)) {
			node <- edges[i, "parent"]
			children <- edges %>% filter(parent == node) %>% pull(child)
			
			if (length(children) >= 2) {
				child_x_vals <- edges %>% filter(child %in% children) %>% pull(coord.x) %>% na.omit()
				child_y_vals <- edges %>% filter(child %in% children) %>% pull(coord.y) %>% na.omit()
				
				if (length(child_x_vals) >= 2) {
					edges[edges$child == node, "coord.x"] <- mean(child_x_vals)
					edges[edges$child == node, "coord.y"] <- mean(child_y_vals)
				}
			}
		}
		iteration <- iteration + 1
	}
	
	####################################################################################################################################
	### Part 3b: Adjust internal node positions for clarity
	####################################################################################################################################
	
	pct_adj <- 0.3
	edges_adj <- edges
	
	for (i in nrow(edges):1) {
		current_node <- edges_adj[i, "child"]
		if (grepl("internal", current_node)) {
			child_x <- edges_adj[edges_adj[, "child"] == current_node, "coord.x"]
			child_y <- edges_adj[edges_adj[, "child"] == current_node, "coord.y"]
			
			parent_node <- edges_adj[edges_adj[, "child"] == current_node, "parent"]
			if (parent_node != "root" && !is.na(child_x)) {
				parent_x <- edges_adj[edges_adj[, "child"] == parent_node, "coord.x"]
				parent_y <- edges_adj[edges_adj[, "child"] == parent_node, "coord.y"]
				
				if (!is.na(parent_x)) {
					adj_x <- child_x + (parent_x - child_x) * pct_adj
					adj_y <- child_y + (parent_y - child_y) * pct_adj
					edges_adj[edges_adj[, "child"] == current_node, "coord.x"] <- adj_x
					edges_adj[edges_adj[, "child"] == current_node, "coord.y"] <- adj_y
				}
			}
		}
	}
	
	# Rename root node to MRCA
	edges_adj <- edges_adj %>%
		mutate(parent = str_replace(parent, paste0(edges_adj[1, "child"], "$"), "MRCA"))
	edges_adj[1, "child"] <- "MRCA"
	
	# Add tree level information
	edges_adj$tree_level <- unlist(lapply(1:nrow(edges_adj), function(n) {
		level <- 0
		child <- edges_adj[n, "child"]
		while (child != "root") {
			child <- edges_adj[edges_adj[, "child"] == child, "parent"]
			level <- level + 1
		}
		return(level)
	}))
	
	####################################################################################################################################
	### Part 3c: Prepare plot data and generate graphics
	####################################################################################################################################
	
	# Create edges for plotting
	edges_ggplot <- rbind(
		edges_adj[, c(2, 4, 5, 6)] %>% mutate(direction = "end"),
		edges_adj[, c(1, 6)] %>%
			left_join(edges_adj[, c(2, 4, 5)], by = c("parent" = "child")) %>%
			select(1, 3, 4, 2) %>%
			rename(child = parent) %>%
			mutate(direction = "start")
	) %>%
		mutate(type = ifelse(str_detect(child, "VER"), "Tumour", ifelse(str_detect(child, "MRCA"), "MRCA", "Node"))) %>%
		filter(str_detect(child, "root", negate = TRUE)) %>%
		left_join(edges_adj[, c(2, 7)], by = "child")
	
	# Adjust MRCA position
	adj_distance <- 15
	orig_mrca_pos <- c(
		edges_ggplot %>% filter(direction == "end", type == "MRCA") %>% pull(coord.x),
		edges_ggplot %>% filter(direction == "end", type == "MRCA") %>% pull(coord.y)
	)
	edges_ggplot[which(edges_ggplot$direction == "end" & edges_ggplot$type == "MRCA"), "coord.y"] <- orig_mrca_pos[2] - adj_distance
	edges_ggplot <- edges_ggplot %>%
		mutate(type = ifelse(type == "MRCA", ifelse(direction == "end", "MRCA", "MRCA_hide"), type))
	
	mrca_line <- data.frame(
		coord.x = c(orig_mrca_pos[1], edges_ggplot %>% filter(direction == "end", type == "MRCA") %>% pull(coord.x)),
		coord.y = c(orig_mrca_pos[2], edges_ggplot %>% filter(direction == "end", type == "MRCA") %>% pull(coord.y))
	)
	
	####################################################################################################################################
	### Part 3d: Plot overlay with colored spots
	####################################################################################################################################
	
	edges_ggplot_col <- edges_ggplot %>%
		mutate(unique_id = ifelse(type == "MRCA", "MRCA", ifelse(type == "Tumour", child, NA)))
	
	point_colors <- c("black", rainbow(length(unique(edges_ggplot_col$unique_id[!is.na(edges_ggplot_col$unique_id)])) - 1))
	names(point_colors) <- na.omit(unique(edges_ggplot_col$unique_id))
	
	# Plot overlay on LCM image if available
	if (!is.null(LCM_image_jpg[[sample_label]])) {
		png(
			filename = paste0("figure_2E_image_", side_suffix, ".png"),
			width = 6000, height = 6000, res = 200
		)
		
		print(edges_ggplot_col %>%
			ggplot(aes(x = coord.x, y = coord.y)) +
			annotation_raster(
				LCM_image_jpg[[sample_label]],
				xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
			) +
			geom_point(size = 30, aes(color = unique_id)) +
			scale_color_manual(values = point_colors) +
			geom_path(
				mapping = aes(group = edge),
				data = edges_ggplot %>% arrange(edge, desc(direction)),
				arrow = arrow(length = unit(30, "pt")),
				colour = "navy", size = 3
			) +
			geom_path(data = mrca_line, colour = "navy", size = 3) +
			coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
			theme(
				text = element_text(size = 32),
				axis.title = element_blank(),
				axis.text = element_blank(),
				axis.ticks = element_blank(),
				legend.position = "none"
			)
		)
		dev.off()
	}
	
	####################################################################################################################################
	### Part 3e: Plot colored tree
	####################################################################################################################################
	
	mrca_node <- as_tibble(medicc_LCM_tree) %>% filter(str_detect(label, "internal")) %>% pull(node) %>% min()
	tree_metadata <- tibble(label = unique(edges_ggplot_col$child[!is.na(edges_ggplot_col$unique_id)]))
	
	lcm_ggtree <- ggtree(tidytree::tree_subset(medicc_LCM_tree, mrca_node, levels_back = 0))
	
	mrca_coords <- lcm_ggtree$data[lcm_ggtree$data$node ==
		lcm_ggtree$data %>% filter(str_detect(label, "internal")) %>% pull(node) %>% min(),
		c("x", "y")
	]
	
	pdf(
		file = paste0("figure_2E_tree_", side_suffix, ".pdf"),
		width = 3, height = 7
	)
	
	print(lcm_ggtree %<+% tree_metadata +
		geom_tippoint(aes(colour = label), size = 5) +
		geom_point(data = mrca_coords, aes(x = x, y = y), shape = 19, size = 5, color = "black") +
		scale_color_manual(values = point_colors) +
		theme_tree2() +
		theme(legend.position = "none")
	)
	dev.off()
	
	cat(paste0("\nGenerated figures for ", sample_label, "\n"))
}

cat("\n=== Figure generation complete ===\n")










