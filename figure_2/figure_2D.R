### title: recreate Figure 2D VAF overlay for R1 back section

####################################################################################################################################
### Part 0: Define input and output directories
####################################################################################################################################

zenodo.dir <- normalizePath("~/Documents/GitHub/MPNST-Zenodo", mustWork = FALSE)
input.dir <- file.path(zenodo.dir, "data")
output.dir <- file.path(zenodo.dir, "results", "figure_2")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

####################################################################################################################################
### Part 1: Load libraries
####################################################################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(jpeg)
  library(scales)
})

####################################################################################################################################
### Part 2: Define and validate required input files
####################################################################################################################################

vaf_cluster_path <- Sys.getenv(
  "R1_VAF_CLUSTER_PATH",
  unset = file.path(input.dir, "LCM", "R1_VAF_cluster.rds")
)

required_files <- c(
  file.path(input.dir, "LCM", "LCM_data_final_manual.csv"),
  file.path(input.dir, "LCM", "barcodes.tsv"),
  file.path(input.dir, "LCM", "Genotype_SNP", "both_haplo_count_combined_by_cell.rds"),
  file.path(input.dir, "images", "R1", "T1.1_B_Overview_cut_positions.jpg"),
  vaf_cluster_path
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing required input file(s):\n",
      paste0(" - ", missing_files, collapse = "\n"),
      "\nUpdate `vaf_cluster_path` in this script or set `R1_VAF_CLUSTER_PATH` ",
      "to the downloaded precomputed matrix path under MPNST-Zenodo/data."
    ),
    call. = FALSE
  )
}

####################################################################################################################################
### Part 3: Load minimal metadata required for plotting
####################################################################################################################################

lcm_locations <- read.delim(file.path(input.dir, "LCM", "LCM_data_final_manual.csv"), sep = ",")
lcm_barcodes <- read.delim(file.path(input.dir, "LCM", "barcodes.tsv")) %>%
  mutate(region = ifelse(str_detect(sample, "006"), "Primary", ifelse(str_detect(sample, "005"), "T1.1", "T4.1")))

haplo_counts <- readRDS(file.path(input.dir, "LCM", "Genotype_SNP", "both_haplo_count_combined_by_cell.rds")) %>%
  as_tibble()
names(haplo_counts) <- trimws(names(haplo_counts))

lcm_ascn_profile_samples <- tibble(
  sample = unique(haplo_counts$Barcode),
  CN_profile = "Yes"
)

lcm_coordinates <- tibble(
  barcode = gsub("\\)", "", gsub(".*\\(", "", lcm_locations$UDF.Index)),
  coord.x = as.numeric(gsub(",.*", "", lcm_locations$Coordinates..X.Y..mm)),
  coord.y = as.numeric(gsub(".*,", "", lcm_locations$Coordinates..X.Y..mm)),
  region = lcm_locations$Region,
  side = lcm_locations$Side,
  sample_number = lcm_locations$Sample.Name.LCM.area
) %>%
  left_join(lcm_barcodes, by = c("barcode", "region")) %>%
  left_join(lcm_ascn_profile_samples, by = "sample") %>%
  mutate(CN_profile = replace_na(CN_profile, "No"))

####################################################################################################################################
### Part 4: Load the precomputed R1 VAF matrix and define figure-specific settings
####################################################################################################################################

vaf_cluster <- readRDS(vaf_cluster_path)
if (is.null(rownames(vaf_cluster)) || is.null(colnames(vaf_cluster))) {
  stop("`R1_VAF_cluster` must have both row names (clusters) and column names (spots).", call. = FALSE)
}

s <- "R1"
section <- "T1.1_B -Back"
section_side <- gsub(".*_", "", section)
expected_clusters <- c("R1 R5", "R1_1", "R1_2", "R1_3")

####################################################################################################################################
### Part 5: Apply the original R1 cluster selection logic and validate inputs
####################################################################################################################################

clusters_in_region <- grep(s, rownames(vaf_cluster), value = TRUE)
if (length(clusters_in_region) < 5) {
  stop(
    paste0(
      "Expected at least 5 R1-related clusters before applying the original selection logic, but found: ",
      paste(clusters_in_region, collapse = ", ")
    ),
    call. = FALSE
  )
}
clusters_in_region <- clusters_in_region[c(2, 3:5)]

if (!identical(clusters_in_region, expected_clusters)) {
  stop(
    paste0(
      "The original R1 cluster selection did not yield the expected Figure 2D facets.\n",
      "Expected: ", paste(expected_clusters, collapse = ", "), "\n",
      "Observed: ", paste(clusters_in_region, collapse = ", ")
    ),
    call. = FALSE
  )
}

missing_spots <- setdiff(colnames(vaf_cluster), lcm_coordinates$sample)
if (length(missing_spots) > 0) {
  stop(
    paste0(
      "The following VAF matrix spot(s) are missing from `LCM_coordinates`: ",
      paste(missing_spots, collapse = ", ")
    ),
    call. = FALSE
  )
}

####################################################################################################################################
### Part 6: Build the plotting data frame for the R1 back section
####################################################################################################################################

plot_df <- lcm_coordinates %>%
  filter(sample %in% colnames(vaf_cluster), side == section_side) %>%
  left_join(
    t(vaf_cluster[clusters_in_region, , drop = FALSE]) %>%
      as.data.frame() %>%
      rownames_to_column("sample"),
    by = "sample"
  ) %>%
  pivot_longer(
    cols = all_of(clusters_in_region),
    names_to = "cluster",
    values_to = "VAF"
  ) %>%
  mutate(cluster = factor(cluster, levels = clusters_in_region))

if (nrow(plot_df) == 0) {
  stop(
    paste0(
      "No section-specific spots were found for section `", section, "` (side `", section_side, "`)."
    ),
    call. = FALSE
  )
}

if (anyNA(plot_df$VAF)) {
  stop("`plot_df$VAF` contains missing values after joining the VAF matrix to coordinates.", call. = FALSE)
}

####################################################################################################################################
### Part 7: Load and crop the histology image using the original paper settings
####################################################################################################################################

lcm_image <- readJPEG(file.path(input.dir, "images", "R1", "T1.1_B_Overview_cut_positions.jpg"))

plot_w <- 1000
plot_h <- 1070
crop_xmin <- 12
crop_xmax <- 92
crop_ymin <- 0
crop_ymax <- 80

row_idx <- seq.int(
  from = max(1, round(nrow(lcm_image) * crop_ymin / 100)),
  to = min(nrow(lcm_image), round(nrow(lcm_image) * crop_ymax / 100))
)
col_idx <- seq.int(
  from = max(1, round(ncol(lcm_image) * crop_xmin / 100)),
  to = min(ncol(lcm_image), round(ncol(lcm_image) * crop_xmax / 100))
)
lcm_image_crop <- lcm_image[row_idx, col_idx, , drop = FALSE]

####################################################################################################################################
### Part 8: Draw and save Figure 2D
####################################################################################################################################

output_file <- file.path(output.dir, "figure_2D.png")
png(filename = output_file, width = plot_w, height = plot_h, res = 200)
print(
  ggplot(plot_df, aes(x = coord.x, y = coord.y, color = VAF)) +
    annotation_raster(lcm_image_crop, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    geom_point(size = 2.8) +
    scale_colour_gradientn(
      colours = c("white", "blue3", "navy"),
      limits = c(0, 0.5),
      breaks = seq(0, 0.5, by = 0.25),
      oob = squish
    ) +
    labs(color = "VAF") +
    coord_cartesian(
      xlim = c(crop_xmin, crop_xmax),
      ylim = c(100 - crop_ymax, 100 - crop_ymin)
    ) +
    theme(
      text = element_text(size = 12),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "top",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-10, -10, -10, -10)
    ) +
    facet_wrap(~cluster)
)
dev.off()

message("Saved Figure 2D to ", output_file)
