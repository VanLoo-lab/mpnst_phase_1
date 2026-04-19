### title: plotting scAIDD doublets
### Code to generate Figure 3b
# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"
#zenodo.dir <- "/Users/ycheng3/Projects/MPNST_phase1/"
# Define input directory with data and output directory to save the figure
input.dir <- paste0(zenodo.dir, "data/scAIDD")
output.dir <- paste0(zenodo.dir, "results/figure_3/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
setwd(output.dir)

### =========================
### 0. Load libraries and data
### =========================
library(dplyr)
library(ggplot2)

both_haplo_count_combined_by_cell <- readRDS(paste0(input.dir,"/both_haplo_count_combined_by_cell.rds"))
both_haplo_count_doublet_by_cell <- readRDS(paste0(input.dir,"/MPNST_D_both_haplo_count_doublet_by_cell.rds"))
haplo_scatter_normal<-readRDS(paste0(input.dir,"/MPNST_C_haplo_scatter_group_1_ordered_ids.rds"))
haplo_scatter_doublet<-readRDS(paste0(input.dir,"/MPNST_C_haplo_scatter_group_2_ordered_ids.rds"))
haplo_scatter_malignant<-readRDS(paste0(input.dir,"/MPNST_C_haplo_scatter_group_3_ordered_ids.rds"))


### =========================
### 1. add combined_type to real cells
### =========================
both_haplo_count_combined_by_cell <- both_haplo_count_combined_by_cell %>%
  mutate(
    combined_type = case_when(
      Barcode %in% haplo_scatter_normal ~ "Normal cells",
      Barcode %in% haplo_scatter_doublet ~ "Suspected doublets",
      Barcode %in% haplo_scatter_malignant ~ "Malignant cells",
      TRUE ~ NA_character_
    )
  )

### =========================
### 2. add combined_type to simulated doublets
### =========================
both_haplo_count_doublet_by_cell <- both_haplo_count_doublet_by_cell %>%
  mutate(combined_type = "Simulated doublets")

### =========================
### 3. merge both tables
### =========================
haplo_plot_df <- bind_rows(
  both_haplo_count_combined_by_cell,
  both_haplo_count_doublet_by_cell
)

### =========================
### 4. set factor order
### =========================
haplo_plot_df$combined_type <- factor(
  haplo_plot_df$combined_type,
  levels = c("Normal cells", "Malignant cells", "Suspected doublets", "Simulated doublets")
)

### =========================
### 5. color palette
### =========================
combined_colors <- c(
  "Normal cells" = "#00BFC4",
  "Malignant cells" = "#F8766D",
  "Suspected doublets" = "blue",
  "Simulated doublets" = "magenta"
)

### =========================
### 6. haplotype scatter plot
### =========================
png(file = paste0(output.dir, "/figure_3B.png"), width = 1850, height = 1600, res = 300)
ggplot(
  haplo_plot_df,
  aes(x = Haplo_1_Count, y = Haplo_2_Count, color = combined_type)
) +
  geom_point(size = 0.1, alpha = 0.7) +
  scale_color_manual(values = combined_colors) +
  scale_x_log10(limits = c(10,50000)) + scale_y_log10(limits = c(10,50000)) +
  theme_classic(base_size = 14) +
  labs(
    x = "Haplotype 1 count",
    y = "Haplotype 2 count",
    color = "",
    title = "Haplotype scatter by cell type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
dev.off()
