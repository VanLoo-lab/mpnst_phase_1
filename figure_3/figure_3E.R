################################################################################
### Generate Figure 3.e
### GnT dosage effect plot
################################################################################
# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"
# Define input directory with data and output directory to save the figure
input.dir <- paste0(zenodo.dir, "data/GnT")
output.dir <- paste0(zenodo.dir, "results/figure_3/")
#dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
setwd(output.dir)

### =========================
### 0. Load libraries and data
### =========================
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(patchwork)
library(Seurat)
library(ggrepel)

### =========================
### 1. Check GnT CN info and scRNAseq data
### =========================
#Get genes with CN info 
cell_seg_exp <- readRDS(paste0(input.dir, "/MPNST_GnT_cell_seg_exp.rds"))

# check expression in GnT scRNAseq
GnT_scrna <- readRDS(paste0(input.dir, "/GnT_scRNA_MPNST_C.rds")) # 300 cells

# only use the 130 tumor cells with shared DNA + RNA
dna_barcodes <- unique(unlist(
  lapply(cell_seg_exp, function(df) df$barcode)
))
rna_barcodes = paste0("EGAF0000", as.numeric(gsub("EGAF", "", dna_barcodes))-1)
### =========================
### 2. subset seurat object
### =========================
seurat_barcodes <- colnames(GnT_scrna)
cells_use <- colnames(GnT_scrna)[
  seurat_barcodes %in% rna_barcodes
]
GnT_scrna_sub <- subset(
  GnT_scrna,
  cells = cells_use
)

### =========================
### 3. Check gene dosage effect
### =========================
# get region
region_barcode_pair <- GnT_scrna_sub@meta.data %>%
  select(region) %>%
  mutate(barcode_rna = rownames(GnT_scrna_sub@meta.data)) %>%
  mutate(barcode = paste0("EGAF0000", as.numeric(gsub("EGAF", "", barcode_rna))+1))

region_barcode_pair_ordered <- region_barcode_pair[
  match(rna_barcodes, region_barcode_pair$barcode_rna),
]

# filter genes
MPNST_SCT_count_mtx <- GnT_scrna@assays[["SCT"]]@counts %>% as.matrix() %>% t()
MPNST_SCT_data_mtx <- GnT_scrna@assays[["SCT"]]@data %>% as.matrix() %>% t() #log-normalized

min_cells <- 10
exp_cutoff <- 10

high_exp_genes <- colnames(MPNST_SCT_count_mtx)[apply(MPNST_SCT_count_mtx, 2, median) >= exp_cutoff]
suffix = paste0("_median_", exp_cutoff)

all_seg_exp <- do.call(rbind, cell_seg_exp) %>% filter(gene_name %in% high_exp_genes, CN > 0) %>%
  left_join(region_barcode_pair, by = "barcode") %>%
  select(-barcode_rna)

### =========================
### 4. Check mean expression of gene with CN in at least min cells
### =========================
median_exp_by_CN <- all_seg_exp %>%
  group_by(gene_name, CN) %>% filter(n() >= min_cells) %>%
  summarize(median_SCT_count = median(SCT_count, na.rm = TRUE), mean_SCT_count = mean(SCT_count, na.rm = TRUE)) 

exp_outliers_mean <- median_exp_by_CN %>% left_join(median_exp_by_CN, by = "gene_name") %>% mutate(expected_SCT_count.x = mean_SCT_count.y*(CN.x/CN.y)) %>% filter(CN.x > CN.y) %>%
  mutate(exp_diff = log(abs(expected_SCT_count.x - mean_SCT_count.x), 10)) %>% arrange(desc(exp_diff)) %>% mutate(highlight_genes = paste0(gene_name, " CN", CN.x, " vs ", CN.y))
genes_to_highlight_mean <- rbind(exp_outliers_mean[1:25,])

### =========================
### 5. Plot gene dosage effect
### =========================
png(file = paste0(output.dir, "/figure_3E.png"), width = 1500, height = 1500, res = 200)
options(scipen = 999)
median_exp_by_CN %>% left_join(median_exp_by_CN, by = "gene_name") %>% mutate(expected_SCT_count.x = mean_SCT_count.y*(CN.x/CN.y)) %>% filter(CN.x > CN.y) %>%
  left_join(genes_to_highlight_mean) %>% mutate(highlight = ifelse(is.na(highlight_genes), "No", "Yes")) %>%
  ggplot(aes(x = expected_SCT_count.x, y = mean_SCT_count.x, colour = highlight)) + geom_point() + 
  geom_abline(intercept = 0, slope = 1, colour = "red") + 
  scale_x_log10(limits = c(1,100000), labels = scales::comma) + scale_y_log10(limits = c(1,100000), labels = scales::comma) +
  xlab("Mean expected raw UMI count (SCT)") + ylab("Mean observed raw UMI count (SCT)") +
  scale_color_manual(values = c("black", "black")) +
  theme_classic(base_size = 16) + theme(aspect.ratio=1, legend.position = "none")
dev.off()


