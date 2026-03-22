#012026 - Haixi GnT - Gene dosage

.libPaths(c("/rsrch6/scratch/genetics/vanloolab/Yidan/SeaDragon2_lib/",.libPaths(),"/rsrch6/scratch/genetics/vanloolab/packages/R/R4.1.0/SeaDragon2/"))

library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(patchwork)
library(Seurat)
library(ggrepel)
setwd("/rsrch6/home/genetics/vanloolab/Yidan/MPNST/phase1_reannalysis/GnT/tmp")

data_dir <- "/rsrch9/home/genetics/vanloolab/data/MPNST/phase1/Haixi/GnT/results"
out_dir <- "/rsrch6/home/genetics/vanloolab/Yidan/MPNST/phase1_reannalysis/GnT/tmp"
metadata_dir <- "/metadata/"
#Get genes with CN info [Haixi generated version]
cell_seg_exp <- readRDS(paste0(data_dir, "/gene_dosage_exp/MPNST_GnT_cell_seg_exp.rds"))
# get scRNA+GnT data
# integrated_scrna <- readRDS(paste0(data_dir, "/C/MPNST_integrated_3k.rds"))

# check expression in GnT scRNAseq
GnT_scrna <- readRDS(paste0(data_dir, "/C/MPNST_C.rds")) # 300 cells

# only use the 130 tumor cells with shared DNA + RNA
dna_barcodes <- unique(unlist(
  lapply(cell_seg_exp, function(df) df$barcode)
))
rna_barcodes = paste0("EGAF0000", as.numeric(gsub("EGAF", "", dna_barcodes))-1)
# sanity check: length(intersect(colnames(GnT_scrna), rna_barcodes)) -> should == 130

# subset seurat object
seurat_barcodes <- colnames(GnT_scrna)
cells_use <- colnames(GnT_scrna)[
  seurat_barcodes %in% rna_barcodes
]
GnT_scrna_sub <- subset(
  GnT_scrna,
  cells = cells_use
)

fos_gene <- "FOS"
fos_nearby_genes <- c(
  "FOS",
  "LINC01220",
  "JDP2-AS1",
  "JDP2",
  "BATF",
  "TMED10",
  "NEK9",
  "ZC2HC1C",
  "ACYP1",
  "MLH3",
  "EIF2B2"
)
fos_nearby_genes <- intersect(fos_nearby_genes, rownames(GnT_scrna_sub))


#####3. check gene dosage [modify Haixi's code]
####################################################################################################################################
###18/08/23 Look at relationship between CN and gene expression
####################################################################################################################################

#### get region (that can be tracked by RNA_barcode)
region_barcode_pair <- GnT_scrna_sub@meta.data %>%
  select(region) %>%
  mutate(barcode_rna = rownames(GnT_scrna_sub@meta.data)) %>%
  mutate(barcode = paste0("EGAF0000", as.numeric(gsub("EGAF", "", barcode_rna))+1))

region_barcode_pair_ordered <- region_barcode_pair[
  match(rna_barcodes, region_barcode_pair$barcode_rna),
]

table(region_barcode_pair_ordered$region)
#  P R2 R3 R4 R5 
# 11 27 23 35 34 

#### filter genes
# MPNST_RNA_count_mtx <- GnT_scrna@assays[["RNA"]]@counts %>% as.matrix() %>% t()
# MPNST_RNA_data_mtx <- GnT_scrna@assays[["RNA"]]@data %>% as.matrix() %>% t() #log-normalized
MPNST_SCT_count_mtx <- GnT_scrna@assays[["SCT"]]@counts %>% as.matrix() %>% t()
MPNST_SCT_data_mtx <- GnT_scrna@assays[["SCT"]]@data %>% as.matrix() %>% t() #log-normalized

min_cells <- 10
exp_cutoff <- 10
# high_exp_genes <- colnames(MPNST_SCT_count_mtx)[colSums(MPNST_SCT_count_mtx) > exp_cutoff]
high_exp_genes <- colnames(MPNST_SCT_count_mtx)[apply(MPNST_SCT_count_mtx, 2, median) >= exp_cutoff]
suffix = paste0("_median_", exp_cutoff)

all_seg_exp <- do.call(rbind, cell_seg_exp) %>% filter(gene_name %in% high_exp_genes, CN > 0) %>%
  left_join(region_barcode_pair, by = "barcode") %>%
  select(-barcode_rna)

#Check mean expression of gene with CN in at least min cells
median_exp_by_CN <- all_seg_exp %>%
  group_by(gene_name, CN) %>% filter(n() >= min_cells) %>%
  summarize(median_SCT_count = median(SCT_count, na.rm = TRUE), mean_SCT_count = mean(SCT_count, na.rm = TRUE)) 

## regenerate figure in manuscript (old version)
# no highlight
png(filename = paste0(out_dir, "/20260115_output/20260115_MPNST_GnT_SCT_count_median_cells_vs_expected_min_cells_", min_cells, suffix, ".png"), width = 2000, height = 2000, res = 200)
median_exp_by_CN %>% left_join(median_exp_by_CN, by = "gene_name") %>% mutate(expected_SCT_count.x = median_SCT_count.y*(CN.x/CN.y)) %>% filter(CN.x > CN.y) %>%
  ggplot(aes(x = expected_SCT_count.x, y = median_SCT_count.x)) + geom_point() + 
  geom_abline(intercept = 0, slope = 1, colour = "red") + 
  scale_x_log10(limits = c(1,100000)) + scale_y_log10(limits = c(1,100000)) +
  xlab("Median expected raw UMI count (SCT)") + ylab("Median raw UMI count (SCT)") +
  theme_classic(base_size = 16) + theme(aspect.ratio=1)
dev.off()

###### use mean ###### <- the current figure 
exp_outliers_mean <- median_exp_by_CN %>% left_join(median_exp_by_CN, by = "gene_name") %>% mutate(expected_SCT_count.x = mean_SCT_count.y*(CN.x/CN.y)) %>% filter(CN.x > CN.y) %>%
  mutate(exp_diff = log(abs(expected_SCT_count.x - mean_SCT_count.x), 10)) %>% arrange(desc(exp_diff)) %>% mutate(highlight_genes = paste0(gene_name, " CN", CN.x, " vs ", CN.y))
genes_to_highlight_mean <- rbind(exp_outliers_mean[1:25,])


############Plot new figure 3E
#pdf(file = paste0(out_dir, "/20260112_MPNST_GnT_SCT_count_mean_cells_vs_expected_lab_min_cells_", min_cells, suffix, ".pdf"), width = 7, height = 7)
png(file = paste0(out_dir, "/20260112_MPNST_GnT_SCT_count_mean_cells_vs_expected_lab_min_cells_", min_cells, suffix, ".png"), width = 1500, height = 1500, res = 200)

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
options(scipen = 0)

# the "median_exp_by_CN" is actually mean_exp_by_CN here. was in a rush and didn't correct it..