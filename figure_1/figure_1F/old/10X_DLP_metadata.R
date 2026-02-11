#Script containing common metadata for 10X_DNA and DLP
# library(tidyverse, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")

####################################################################################################################################
#10X info
####################################################################################################################################
samples = c("R1", "R2", "R3", "R4", "R5", "P")
names(samples) = c("FIT208A3", "FIT208A4", "FIT208A5", "FIT208A6", "FIT208A7", "FIT208A8")
# SUBSETS = c("P", "R1", "R2", "R3", "R4", "R5", "G2MS", "Ribo")

bulk_samples = c("VER236A1", "VER236A2", "VER236A3", "VER236A4", "VER236A5", "VER236A6")
names(bulk_samples) <- samples

#Colours
region_colours <- c("#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3", "#F8766D")[c(6,1:5)]
names(region_colours) <- samples[c(6,1:5)]

#File locations
raw.scDNA.dir = "/camp/project/proj-vanloo/analyses/averfaillie/10x_DNA/"
region.folders = c("SC18143_1/FIT208A3_b/",
                   "SC18143_1/FIT208A4_b/",
                   "sc18143_2b/FIT208A5_combo/",
                   "sc18143_2b/FIT208A6_combo/",
                   "sc18143_2b/FIT208A7_combo/",
                   "sc18143_2b/FIT208A8_combo/")
input.dir = c(paste0(raw.scDNA.dir, region.folders,"outs/"))
names(input.dir) <- samples

#R1:150gb 2220 cells
#R2:122.1gb 1440 cells
#R3:14.3gb 106 cells
#R4:50gb 363 cells
#R5:139.5gb 1513 cells
#P:99.7gb 1154 cells

barcode.dir <- paste0(input.dir,"per_cell_summary_metrics.csv")
barcodes <- list()
for (s in 1:length(samples)) {
  barcodes[[s]] <- as.character(read.csv(barcode.dir[s])[,"barcode"])
  barcodes[[s]] <- paste0(samples[s], "_10X_", barcodes[[s]])
}
# for (s in 1:length(samples)) {
#   barcodes[[s]] <- as.character(read.csv(barcode.dir[s])[,"barcode"])
#   barcodes[[s]] <- gsub(paste0("(................)-1"), paste0(samples[s], "_\\1"),barcodes[[s]])
# }

chrom_numbers <- c(1:23)
chrom_names <- c(1:22,"X")

region_coverage.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/results/Genotype_de_novo_SNVs/")
pcf.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/results/ASCAT.sc/")
asCN.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/results/ASCAT.sc/asCN/")
medicc.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/results/MEDICC2/")

####################################################################################################################################
#DLP info
####################################################################################################################################
#R1 576 cells
#R2 627 cells
#R3 417 cells
#R4 584 cells
#R5 488 cells
#P 308 cells

fastq.dir.1 = "/camp/project/proj-vanloo/analyses/hyan/mpnst/DLP_plus/data/PM22040/220401_A01366_0166_AH2CNFDRX2/fastq" #R1R2
fastq.dir.2 = "/camp/project/proj-vanloo/analyses/hyan/mpnst/DLP_plus/data/PM22040/220228_A01366_0152_AHTTMYDRXY/fastq" #R3R4
fastq.dir.3 = "/camp/project/proj-vanloo/analyses/hyan/mpnst/DLP_plus/data/PM22040/220909_A01366_0277_AHKF5CDRX2/fastq" #R5P

#Load location
chip_sample_sheets <- list(LES4677A1 = read.csv("/camp/project/proj-vanloo/analyses/hyan/mpnst/DLP_plus/data/PM22040_LES4677A1_124574_sample_sheet.csv", stringsAsFactors = F),
                           LES4677A2 = read.csv("/camp/project/proj-vanloo/analyses/hyan/mpnst/DLP_plus/data/PM22040_LES4677A2_124576_sample_sheet.csv", stringsAsFactors = F),
                           LES4677A4 = read.csv("/camp/project/proj-vanloo/analyses/hyan/mpnst/DLP_plus/data/PM22040_LES4677A4_124668_sample_sheet.csv", stringsAsFactors = F))
chip_layouts <- lapply(chip_sample_sheets, function(s) {
  data.frame(x = as.numeric(gsub(".*([[:digit:]][[:digit:]])x.*", "\\1", s$Sample_ID)),
             y = as.numeric(gsub(".*([[:digit:]][[:digit:]])y", "\\1", s$Sample_ID)),
             sample = s$User_Sample_Name, stringsAsFactors = F)
})

chip_layouts_mtx <- lapply(chip_layouts, function(s) {
  mtx <- matrix(nrow = 72, ncol = 72)
  for (b in 1:nrow(s)) {
    mtx[s$x[b], s$y[b]] <- which(s$sample[b] == unique(s$sample))
  }
  return(mtx)
})
# heatmap(chip_layouts_mtx[[1]], Colv = NA, Rowv = NA, revC = T)

#Get barcodes/fastqs
alt_names <- c("T1.1", "T2.2", "T3.2", "T4.1", "T5.1", "P")
names(alt_names) <- samples

DLP_barcode.dir <- "/camp/project/proj-swanton-vanloo/working/DLP_plus/barcodes/"
DLP_barcodes <- lapply(samples[1:4], function(b) {
  barcode_tsv <- read.delim(paste0(DLP_barcode.dir, alt_names[b], "_barcodes.tsv"))
  barcode_tsv$FASTQ.DIR = ifelse(grepl("LES4677A1", barcode_tsv$sample_id), fastq.dir.1, ifelse(grepl("LES4677A2", barcode_tsv$sample_id), fastq.dir.2, fastq.dir.3))
  return(barcode_tsv)
})
# DLP_barcodes <- c(DLP_barcodes,
#                   list(data.frame(sample_id = chip_sample_sheets[[3]] %>% filter(User_Sample_Name == "T5.1") %>% pull(Sample_ID),
#                                   BARCODE = DLP_barcodes[[2]] %>% arrange(sample_id) %>% pull(BARCODE),#using R2 DLP_barcodes
#                                   FASTQ.DIR = fastq.dir.3, stringsAsFactors = F)),
#                   list(data.frame(sample_id = chip_sample_sheets[[3]] %>% filter(User_Sample_Name == "primary ") %>% pull(Sample_ID),
#                                   BARCODE = DLP_barcodes[[1]] %>% arrange(sample_id) %>% pull(BARCODE),#using R1 DLP_barcodes
#                                   FASTQ.DIR = fastq.dir.3, stringsAsFactors = F)))
# names(DLP_barcodes) = samples

DLP_barcodes <- c(DLP_barcodes,
                  list(data.frame(sample_id = chip_sample_sheets[[3]][chip_sample_sheets[[3]]$User_Sample_Name == "T5.1","Sample_ID"],
                                  BARCODE = DLP_barcodes[[2]][order(DLP_barcodes[[2]]$sample_id),"BARCODE"],#using R2 DLP_barcodes
                                  FASTQ.DIR = fastq.dir.3, stringsAsFactors = F)),
                  list(data.frame(sample_id = chip_sample_sheets[[3]][chip_sample_sheets[[3]]$User_Sample_Name == "primary ","Sample_ID"],
                                  BARCODE = DLP_barcodes[[1]][order(DLP_barcodes[[2]]$sample_id),"BARCODE"],#using R1 DLP_barcodes
                                  FASTQ.DIR = fastq.dir.3, stringsAsFactors = F)))
names(DLP_barcodes) = samples

all_barcodes <- unlist(lapply(1:length(samples), function(s){
  return(paste0(samples[s], "_", as.character(DLP_barcodes[[samples[s]]]$BARCODE)))
}))

all_DLP_barcodes <- unlist(lapply(1:length(samples), function(s){
  return(paste0(samples[s], "_DLP_", as.character(DLP_barcodes[[samples[s]]]$sample_id), ".bam"))
}))

#File locations
DLP.hg19.dir = "/camp/project/proj-swanton-vanloo/working/DLP_plus/bam/"
DLP.dir = "/camp/project/proj-vanloo/analyses/hyan/mpnst/DLP_plus/data/hg38_bam/"

DLP.singlebam.dir <- paste0(DLP.dir, samples)
names(DLP.singlebam.dir) <- samples

DLP.input.dir <- paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/DLP_plus/data/modified_BAMs/", samples,"/")
names(DLP.input.dir) <- samples

####################################################################################################################################
#Cluster info generated after running kmeans in ASCAT.sc.R then cleaned up
####################################################################################################################################
tenX_DLP_ASCAT.dir <- "/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/results/ASCAT.sc/"

K = 22
run <- "22_kmeans_10X_DLP"
remove_kmeans2 = T
scDNA_CN_cluster_ids <- readRDS(paste0(tenX_DLP_ASCAT.dir, "MPNST_all_k_means_K",K,"_clusters.rds"))
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


if (F) {
  allele_freq_all_regions <- list()
  #Note very slow
  for (s in samples) {
    print(s)
    allele_freq_all_regions[[s]] <- do.call(rbind, lapply(chrom_numbers, function(i) {
      readRDS(paste0("../Genotype_SNPs/",s,"/haplotype_count_",s,"_chr",i,".rds"))[[i]]
    }))
  }
  
  DLP_clusters_coverage <- lapply(scDNA_CN_cluster_ids, function(k) {
    cluster_coverage <- lapply(samples, function(s){
      return(allele_freq_all_regions[[s]][!duplicated(allele_freq_all_regions[[s]][,1:2]),] %>% 
               left_join() %>% filter() %>% pull(Coverage) %>% median()) #Use Coverage column to calculate median coverage at positions
    })
    names(region_coverage) <- samples
  })
  saveRDS(DLP_clusters_coverage, "MPNST_DLP_clusters_coverage.rds")
}

