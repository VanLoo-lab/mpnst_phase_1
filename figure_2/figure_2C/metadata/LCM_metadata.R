#Script containing common metadata for LCM
library(tidyverse)
library(jpeg)

####################################################################################################################################
### 02/09/2021 First attempt at Shiny
####################################################################################################################################
LCM_samples = c("R1", "R4", "P")
names(LCM_samples) <- c("T1.1", "T4.1", "Primary")

regions <- list(c("Primary", "A -Front"),
                c("Primary", "B -Back"),
                c("Primary", "C -Side"),
                c("T1.1", "A -Front"),
                c("T1.1", "B -Back"),
                c("T1.1", "C -Side"),
                c("T4.1", "A -Front"),
                c("T4.1", "B -Back"),
                c("T4.1", "C -Side"))

#006 is primary (48), #005 is T1 (44), #007 is T4 (42)
# pcf.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DNA/results/pcf/")
# LCM.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/LCM/")
# LCM.bam.dir <- "/camp/lab/vanloop/working/mtarabichi/projects/spatial_genomics/input/"
# scDNA_medicc.dir = "/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DNA/results/MEDICC2/scDNA_20_kmeans_R3mod_1.5/"
# GnT_medicc.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/GnT/results/MEDICC2/MPNST_1_10X_LCM_GnT/")
# asCN.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DNA/results/pcf/asCN/")

LCM.dir <- paste0(INPUTDIR, "LCM/")

# LCM_locations <- lapply(regions, function(s) {
#   read.delim(paste0(LCM.dir, "data/LCM_data_final.csv"), sep = ",") %>% filter(Region == s[1], Side == s[2])
# })

LCM_locations <- read.delim(paste0(LCM.dir, "LCM_data_final_manual.csv"), sep = ",")

LCM_barcodes <- read.delim(paste0(LCM.dir, "barcodes.tsv")) %>% 
  mutate(region = ifelse(str_detect(sample, "006"), "Primary", ifelse(str_detect(sample, "005"), "T1.1", "T4.1")))

### Alex Stein: Problem, I cannot find the files with "_ASCAT.sc.txt"

# LCM_asCN_profiles.dir <- "/camp/lab/vanloop/working/mtarabichi/projects/spatial_genomics/lcm_profiles/"
# LCM_asCN_profile_samples <- data.frame(sample = gsub("_ASCAT.sc.txt", "", list.files(LCM_asCN_profiles.dir)),
#                                        CN_profile = "Yes")

### Alex Stein: Here is an attempt to retrieve the information from the allele_freq files
base_dir <- "/Users/alexanderstein/Hyperion/people/mtarabichi/MPNST/pvanloo_20250618/Haixi/LCM/results/Genotype_SNPs/Allele_Freq/"   # folder containing P, R1, R2
# list all .rds files ending with L00x.rds in subfolders
files <- list.files(base_dir, pattern = "L00[567]\\.rds$", recursive = TRUE, full.names = FALSE)
# extract sample name (between 'allele_freq_' and '_L00x.rds')
sample_names <- sub("^allele_freq_(.*)_L00[567]\\.rds$", "\\1", basename(files))

# replace old logic with these sample names
LCM_asCN_profile_samples <- data.frame(
  sample = sample_names,
  CN_profile = "Yes"
)


### Alex Stein: Not needed here
#Get images
# LCM_images <- data.frame(region = c(rep("Primary",3), rep("T1.1",3), rep("T4.1",3)),
#                          side = rep(c("A -Front", "B -Back", "C -Side"),3),
#                          images = c(paste0(LCM.dir, "data/images/Primary/", c("Capture_12-24POSITIONS_WRITTEN.JPG", "Primary_B_OV20X_cut_positions.jpg", "Primary_C_side_OV20x_CUT_positions.jpg")),
#                                     paste0(LCM.dir, "data/images/R1/", c("OVERVIEW_1.25X_CUT_24.JPG", "T1.1_B_Overview_cut_positions.jpg", "T1.1_side_OV_20x_CUT_positions.jpg")),
#                                     paste0(LCM.dir, "data/images/R4/", c("T4.1_Front_OVERVIEW_CUT_6.3x_positions.jpeg", "T4.1_Back_OV_20x_cut.jpg", "T4.1_side_OV_cut_20x.jpg"))), stringsAsFactors = F)
# 
# LCM_image_jpg <- lapply(LCM_images$images, function(j) {
#   readJPEG(j)
# })
# names(LCM_image_jpg) <- unlist(lapply(regions, function(s) paste0(s, collapse = "_")))

#Note originally tried cristina's coordinates but remapped manually by setting x y ranges to 0 to 100 and clicking on spots


LCM_coordinates <- data.frame(barcode = gsub("\\)", "", gsub(".*\\(", "", LCM_locations$UDF.Index)),
                              coord.x = gsub(",.*", "", LCM_locations$Coordinates..X.Y..mm) %>% as.numeric(),
                              coord.y = gsub(".*,", "", LCM_locations$Coordinates..X.Y..mm) %>% as.numeric(),
                              region = LCM_locations$Region,
                              side = LCM_locations$Side,
                              sample_number = LCM_locations$Sample.Name.LCM.area) %>% #filter(barcode != "") %>% 
  left_join(LCM_barcodes, by = c("barcode", "region")) %>%
  left_join(LCM_asCN_profile_samples, by = "sample") %>% mutate(CN_profile = ifelse(is.na(CN_profile), "No", "Yes"))

### Alex Stein: Problem, we need to fix LCM_coordinates first

#Load in genotype SNP data
both_haplo_count_combined_by_cell <- readRDS(paste0(LCM.dir, "/Genotype_SNP/both_haplo_count_combined_by_cell.rds"))

both_haplo_count_combined_by_cell <- both_haplo_count_combined_by_cell %>% mutate(Haplotype_Ratio = Haplo_1_Count/Haplo_2_Count)
LCM_coordinates_tumour <- LCM_coordinates %>% left_join(both_haplo_count_combined_by_cell, by = c("barcode", "coord.x", "coord.y", "region", "side", "sample_number", "CN_profile", "sample" = "Barcode")) %>%
  filter(Haplotype_Ratio > 2) %>% filter(sample !="VER315A36_S270_L006") %>% filter(sample !="VER315A110_S319_L007")




