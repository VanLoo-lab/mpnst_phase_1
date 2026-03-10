#Script containing common metadata for Visium
library(tidyverse)

visium_samples = c("Pa", "R1", "R3", "Pb", "R2", "R4", "R5a", "R5b")
names(visium_samples) = c("VER683A1", "VER683A2", "VER683A3", "VER683A4", "VER683A5", "VER683A6", "VER683A7", "VER683A8")
# sample = "C"
color_scale <- colorRampPalette(c("darkblue", "white", "darkred"))(100)

combos <- combn(visium_samples[c(2:8,1)], 2)
pairs <- paste0(combos[1,], "_", combos[2,])

#Load bulk CN profiles
bb.input.dir = paste0(INPUTDIR,"bulk/BB_spiked/All_with_Ext/")
bulk_samples <- c("R1", "R2", "R3", "R4", "R5", "P")
bb_files_orig <- paste0(bb.input.dir,bulk_samples,"_subclones.txt")
names(bb_files_orig) <- bulk_samples
bb_subclones <- lapply(bulk_samples, function(s) {read.delim(bb_files_orig[s])})
names(bb_subclones) <- bulk_samples

#raw.visium.input.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_spatial/data/run1_2_3_man_align/")
#input.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_spatial/results/separate/")
#r_p_files <- paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/bulk/results/BB/All_with_Ext/",bulk_samples,"_rho_and_psi.txt")
#gencode <- read_delim("/camp/project/proj-vanloo/analyses/hyan/ref_files/inferCNV/gencode_v21_gen_pos.complete.clean.txt", delim = "\t", col_names = c("gene", "chr", "start", "end"))

