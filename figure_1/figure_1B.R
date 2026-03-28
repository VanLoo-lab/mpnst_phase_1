### title: plotting subclonal CN profiles from the Battenberg algorithm

# Comment: The figure was generated based on SNP counts, 
#   which must be treated as sensitive data. 
#   --> Part 1a generates the figure with SNP count data file
#   --> Part 1b generates the figure without SNP count data file

# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"

# Define input directory with data and output directory to save the figure
input.dir <- paste0(zenodo.dir, "data/bulk/battenberg/")
output.dir <- paste0(zenodo.dir, "results/figure_1/")
setwd(output.dir)

####################################################################################################################################
### Part 0: Load libraries and data
####################################################################################################################################

### Load libraries
library(tidyverse)

### Load data

# Copy number states
prefix <- "R1_R2_R3_R4_R5_P_spiked_BattenbergProfile"

pos_min <- readRDS(paste0(input.dir, prefix, "_pos_min.rds"))
pos_max <- readRDS(paste0(input.dir, prefix, "_pos_max.rds"))

segment_states_min <- readRDS(paste0(input.dir, prefix, "_segment_states_min.rds"))
segment_states_tot <- readRDS(paste0(input.dir, prefix, "_segment_states_tot.rds"))


####################################################################################################################################
### Part 1a: Create Figure 1B (if BAFvals_all is given)
####################################################################################################################################

if(F) {
  ## Load BAF values data
  BAFvals_all <- readRDS(paste0(input.dir, "sensitive/BAFvals_all.rds"))
  
  ### Create chromosome labels
  chr_names = c(1:22,"X")
  chr.segs = lapply(1:length(chr_names), function(ch) { which(BAFvals_all[[1]]$Chromosome==chr_names[ch]) })
  
  ### Define maximum ploidy status to be displayed
  ylim <- 8
  
  ### Figure generation
  
  # Define pdf file to save the figure
  pdf(file = "figure_1B.pdf", width = 14, height = 4)
  
  # Plot Battenberg profile
  par(mar = c(3,5,5,0.5), cex = 0.8, cex.main=2, cex.axis = 2, lwd = 1.5)
  maintitle = paste0("Battenberg profile across regions")
  plot(c(20000,nrow(BAFvals_all[[1]])-20000), c(0,ylim), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
  abline(v=0,lty=1,col="lightgrey")
  # Horizontal lines for y=0 to y=5
  abline(h=c(0:ylim),lty=1,col="lightgrey")
  #Plot CN for each region with unique colours
  region_colours <- c("#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3", "#F8766D")
  for (r in c((length(BAFvals_all)-1):1,length(BAFvals_all))) {
    # Minor allele in gray, total CN in orange "#E69F00" (original colour)
    segments(x0=pos_min[[r]], y0=segment_states_min[[r]], x1=pos_max[[r]], y1=segment_states_min[[r]], col="#2f4f4f", pch="|", lwd=5, lend=1)
    segments(x0=pos_min[[r]], y0=segment_states_tot[[r]], x1=pos_max[[r]], y1=segment_states_tot[[r]], col=region_colours[r], pch="|", lwd=5, lend=1)
  }
  
  # Plot vertical lines showing start/end of chromosomes
  chrk_tot_len = 0
  for (i in 1:length(chr.segs)) {
    chrk = chr.segs[[i]];
    chrk_hetero = names(BAFvals_all[[1]])[chrk]
    chrk_tot_len_prev = chrk_tot_len
    chrk_tot_len = chrk_tot_len + length(chrk_hetero)
    vpos = chrk_tot_len;
    tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
    abline(v=vpos,lty=1,col="lightgrey")
    text(tpos,ylim-0.4*(i %% 2),chr_names[i], pos = 1, cex = 1.5)
  }
  
  # Add the legend
  legend("bottom", inset=c(0,-0.16), xpd = TRUE, legend=c("P", paste0("R", 1:5)), horiz = T, bty = "n",
         col=region_colours[c(6,1:5)], bg = "white", pch=16, cex=1.5)
  
  dev.off()
}

####################################################################################################################################
### Part 1b: Create Figure 1B (if BAFvals_all is *not* given)
####################################################################################################################################

if(T) {
  ### Define metadata
  n_SNPs <- 2009404
  n_samples <- 6
  
  chr_names = c(1:22,"X")
  chr.segs <- readRDS(paste0(input.dir,"chr_segs.rds"))
  chrk_hetero <- readRDS(paste0(input.dir,"chrk_hetero.rds"))
  
  ### Define maximum ploidy status to be displayed
  ylim <- 8
  
  ### Figure generation
  
  # Define pdf file to save the figure
  pdf(file = "figure_1B.pdf", width = 14, height = 4)
  
  # Plot Battenberg profile
  par(mar = c(3,5,5,0.5), cex = 0.8, cex.main=2, cex.axis = 2, lwd = 1.5)
  maintitle = paste0("Battenberg profile across regions")
  plot(c(20000,n_SNPs-20000), c(0,ylim), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
  abline(v=0,lty=1,col="lightgrey")
  # Horizontal lines for y=0 to y=5
  abline(h=c(0:ylim),lty=1,col="lightgrey")
  #Plot CN for each region with unique colours
  region_colours <- c("#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3", "#F8766D")
  for (r in c((n_samples-1):1,n_samples)) {
    # Minor allele in gray, total CN in orange "#E69F00" (original colour)
    segments(x0=pos_min[[r]], y0=segment_states_min[[r]], x1=pos_max[[r]], y1=segment_states_min[[r]], col="#2f4f4f", pch="|", lwd=5, lend=1)
    segments(x0=pos_min[[r]], y0=segment_states_tot[[r]], x1=pos_max[[r]], y1=segment_states_tot[[r]], col=region_colours[r], pch="|", lwd=5, lend=1)
  }
  
  # Plot vertical lines showing start/end of chromosomes
  chrk_tot_len = 0
  chrk_hetero_lengths <- c(151991, 173568, 148854, 150299, 138612, 131345, 115305, 115952,
           85144, 100690, 101462, 94704, 75892, 65222, 58190, 57598,
           46853, 60618, 35497, 46195, 28414, 21177, 5822)
  for (i in 1:length(chr.segs)) {
    chrk = chr.segs[[i]];
    chrk_hetero_length = chrk_hetero_lengths[i]
    chrk_tot_len_prev = chrk_tot_len
    chrk_tot_len = chrk_tot_len + chrk_hetero_length
    vpos = chrk_tot_len;
    tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
    abline(v=vpos,lty=1,col="lightgrey")
    text(tpos,ylim-0.4*(i %% 2),chr_names[i], pos = 1, cex = 1.5)
  }
  
  # Add the legend
  legend("bottom", inset=c(0,-0.16), xpd = TRUE, legend=c("P", paste0("R", 1:5)), horiz = T, bty = "n",
         col=region_colours[c(6,1:5)], bg = "white", pch=16, cex=1.5)
  
  dev.off()
}


