### title: plotting subclonal CN profiles from the Battenberg algorithm
# the code generates Figure 1B in the manuscript (as well as other outputs)
# this file was previously named "BB_scDNA_analysis.R"

### To run the code:
# (1) Make sure to load SessionInfo() for correct package version 
# (2) Set both INPUTDIR and OUTPUTDIR to correct path
# (3) Part 4 contains essential code to recreate Figure 1B.
#     Parts can be ignored by setting False in the "if" commands, i.e. "if(F)".

# Computational cost is not too high such that the code should be feasible to 
#   run on local machines.

####################################################################################################################################
### Part 0: Code preperation
####################################################################################################################################

# In case plots look weird, the following setting may fix it
#options(bitmapType='cairo') 

### Libraries
library(tidyverse)
library(GenomicRanges)
library(VariantAnnotation)
library(qqman)
library(grid)
library(gridExtra)
library(pheatmap)
library(pbapply)
library(Battenberg)
library(ASCAT)
library(doParallel)
library(foreach)


TUMOURNAME = paste0("VER236A",1:6)
names(TUMOURNAME) <- c(paste0("R",1:5),"P")

INPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_1/data/"
OUTPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_1/results/figure_1b/"
setwd(OUTPUTDIR)

BB.dir        <- paste0(INPUTDIR, "bulk/BB/All_with_Ext/", sep="")
BB_spiked.dir <- paste0(INPUTDIR, "bulk/BB_spiked/All_with_Ext/", sep="")
scDNA.dir     <- paste0(INPUTDIR, "10X/pcf/fittedCN/", sep="")

n = 12
gg_colours <- hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]
# barplot(1:n, col = gg_colours)

####################################################################################################################################
### Part 1: Originally run with normal cells included (in folder with_normal) then rerun by excluding normal cells
####################################################################################################################################

#Split fitted CN profiles by region
if (T) {
  MPNST_scDNA_CN_mtx <- readRDS(paste0(scDNA.dir, "MPNST_scDNA_CN_mtx.rds"))
  ##Second run remove normal cells
  merge_cluster_ids <- readRDS(paste0(scDNA.dir, "MPNST_scDNA_CN_man_ward_ids_125.rds"))
  MPNST_scDNA_CN_mtx <- MPNST_scDNA_CN_mtx[!rownames(MPNST_scDNA_CN_mtx) %in% unlist(merge_cluster_ids[102]),]
  
  MPNST_scDNA_CN <- lapply(names(TUMOURNAME), function(r) {
    bc_by_region <- rownames(MPNST_scDNA_CN_mtx)[grep(r, rownames(MPNST_scDNA_CN_mtx))]
    return(MPNST_scDNA_CN_mtx[bc_by_region,])
  })
  
  #Plot for each region
  system(paste0("mkdir -p ", scDNA.dir, "by_region"))
  # system(paste0("mkdir -p ", scDNA.dir, "by_region_with_normal"))
  
  chr_probes <- readRDS(paste0(scDNA.dir, "../chr_probes.rds"))
  scDNA_CN_colour <- c("#00008B", "#7B7BC3", "#FFFFFF", "#FFCCCC", "#FF8080", "#FF3333", "#E60000", "#B31000", "#8B0000", "#660000", "#330000")
  scDNA_CN_breaks = seq(-0.5,10.5, by = 1)
  
  for (i in 1:length(TUMOURNAME)) {
    print("Clustering")
    saveRDS(hclust(dist(MPNST_scDNA_CN[[i]], method = "manhattan"), method = "ward.D2"), 
            paste0(scDNA.dir, "by_region/MPNST_scDNA_CN_",names(TUMOURNAME)[i],"_hclust_man_ward.rds"))#Use manhattan distance
    MPNST_CN_region_hclust_man_ward <- readRDS(paste0(scDNA.dir, "by_region/MPNST_scDNA_CN_",names(TUMOURNAME)[i],"_hclust_man_ward.rds"))
  
    #Plot CN profiles per region
    png(filename = paste0(scDNA.dir, "by_region/MPNST_scDNA_CN_",names(TUMOURNAME)[i],".png"), width = 4000, height = 2000, res = 200)
    pheatmap(mat = MPNST_scDNA_CN[[i]], cluster_rows = MPNST_CN_region_hclust_man_ward, cluster_cols = F, 
             show_rownames = F, show_colnames = F, color = scDNA_CN_colour, breaks = scDNA_CN_breaks, fontsize = 14, main = paste0("MPNST ", names(TUMOURNAME)[i], " scDNA Fitted CN Heatmap"))
    grid.lines(x=chr_probes$cum.probes[1]/ncol(MPNST_scDNA_CN[[i]])*0.933+0.038, y=c(0.004,0.984), gp=gpar(col="black", lwd=2))
    for (k in chr_probes$chrom[-1]) {
      grid.lines(x=chr_probes$cum.probes[k+1]/ncol(MPNST_scDNA_CN[[i]])*0.933+0.038, y=c(0.004,0.984), gp=gpar(col="black", lwd=2))
      grid.text(chr_probes$chrom[k+1], x=chr_probes$cum.probes[k]/ncol(MPNST_scDNA_CN[[i]])*0.933+0.048, y=0.01, gp=gpar(cex = 2))
    }
    dev.off()
  }
}

#Generate freq table for top 5 CN integers
if (T) {
  for (i in 1:length(TUMOURNAME)) {
    MPNST_scDNA_sc <- lapply(1:ncol(MPNST_scDNA_CN[[i]]), function(b) {
      CN_freq <- MPNST_scDNA_CN[[i]][,b] %>% table %>% as_tibble %>% rename('.' = "CN") %>% arrange(-n) %>% dplyr::slice(1:5) %>% 
        add_row(CN = "rest", n = nrow(MPNST_scDNA_CN[[i]])-sum(.$n))
      if (length(CN_freq$n) == 6) {
        return(CN_freq$n)
      } else {
        return(c(CN_freq$n, rep(0, 6-length(CN_freq$n))))
      }
    })
    MPNST_scDNA_sc <- do.call(cbind, MPNST_scDNA_sc)
    saveRDS(MPNST_scDNA_sc, paste0("MPNST_scDNA_",names(TUMOURNAME)[i],"_sc_freq.rds"))
    
    MPNST_scDNA_sc <- cbind(c(1:6), MPNST_scDNA_sc) %>% as_tibble
    colnames(MPNST_scDNA_sc) <- c("rank",1:ncol(MPNST_scDNA_sc))
    
    png(filename = paste0("MPNST_",names(TUMOURNAME)[[i]],"_sc_freq.png"), width = 4000, height = 2000, res = 200)
    print(MPNST_scDNA_sc %>% pivot_longer(-rank, names_to = "position", values_to = "count") %>% mutate(position = as.integer(position)) %>%
            ggplot(aes(fill = as.factor(rank), y = count, x = position)) + geom_bar(position="fill", stat="identity", width = 1) +
            labs(title = paste0("Subclones: MPNST ", names(TUMOURNAME)[[i]])) +
            theme_minimal() +
            theme(legend.position = "none",
                  plot.title = element_text(size = 20),
                  axis.text.x=element_blank(),
                  axis.title.x=element_blank()))
    grid.lines(x=chr_probes$cum.probes[1]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.072, y=c(0.004,0.917), gp=gpar(col="black", lwd=2))
    for (k in chr_probes$chrom[-1]) {
      grid.lines(x=chr_probes$cum.probes[k+1]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.072, y=c(0.004,0.917), gp=gpar(col="black", lwd=2))
      grid.text(chr_probes$chrom[k+1], x=chr_probes$cum.probes[k]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.082, y=0.01, gp=gpar(cex = 2))
    }
    dev.off()
  }
}

#Generate subclone % across genome
if (T) {
  BB_sc <- list()
  for (i in 1:length(TUMOURNAME)) {
    BB_sc[[i]] <- read_tsv(paste0(BB.dir, names(TUMOURNAME)[[i]], "_subclones.txt"))
  }
  
  BB_spiked_sc <- list()
  for (i in 1:length(TUMOURNAME)) {
    BB_spiked_sc[[i]] <- read_tsv(paste0(BB_spiked.dir, names(TUMOURNAME)[[i]], "_subclones.txt"))
  }
  
  chr_probes <- readRDS(paste0(scDNA.dir, "../chr_probes.rds"))
  chr_probes_positions <- lapply(2:nrow(chr_probes), function(r) {
    return((1:chr_probes$total.probes[r] * 500000) - 250000)
  })
  chr_probes_positions <- do.call(rbind, lapply(1:length(chr_probes_positions), function(c) {
    do.call(rbind, lapply(chr_probes_positions[[c]], function (r) {data.frame("CHR" = c, "POS" = r)}))
  }))
  
  BB_sc_frac <- list()
  for (i in 1:length(TUMOURNAME)) {
    #Load in scDNA subclone data
    MPNST_scDNA_sc <- readRDS(paste0("MPNST_scDNA_",names(TUMOURNAME)[i],"_sc_freq.rds"))
    MPNST_scDNA_sc <- cbind(c(1:6), MPNST_scDNA_sc) %>% as_tibble
    colnames(MPNST_scDNA_sc) <- c("rank",1:ncol(MPNST_scDNA_sc))
    
    #Get frac1 and frac2 from BB calls for each position
    BB_sc_frac[[i]] <- lapply(1:nrow(chr_probes_positions), function(p) {
      pos_frac <- BB_sc[[i]] %>% filter(chr == chr_probes_positions[p,1] & 
                                          startpos < chr_probes_positions[p,2] &
                                          endpos > chr_probes_positions[p,2]) %>%
        dplyr::select(frac1_A, frac2_A)
      return(as.data.frame(pos_frac[1,]))
    })
    BB_sc_frac[[i]] <- do.call(rbind, BB_sc_frac[[i]])
    #Remove NAs and plot larger subclone
    BB_sc_frac[[i]]$frac1_A[is.na(BB_sc_frac[[i]]$frac1_A)] <- 0
    BB_sc_frac[[i]] <- BB_sc_frac[[i]] %>% mutate(highest_frac = ifelse(frac1_A>0.5, 1-frac1_A,frac1_A),
                                                  position = as.integer(rownames(.)))
    #Plot subclone percent across genome
    png(filename = paste0("MPNST_",names(TUMOURNAME)[[i]],"_sc_freq_BB.png"), width = 4000, height = 2000, res = 200)
    print(MPNST_scDNA_sc %>% pivot_longer(-rank, names_to = "position", values_to = "count") %>% mutate(position = as.integer(position)) %>%
            left_join(BB_sc_frac[[i]], by = "position") %>%
            ggplot(aes(fill = as.factor(rank), y = count, x = position)) + geom_bar(position="fill", stat="identity", width = 1) +
            geom_point(aes(x = position, y = highest_frac, fill = "black")) +
            labs(title = paste0("Subclones: MPNST ", names(TUMOURNAME)[[i]])) +
            theme_minimal() +
            theme(legend.position = "none",
                  plot.title = element_text(size = 20),
                  axis.text.x=element_blank(),
                  axis.title.x=element_blank()))
    grid.lines(x=chr_probes$cum.probes[1]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.072, y=c(0.004,0.917), gp=gpar(col="black", lwd=2))
    for (k in chr_probes$chrom[-1]) {
      grid.lines(x=chr_probes$cum.probes[k+1]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.072, y=c(0.004,0.917), gp=gpar(col="black", lwd=2))
      grid.text(chr_probes$chrom[k+1], x=chr_probes$cum.probes[k]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.082, y=0.01, gp=gpar(cex = 2))
    }
    dev.off()
  }
  saveRDS(BB_sc_frac, "BB_sc_frac.rds")
  
  BB_spiked_sc_frac <- list()
  for (i in 1:length(TUMOURNAME)) {
    #Load in scDNA subclone data
    MPNST_scDNA_sc <- readRDS(paste0("MPNST_scDNA_",names(TUMOURNAME)[i],"_sc_freq.rds"))
    MPNST_scDNA_sc <- cbind(c(1:6), MPNST_scDNA_sc) %>% as_tibble
    colnames(MPNST_scDNA_sc) <- c("rank",1:ncol(MPNST_scDNA_sc))
    
    #Get frac1 and frac2 from BB calls for each position
    BB_spiked_sc_frac[[i]] <- lapply(1:nrow(chr_probes_positions), function(p) {
      pos_frac <- BB_spiked_sc[[i]] %>% filter(chr == chr_probes_positions[p,1] & 
                                          startpos < chr_probes_positions[p,2] &
                                          endpos > chr_probes_positions[p,2]) %>%
        dplyr::select(frac1_A, frac2_A)
      return(as.data.frame(pos_frac[1,]))
    })
    BB_spiked_sc_frac[[i]] <- do.call(rbind, BB_spiked_sc_frac[[i]])
    #Remove NAs and plot larger subclone
    BB_spiked_sc_frac[[i]]$frac1_A[is.na(BB_spiked_sc_frac[[i]]$frac1_A)] <- 0
    BB_spiked_sc_frac[[i]] <- BB_spiked_sc_frac[[i]] %>% mutate(highest_frac = ifelse(frac1_A>0.5, 1-frac1_A,frac1_A),
                                                  position = as.integer(rownames(.)))
    #Plot subclone percent across genome
    png(filename = paste0("MPNST_",names(TUMOURNAME)[[i]],"_sc_freq_BB_spiked.png"), width = 4000, height = 2000, res = 200)
    print(MPNST_scDNA_sc %>% pivot_longer(-rank, names_to = "position", values_to = "count") %>% mutate(position = as.integer(position)) %>%
            left_join(BB_spiked_sc_frac[[i]], by = "position") %>%
            ggplot(aes(fill = as.factor(rank), y = count, x = position)) + geom_bar(position="fill", stat="identity", width = 1) +
            geom_point(aes(x = position, y = highest_frac, fill = "black")) +
            labs(title = paste0("Subclones: MPNST ", names(TUMOURNAME)[[i]])) +
            theme_minimal() +
            theme(legend.position = "none",
                  plot.title = element_text(size = 20),
                  axis.text.x=element_blank(),
                  axis.title.x=element_blank()))
    grid.lines(x=chr_probes$cum.probes[1]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.072, y=c(0.004,0.917), gp=gpar(col="black", lwd=2))
    for (k in chr_probes$chrom[-1]) {
      grid.lines(x=chr_probes$cum.probes[k+1]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.072, y=c(0.004,0.917), gp=gpar(col="black", lwd=2))
      grid.text(chr_probes$chrom[k+1], x=chr_probes$cum.probes[k]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.082, y=0.01, gp=gpar(cex = 2))
    }
    dev.off()
  }
  saveRDS(BB_spiked_sc_frac, "BB_spiked_sc_frac.rds")
}

####################################################################################################################################
### Part 2: Run with CN values rather than just top CNs
####################################################################################################################################
#Generate freq table for CN 0:10 and rest
if (T) {
  for (i in 1:length(TUMOURNAME)) {
    MPNST_scDNA_sc_cn <- lapply(1:ncol(MPNST_scDNA_CN[[i]]), function(b) {
      
      CN_freq <- MPNST_scDNA_CN[[i]][,b] %>% as_tibble %>% mutate(CN = as.character(value)) %>% group_by(CN) %>% summarise(n = n()) %>% 
        right_join(tibble(CN = as.character(0:10)), by = "CN") %>% replace_na(list(n = 0)) %>% add_row(CN = "rest", n = nrow(MPNST_scDNA_CN[[i]])-sum(.$n))
      return(CN_freq$n)
      # if (length(CN_freq$n) == 6) {
      #   return(CN_freq$n)
      # } else {
      #   return(c(CN_freq$n, rep(0, 6-length(CN_freq$n))))
      # }
    })
    MPNST_scDNA_sc_cn <- do.call(cbind, MPNST_scDNA_sc_cn)
    saveRDS(MPNST_scDNA_sc_cn, paste0("MPNST_scDNA_",names(TUMOURNAME)[i],"_sc_cn_freq.rds"))
    
    MPNST_scDNA_sc_cn <- cbind(c(0:11), MPNST_scDNA_sc_cn) %>% as_tibble
    colnames(MPNST_scDNA_sc_cn) <- c("rank",1:ncol(MPNST_scDNA_sc_cn))
    
    png(filename = paste0("MPNST_",names(TUMOURNAME)[[i]],"_sc_cn_freq.png"), width = 4000, height = 2000, res = 200)
    print(MPNST_scDNA_sc_cn %>% pivot_longer(-rank, names_to = "position", values_to = "count") %>% mutate(position = as.integer(position)) %>%
            ggplot(aes(fill = as.factor(rank), y = count, x = position)) + geom_bar(position="fill", stat="identity", width = 1) +
            scale_fill_manual(name = "CN", values = gg_colours[c(1,5,9,2,6,10,3,7,11,4,8,12)]) +
            labs(title = paste0("Subclones: MPNST ", names(TUMOURNAME)[[i]])) +
            theme_minimal() +
            theme(legend.position = c(0.97, 0.783),
                  legend.box.background = element_blank(),
                  legend.background = element_rect(fill = "white"),
                  plot.title = element_text(size = 20),
                  axis.text.x=element_blank(),
                  axis.title.x=element_blank()))
    grid.lines(x=chr_probes$cum.probes[1]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.072, y=c(0.004,0.917), gp=gpar(col="black", lwd=2))
    for (k in chr_probes$chrom[-1]) {
      grid.lines(x=chr_probes$cum.probes[k+1]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.072, y=c(0.004,0.917), gp=gpar(col="black", lwd=2))
      grid.text(chr_probes$chrom[k+1], x=chr_probes$cum.probes[k]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.082, y=0.01, gp=gpar(cex = 2))
    }
    dev.off()
  }
}

#Generate subclone % across genome
if (T) {
  BB_sc <- list()
  for (i in 1:length(TUMOURNAME)) {
    BB_sc[[i]] <- read_tsv(paste0(BB.dir, names(TUMOURNAME)[[i]], "_subclones.txt"))
  }
  
  BB_spiked_sc <- list()
  for (i in 1:length(TUMOURNAME)) {
    BB_spiked_sc[[i]] <- read_tsv(paste0(BB_spiked.dir, names(TUMOURNAME)[[i]], "_subclones.txt"))
  }
  
  chr_probes <- readRDS(paste0(scDNA.dir, "../chr_probes.rds"))
  chr_probes_positions <- lapply(2:nrow(chr_probes), function(r) {
    return((1:chr_probes$total.probes[r] * 500000) - 250000)
  })
  chr_probes_positions <- do.call(rbind, lapply(1:length(chr_probes_positions), function(c) {
    do.call(rbind, lapply(chr_probes_positions[[c]], function (r) {data.frame("CHR" = c, "POS" = r)}))
  }))
  
  BB_sc_frac <- list()
  for (i in 1:length(TUMOURNAME)) {
    #Load in scDNA subclone data
    MPNST_scDNA_sc_cn <- readRDS(paste0("MPNST_scDNA_",names(TUMOURNAME)[i],"_sc_cn_freq.rds"))
    MPNST_scDNA_sc_cn <- cbind(c(0:11), MPNST_scDNA_sc_cn) %>% as_tibble
    colnames(MPNST_scDNA_sc_cn) <- c("rank",1:ncol(MPNST_scDNA_sc_cn))

    #Get frac1 and frac2 from BB calls for each position
    BB_sc_frac[[i]] <- lapply(1:nrow(chr_probes_positions), function(p) {
      pos_frac <- BB_sc[[i]] %>% filter(chr == chr_probes_positions[p,1] & 
                                          startpos < chr_probes_positions[p,2] &
                                          endpos > chr_probes_positions[p,2]) %>%
        dplyr::select(frac1_A, frac2_A)
      return(as.data.frame(pos_frac[1,]))
    })
    BB_sc_frac[[i]] <- do.call(rbind, BB_sc_frac[[i]])
    #Remove NAs and plot larger subclone
    BB_sc_frac[[i]]$frac1_A[is.na(BB_sc_frac[[i]]$frac1_A)] <- 0
    BB_sc_frac[[i]] <- BB_sc_frac[[i]] %>% mutate(highest_frac = ifelse(frac1_A>0.5, 1-frac1_A,frac1_A),
                                                  position = as.integer(rownames(.)))
    #Plot subclone percent across genome
    png(filename = paste0("MPNST_",names(TUMOURNAME)[[i]],"_sc_cn_freq_BB.png"), width = 4000, height = 2000, res = 200)
    print(MPNST_scDNA_sc_cn %>% pivot_longer(-rank, names_to = "position", values_to = "count") %>% mutate(position = as.integer(position)) %>%
            left_join(BB_sc_frac[[i]], by = "position") %>%
            ggplot(aes(fill = as.factor(rank), y = count, x = position)) + geom_bar(position="fill", stat="identity", width = 1) +
            geom_point(aes(x = position, y = highest_frac)) +
            scale_fill_manual(name = "CN", values = c(gg_colours[c(1,5,9,2,6,10,3,7,11,4,8,12)], "black")) +
            labs(title = paste0("Subclones: MPNST ", names(TUMOURNAME)[[i]])) +
            theme_minimal() +
            theme(legend.position = c(0.97, 0.783),
                  legend.box.background = element_blank(),
                  legend.background = element_rect(fill = "white"),
                  plot.title = element_text(size = 20),
                  axis.text.x=element_blank(),
                  axis.title.x=element_blank()))
    grid.lines(x=chr_probes$cum.probes[1]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.072, y=c(0.004,0.917), gp=gpar(col="black", lwd=2))
    for (k in chr_probes$chrom[-1]) {
      grid.lines(x=chr_probes$cum.probes[k+1]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.072, y=c(0.004,0.917), gp=gpar(col="black", lwd=2))
      grid.text(chr_probes$chrom[k+1], x=chr_probes$cum.probes[k]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.082, y=0.01, gp=gpar(cex = 2))
    }
    dev.off()
  }
  saveRDS(BB_sc_frac, "BB_sc_cn_frac.rds")
  
  BB_spiked_sc_frac <- list()
  for (i in 1:length(TUMOURNAME)) {
    #Load in scDNA subclone data
    MPNST_scDNA_sc_cn <- readRDS(paste0("MPNST_scDNA_",names(TUMOURNAME)[i],"_sc_cn_freq.rds"))
    MPNST_scDNA_sc_cn <- cbind(c(0:11), MPNST_scDNA_sc_cn) %>% as_tibble
    colnames(MPNST_scDNA_sc_cn) <- c("rank",1:ncol(MPNST_scDNA_sc_cn))
    
    #Get frac1 and frac2 from BB calls for each position
    BB_spiked_sc_frac[[i]] <- lapply(1:nrow(chr_probes_positions), function(p) {
      pos_frac <- BB_spiked_sc[[i]] %>% filter(chr == chr_probes_positions[p,1] & 
                                                 startpos < chr_probes_positions[p,2] &
                                                 endpos > chr_probes_positions[p,2]) %>%
        dplyr::select(frac1_A, frac2_A)
      return(as.data.frame(pos_frac[1,]))
    })
    BB_spiked_sc_frac[[i]] <- do.call(rbind, BB_spiked_sc_frac[[i]])
    #Remove NAs and plot larger subclone
    BB_spiked_sc_frac[[i]]$frac1_A[is.na(BB_spiked_sc_frac[[i]]$frac1_A)] <- 0
    BB_spiked_sc_frac[[i]] <- BB_spiked_sc_frac[[i]] %>% mutate(highest_frac = ifelse(frac1_A>0.5, 1-frac1_A,frac1_A),
                                                                position = as.integer(rownames(.)))
    #Plot subclone percent across genome
    png(filename = paste0("MPNST_",names(TUMOURNAME)[[i]],"_sc_cn_freq_BB_spiked.png"), width = 4000, height = 2000, res = 200)
    print(MPNST_scDNA_sc_cn %>% pivot_longer(-rank, names_to = "position", values_to = "count") %>% mutate(position = as.integer(position)) %>%
            left_join(BB_spiked_sc_frac[[i]], by = "position") %>%
            ggplot(aes(fill = as.factor(rank), y = count, x = position)) + geom_bar(position="fill", stat="identity", width = 1) +
            geom_point(aes(x = position, y = highest_frac)) +
            scale_fill_manual(name = "CN", values = c(gg_colours[c(1,5,9,2,6,10,3,7,11,4,8,12)], "black")) +
            labs(title = paste0("Subclones: MPNST ", names(TUMOURNAME)[[i]])) +
            theme_minimal() +
            theme(legend.position = c(0.97, 0.783),
                  legend.box.background = element_blank(),
                  legend.background = element_rect(fill = "white"),
                  plot.title = element_text(size = 20),
                  axis.text.x=element_blank(),
                  axis.title.x=element_blank()))
    grid.lines(x=chr_probes$cum.probes[1]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.072, y=c(0.004,0.917), gp=gpar(col="black", lwd=2))
    for (k in chr_probes$chrom[-1]) {
      grid.lines(x=chr_probes$cum.probes[k+1]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.072, y=c(0.004,0.917), gp=gpar(col="black", lwd=2))
      grid.text(chr_probes$chrom[k+1], x=chr_probes$cum.probes[k]/ncol(MPNST_scDNA_CN[[i]])*0.88+0.082, y=0.01, gp=gpar(cex = 2))
    }
    dev.off()
  }
  saveRDS(BB_spiked_sc_frac, "BB_spiked_sc_cn_frac.rds")
}

####################################################################################################################################
### Part 3: Plot log R of bulk and sc (chromplot)
####################################################################################################################################

if (T) {
  ##Load in sc CN data
  MPNST_scDNA_CN_mtx <- readRDS(paste0(scDNA.dir, "MPNST_scDNA_CN_mtx.rds"))
  ##Second run remove normal cells
  merge_cluster_ids <- readRDS(paste0(scDNA.dir, "MPNST_scDNA_CN_man_ward_ids_125.rds"))
  MPNST_scDNA_CN_mtx <- MPNST_scDNA_CN_mtx[!rownames(MPNST_scDNA_CN_mtx) %in% unlist(merge_cluster_ids[102]),]
  
  MPNST_scDNA_CN <- lapply(names(TUMOURNAME), function(r) {
    bc_by_region <- rownames(MPNST_scDNA_CN_mtx)[grep(r, rownames(MPNST_scDNA_CN_mtx))]
    return(MPNST_scDNA_CN_mtx[bc_by_region,])
  })
  
  chr_probes <- readRDS(paste0(scDNA.dir, "../chr_probes.rds"))
  chr_probes_positions <- lapply(2:nrow(chr_probes), function(r) {
    return((1:chr_probes$total.probes[r] * 500000) - 250000)
  })
  chr_probes_positions <- do.call(rbind, lapply(1:length(chr_probes_positions), function(c) {
    do.call(rbind, lapply(chr_probes_positions[[c]], function (r) {data.frame("CHR" = c, "POS" = r)}))
  }))
  
  ##Calculate average CN in sc data then convert to positions
  MPNST_scDNA_CN_mean <- lapply(MPNST_scDNA_CN, function(r){
    CN_mean <- apply(r, 2, mean)
    return(cbind(chr_probes_positions, CN_mean))
  })
  
  ### Now comes the Battenberg analysis
  
  BB_scDNA_logr_plot = function(samplename, subclones, logr, outputfile, purity, scCN_data, plot.bulk = F, plot.sc = F) {
    # Smooth the logR
    colnames(logr)[3] = "raw_logr"
    logr$logr_smoothed = Battenberg:::runmed_data(logr$Chromosome, logr$raw_logr, 101)
    
    # Prepare subclones data
    subclones$len = subclones$endpos/1000-subclones$startpos/1000
    subclones$total_major = Battenberg:::calc_total_cn_major(subclones)
    subclones$total_minor = Battenberg:::calc_total_cn_minor(subclones)
    subclones$total_cn = subclones$total_minor + subclones$total_major
    subclones$is_subclonal = subclones$frac1_A < 1
    subclones$is_50_50 = subclones$frac1_A >= 0.48 & subclones$frac1_A <= 0.52
    
    # Calculate psi from the data
    ploidy = Battenberg:::calc_ploidy(subclones)
    psi = Battenberg:::psit2psi(purity, ploidy)
    
    # Estimate total CN for each segment based on the logR
    logr$total_cn = NA
    logr$total_cn_psi = NA
    for (j in (1:nrow(subclones))) {
      print(j)
      sel = which(logr$Chromosome == subclones$chr[j] & logr$Position >= subclones$startpos[j] & logr$Position <= subclones$endpos[j])
      tumour_cn = Battenberg:::calculate_bb_total_cn(subclones[j,,drop=F])
      total_cn = purity*tumour_cn + 2*(1-purity)
      logr$total_cn[sel] = Battenberg:::logr2tumcn(purity, total_cn, logr$logr_smoothed[sel])
      logr$total_cn_psi[sel] = Battenberg:::logr2tumcn(purity, psi, logr$logr_smoothed[sel])
    }
    
    # Plot every 100 data point, there are too many for them all to be seen
    logr_plot = logr[seq(1, nrow(logr), 100),]
    
    # Sync the levels for chromosome so that all corresponding data ends up in the same plot
    logr_plot$Chromosome = factor(logr_plot$Chromosome, levels=gtools::mixedsort(unique(logr_plot$Chromosome)))
    subclones$Chromosome = factor(subclones$chr, levels=levels(logr_plot$Chromosome))
    
    # Set plot boundaries for x and y - take as y value the maximum between the data and the fit
    max_cn_plot_data = ceiling(quantile(logr_plot$total_cn_psi, c(.98), na.rm=T))
    max_cn_plot_fit = ceiling(quantile(unlist(lapply(1:nrow(subclones), function(i) rep(subclones$total_cn[i], subclones$len[i]))), c(.98), na.rm=T))
    max_cn_plot = ifelse(max_cn_plot_fit > max_cn_plot_data, max_cn_plot_fit, max_cn_plot_data)
    maxpos = max(logr$Position)
    
    # catch case when there is no clonal CNA called
    if (is.na(max_cn_plot) | max_cn_plot < 4) {
      max_cn_plot = 4
    }
    
    # These are the grey lines in the background
    background = data.frame(xmin=rep(0, (max_cn_plot/2)+1), 
                            xmax=rep(max(logr$Position), (max_cn_plot/2)+1),
                            ymin=seq(0, max_cn_plot, 2)+0.5, 
                            ymax=seq(0, max_cn_plot, 2)+1.5)
    
    # Calc a couple of stats for the plot title
    genome_50_50 = sum(subclones$len[subclones$is_50_50]/1000)
    prop_subclonal = round(sum(subclones$len[subclones$is_subclonal]) / sum(subclones$len), 2)
    homdel = sum(subclones$len[subclones$total_cn == 0]/1000)
    plot_title = samplename
    plot_subtitle = paste0("Purity: ", round(purity, 2), " - Ploidy: ", round(ploidy, 2), " - Hom del: ", round(homdel, 2), "Mb - Prop. subclonal: ", prop_subclonal, " - Subclonal 50/50: ", round(genome_50_50, 2), "Mb")
    
    rect_height_padding = 0.2
    
    # Build the actual plot - CNA segments are drawn separately depending on their category as categories have different colours
    p = ggplot() + 
      geom_rect(data=background, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.5) +
      geom_point(data=logr_plot, mapping=aes(x=Position, y=total_cn_psi), size=0.5) + 
      ylab("Copy Number") +
      scale_y_continuous(breaks=seq(0, max_cn_plot, 2)) +  #, limits=c(-rect_height_padding, max_cn_plot+rect_height_padding)
      # Axis ticks every 10Mb
      scale_x_continuous(breaks=seq(1, max(logr$Position), 10000000)[-1], labels=round(seq(0, maxpos, 10000000) / 1000000)[-1], expand=c(0, 0)) +
      # Don't restrict the plotting area, zoom. that way segments that go outside the limits are partially plotted still
      coord_cartesian(ylim=c(-rect_height_padding, max_cn_plot+rect_height_padding)) +
      facet_wrap(~Chromosome, ncol=2, strip.position="right") + 
      # ggtitle(plot_title) +
      ggtitle(bquote(atop(.(plot_title), atop(.(plot_subtitle), "")))) +
      theme_bw() + theme(axis.title.x=element_blank(),
                         axis.text.x=element_text(colour="black",size=16,face="plain"),
                         axis.text.y = element_text(colour="black",size=16,face="plain"), 
                         axis.title.y = element_text(colour="black",size=20,face="plain"),
                         strip.text.y = element_text(colour="black",size=20,face="plain"),
                         plot.title = element_text(colour="black",size=36,face="plain",hjust = 0.5))
    
    if (plot.bulk == T) {
      # Plot the copy number segments - some of the data.frames may be empty, so check for that first before adding to the plot
      sel = !subclones$is_subclonal
      if (any(sel)) {
        # Minor allele - Normal clonal copy number
        p = p + geom_rect(data=subclones[sel, ], mapping=aes(xmin=startpos, xmax=endpos, ymin=total_minor-rect_height_padding, ymax=total_minor+rect_height_padding), fill="#2f4f4f")
      }
      sel = subclones$is_subclonal & !subclones$is_50_50
      if (any(sel)) {
        # Minor allele - Normal subclonal copy number
        p = p + geom_rect(data=subclones[sel, ], mapping=aes(xmin=startpos, xmax=endpos, ymin=total_minor-rect_height_padding, ymax=total_minor+rect_height_padding), fill="#2f3f4f")
      }
      sel = subclones$is_subclonal & subclones$is_50_50
      if (any(sel)) {
        # Minor allele - Subclonal segments right in between two clonal states
        p = p + geom_rect(data=subclones[sel, ], mapping=aes(xmin=startpos, xmax=endpos, ymin=total_minor-rect_height_padding, ymax=total_minor+rect_height_padding), fill="#2f3f4f", colour="red")
      }
      sel = !subclones$is_subclonal
      if (any(sel)) {
        # Major allele - clonal copy number
        p = p + geom_rect(data=subclones[sel, ], mapping=aes(xmin=startpos, xmax=endpos, ymin=total_cn-rect_height_padding, ymax=total_cn+rect_height_padding), fill="#E69F00")
      }
      sel = subclones$is_subclonal
      if (any(sel)) {
        # Major allele - subclonal copy number
        p = p + geom_rect(data=subclones[sel, ], mapping=aes(xmin=startpos, xmax=endpos, ymin=total_cn-rect_height_padding, ymax=total_cn+rect_height_padding), fill="#E55300")
      }
    }
    
    if (plot.sc == T) {
      #Add in sc CN lines
      scDNA_CN_mean <- scCN_data %>% mutate(Chromosome = gsub("23", "X", as.character(CHR))) %>%
        mutate(Chromosome = factor(Chromosome, levels=levels(logr_plot$Chromosome)))
      
      p = p + geom_rect(data=scDNA_CN_mean, 
                        mapping=aes(xmin=POS-249999, xmax=POS+249999, ymin=CN_mean-rect_height_padding, ymax=CN_mean+rect_height_padding), 
                        fill="#C77CFF")
    } 
    png(outputfile, width=2000, height=1300)
    print(p)
    dev.off()
  }
  
  ##For non spiked 
  for (i in 1:length(TUMOURNAME)) {
    #Code from Battenberg::make_posthoc_plots
    samplename = names(TUMOURNAME)[i]
    logr_file = paste0(BB.dir, names(TUMOURNAME)[i], "_mutantLogR_gcCorrected.tab")
    subclones_file = paste0(BB.dir, names(TUMOURNAME)[i], "_subclones.txt")
    rho_psi_file = paste0(BB.dir, names(TUMOURNAME)[i], "_rho_and_psi.txt")
    
    logr = Battenberg::read_table_generic(logr_file)
    subclones = Battenberg::read_table_generic(subclones_file)
    rho_psi = read.table(rho_psi_file, header=T, stringsAsFactors=F)
    purity = rho_psi["FRAC_GENOME", "rho"]
    
    BB_scDNA_logr_plot(samplename, subclones, logr, paste0(samplename, "_BB_totalcn_chrom_plot.png"), purity, 
                       scCN_data = MPNST_scDNA_CN_mean[[i]], plot.bulk = T, plot.sc = F)
    BB_scDNA_logr_plot(samplename, subclones, logr, paste0(samplename, "_scDNA_totalcn_chrom_plot.png"), purity,
                       scCN_data = MPNST_scDNA_CN_mean[[i]], plot.bulk = F, plot.sc = T)
    BB_scDNA_logr_plot(samplename, subclones, logr, paste0(samplename, "_BB_scDNA_totalcn_chrom_plot.png"), purity,
                       scCN_data = MPNST_scDNA_CN_mean[[i]], plot.bulk = T, plot.sc = T)
  }
  
  ##For spiked 
  for (i in 1:length(TUMOURNAME)) {
    #Code from Battenberg::make_posthoc_plots
    samplename = names(TUMOURNAME)[i]
    logr_file = paste0(BB_spiked.dir, names(TUMOURNAME)[i], "_mutantLogR_gcCorrected.tab")
    subclones_file = paste0(BB_spiked.dir, names(TUMOURNAME)[i], "_subclones.txt")
    rho_psi_file = paste0(BB_spiked.dir, names(TUMOURNAME)[i], "_rho_and_psi.txt")
    
    logr = Battenberg::read_table_generic(logr_file)
    subclones = Battenberg::read_table_generic(subclones_file)
    rho_psi = read.table(rho_psi_file, header=T, stringsAsFactors=F)
    purity = rho_psi["FRAC_GENOME", "rho"]
    
    BB_scDNA_logr_plot(samplename, subclones, logr, paste0(samplename, "_spiked_BB_totalcn_chrom_plot.png"), purity, 
                       scCN_data = MPNST_scDNA_CN_mean[[i]], plot.bulk = T, plot.sc = F)
    BB_scDNA_logr_plot(samplename, subclones, logr, paste0(samplename, "_spiked_scDNA_totalcn_chrom_plot.png"), purity,
                       scCN_data = MPNST_scDNA_CN_mean[[i]], plot.bulk = F, plot.sc = T)
    BB_scDNA_logr_plot(samplename, subclones, logr, paste0(samplename, "_spiked_BB_scDNA_totalcn_chrom_plot.png"), purity,
                       scCN_data = MPNST_scDNA_CN_mean[[i]], plot.bulk = T, plot.sc = T)
  }
}

####################################################################################################################################
### Part 4: Plot BB CN profiles across regions and look at censensus amount across genome
####################################################################################################################################
if (T) {
  #Define functions (modified from Battenberg)
  plot.gw.subclonal.cn = function(subclones, BAFvals, output.gw.figures.prefix, chr.names) {
    #@hxy265 make outputs lists
    pos_min <- list()
    pos_max <- list()
    segment_states_min <- list()
    segment_states_maj <- list()
    segment_states_tot <- list()
    
    for (s in 1:length(subclones)) {
      # Map start and end of each segment into the BAF values. The plot uses the index of this BAF table as x-axis
      pos_min[[s]] = array(NA, nrow(subclones[[s]]))
      pos_max[[s]] = array(NA, nrow(subclones[[s]]))
      for (i in 1:nrow(subclones[[s]])) {
        segm_chr = subclones[[s]]$chr[i] == BAFvals[[s]]$Chromosome & subclones[[s]]$startpos[i] < BAFvals[[s]]$Position & subclones[[s]]$endpos[i] >= BAFvals[[s]]$Position
        pos_min[[s]][i] = min(which(segm_chr))
        pos_max[[s]][i] = max(which(segm_chr))
      }
      
      # For those segments that are subclonal, Obtain the second state.
      is_subclonal = which(subclones[[s]]$frac1_A < 1)
      subcl_min = array(NA, length(is_subclonal))
      subcl_max = array(NA, length(is_subclonal))
      for (i in 1:length(is_subclonal)) {
        segment_index = is_subclonal[i]
        segm_chr = subclones[[s]]$chr[segment_index] == BAFvals[[s]]$Chromosome & subclones[[s]]$startpos[segment_index] < BAFvals[[s]]$Position & subclones[[s]]$endpos[segment_index] >= BAFvals[[s]]$Position
        subcl_min[i] = min(which(segm_chr))
        subcl_max[i] = max(which(segm_chr))
      }
      
      # Determine whether it's the major or the minor allele that is represented by two states
      is_subclonal_maj = abs(subclones[[s]]$nMaj1_A - subclones[[s]]$nMaj2_A) > 0
      is_subclonal_min = abs(subclones[[s]]$nMin1_A - subclones[[s]]$nMin2_A) > 0
      is_subclonal_maj[is.na(is_subclonal_maj)] = F
      is_subclonal_min[is.na(is_subclonal_min)] = F
      
      # BB represents subclonal CN as a mixture of two CN states. Calculate this mixture for both minor allele and total CN.
      #segment_states_min[[s]] = subclones[[s]]$nMin1_A * ifelse(is_subclonal_min, subclones[[s]]$frac1_A, 1)  + ifelse(is_subclonal_min, subclones[[s]]$nMin2_A, 0) * ifelse(is_subclonal_min, subclones[[s]]$frac2_A, 0)
      #segment_states_tot[[s]] = (subclones[[s]]$nMaj1_A+subclones[[s]]$nMin1_A) * ifelse(is_subclonal_maj, subclones[[s]]$frac1_A, 1) + ifelse(is_subclonal_maj, subclones[[s]]$nMaj2_A+subclones[[s]]$nMin2_A, 0) * ifelse(is_subclonal_maj, subclones[[s]]$frac2_A, 0)
      
      segment_states_min[[s]] = subclones[[s]]$nMin1_A * ifelse(is_subclonal_min, subclones[[s]]$frac1_A, 1)  + ifelse(is_subclonal_min, subclones[[s]]$nMin2_A, 0) * ifelse(is_subclonal_min, subclones[[s]]$frac2_A, 0)
      segment_states_maj[[s]] = subclones[[s]]$nMaj1_A * ifelse(is_subclonal_maj, subclones[[s]]$frac1_A, 1)  + ifelse(is_subclonal_maj, subclones[[s]]$nMaj2_A, 0) * ifelse(is_subclonal_maj, subclones[[s]]$frac2_A, 0)
      segment_states_tot[[s]] = segment_states_maj[[s]] + segment_states_min[[s]]
    }
    
    # Determine which SNPs are on which chromosome, to be used as a proxy for chromosome size in the plots
    chr.segs = lapply(1:length(chr.names), function(ch) { which(BAFvals[[1]]$Chromosome==chr.names[ch]) })
    
    # Plot subclonal copy number as mixtures of two states
    ###@hxy265 added all to png name and create.bb.plot.average.all function
    saveRDS(pos_min, paste0(output.gw.figures.prefix, "_pos_min.rds"))
    saveRDS(pos_max, paste0(output.gw.figures.prefix, "_pos_max.rds"))
    saveRDS(segment_states_min, paste0(output.gw.figures.prefix, "_segment_states_min.rds"))
    saveRDS(segment_states_tot, paste0(output.gw.figures.prefix, "_segment_states_tot.rds"))

    png(filename = paste(output.gw.figures.prefix, "_average_all.png", sep=""), width = 2000, height = 500, res = 200)
    create.bb.plot.average.all(bafsegmented=BAFvals,
                              pos_min=pos_min,
                              pos_max=pos_max,
                              segment_states_min=segment_states_min,
                              segment_states_tot=segment_states_tot,
                              chr.segs=chr.segs,
                              chr.names=chr.names)
    dev.off()
  }
  
  create.bb.plot.average.all = function(bafsegmented, pos_min, pos_max, segment_states_min, segment_states_tot, chr.segs, chr.names, ylim=8) {
    # Plot main frame and title
    par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
    maintitle = paste0("Battenberg Profile ")
    plot(c(1,nrow(bafsegmented[[1]])), c(0,ylim), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    abline(v=0,lty=1,col="lightgrey")
    # Horizontal lines for y=0 to y=5
    abline(h=c(0:ylim),lty=1,col="lightgrey")
    
    #Plot CN for each region
    for (r in 1:length(bafsegmented)) {
      # Minor allele in gray, total CN in orange
      segments(x0=pos_min[[r]], y0=segment_states_min[[r]], x1=pos_max[[r]], y1=segment_states_min[[r]], col="#2f4f4f", pch="|", lwd=2, lend=1)
      segments(x0=pos_min[[r]], y0=segment_states_tot[[r]], x1=pos_max[[r]], y1=segment_states_tot[[r]], col="#E69F00", pch="|", lwd=2, lend=1)
    }
    
    # Plot the vertical lines that show start/end of a chromosome
    chrk_tot_len = 0
    for (i in 1:length(chr.segs)) {
      chrk = chr.segs[[i]];
      chrk_hetero = names(bafsegmented[[1]])[chrk]
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk_hetero)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,ylim,chr.names[i], pos = 1, cex = 2)
      abline(v=vpos,lty=1,col="lightgrey")
    }
  }
  
  ###############################
  ###Plot BB CN profiles across regions for non spiked/spiked 
  ###############################
  spiked_data = T
  if (spiked_data == T) {spiked_file = "_spiked"} else {spiked_file = ""}
  
  #For multiple regions
  samplename = names(TUMOURNAME)
  # samplename = names(TUMOURNAME)[c(1,2,3,5,6)]
  baf_segmented_file = paste0(if (spiked_data == T) {BB_spiked.dir} else {BB.dir}, samplename, ".BAFsegmented.txt")
  subclones_file = paste0(if (spiked_data == T) {BB_spiked.dir} else {BB.dir}, samplename, "_subclones.txt")
  output.gw.figures.prefix = paste0(paste0(samplename, collapse = "_"), spiked_file, "_BattenbergProfile")
  chr_names = c(1:22,"X")
  
  subclones_all <- lapply(1:length(samplename), function(r) {
    subcloneres = Battenberg::read_table_generic(subclones_file[r])
    subclones = as.data.frame(subcloneres)
    subclones[,2:ncol(subclones)] = sapply(2:ncol(subclones), function(x) { as.numeric(as.character(subclones[,x])) })
    return(subclones)
  })
  BAFvals_all <- lapply(1:length(samplename), function(r) {
    return(as.data.frame(Battenberg:::read_bafsegmented(baf_segmented_file[r])))
  })
  
  # Code from Battenberg::make_posthoc_plots
  plot.gw.subclonal.cn(subclones=subclones_all, BAFvals=BAFvals_all, output.gw.figures.prefix=output.gw.figures.prefix, chr.names=chr_names)
  
  #02/04/21 Add multicolour plot
  if (T) {
    pos_min <- readRDS(paste0(output.gw.figures.prefix, "_pos_min.rds"))
    pos_max <- readRDS(paste0(output.gw.figures.prefix, "_pos_max.rds"))
    segment_states_min <- readRDS(paste0(output.gw.figures.prefix, "_segment_states_min.rds"))
    segment_states_tot <- readRDS(paste0(output.gw.figures.prefix, "_segment_states_tot.rds"))
    chr.segs = lapply(1:length(chr_names), function(ch) { which(BAFvals_all[[1]]$Chromosome==chr_names[ch]) })
    ylim <- 8
    
    pdf(file = paste0("MPNST", spiked_file, "_average_all_multicolour.pdf"), width = 14, height = 4)
    # Plot main frame and title
    par(mar = c(3,5,5,0.5), cex = 0.8, cex.main=2, cex.axis = 2, lwd = 1.5)
    maintitle = paste0("Battenberg profile across regions")
    plot(c(20000,nrow(BAFvals_all[[1]])-20000), c(0,ylim), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    abline(v=0,lty=1,col="lightgrey")
    # Horizontal lines for y=0 to y=5
    abline(h=c(0:ylim),lty=1,col="lightgrey")
    #Plot CN for each region with unique colours
    region_colours <- c("#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3", "#F8766D")
    for (r in c((length(BAFvals_all)-1):1,length(BAFvals_all))) {
      # for (r in 1:length(BAFvals_all)) {
      # Minor allele in gray, total CN in orange "#E69F00" (original colour)
      segments(x0=pos_min[[r]], y0=segment_states_min[[r]], x1=pos_max[[r]], y1=segment_states_min[[r]], col="#2f4f4f", pch="|", lwd=5, lend=1)
      segments(x0=pos_min[[r]], y0=segment_states_tot[[r]], x1=pos_max[[r]], y1=segment_states_tot[[r]], col=region_colours[r], pch="|", lwd=5, lend=1)
    }
    # Plot the vertical lines that show start/end of a chromosome
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
    legend("bottom", inset=c(0,-0.16), xpd = TRUE, legend=c("P", paste0("R", 1:5)), horiz = T, bty = "n",
           col=region_colours[c(6,1:5)], bg = "white", pch=16, cex=1.5)
    dev.off()
  }
  
  #20/02/23 Make pdf of battenberg plots
  if (T) {
    #Load purity
    purities <- lapply(names(TUMOURNAME), function(s) {
      read.delim(paste0(BB.dir, s, "_rho_and_psi.txt")) #use non spiked
    })
    
    #Load ploidy
    ploidies <- lapply(names(TUMOURNAME), function(s) {
      read.delim(paste0(BB_spiked.dir, s, "_rho_and_psi.txt")) #use spiked
    })
    
    pos_min <- readRDS(paste0(paste0(paste0(samplename, collapse = "_"), spiked_file, "_BattenbergProfile"), "_pos_min.rds"))
    pos_max <- readRDS(paste0(paste0(paste0(samplename, collapse = "_"), spiked_file, "_BattenbergProfile"), "_pos_max.rds"))
    segment_states_min <- readRDS(paste0(paste0(paste0(samplename, collapse = "_"), spiked_file, "_BattenbergProfile"), "_segment_states_min.rds"))
    segment_states_tot <- readRDS(paste0(paste0(paste0(samplename, collapse = "_"), spiked_file, "_BattenbergProfile"), "_segment_states_tot.rds"))
    chr.segs = lapply(1:length(chr_names), function(ch) { which(BAFvals_all[[1]]$Chromosome==chr_names[ch]) })
    ylim <- 8
    
    for (s in 1:length(TUMOURNAME)) {
      subclones = subclones_all[[s]]
      BAFvals = BAFvals_all[[s]]
      rho = purities[[s]]$rho[rownames(purities[[s]])=="FRAC_GENOME"]
      ploidy = ploidies[[s]]$ploidy[rownames(ploidies[[s]])=="FRAC_GENOME"]
      output.gw.figures.prefix = paste0(names(TUMOURNAME)[s], spiked_file, "_BattenbergProfile")
      chr.names = chr_names
      tumourname = names(TUMOURNAME)[s]
      
      pdf(file = paste0("MPNST_", names(TUMOURNAME)[s], spiked_file, "_average.pdf"), width = 14, height = 3)
      par(mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main = 3, cex.axis = 2.5)
      maintitle = paste0(substring(tumourname, 40, first = T), 
                         ", Ploidy: ", sprintf("%1.2f", ploidy), ", Purity: ", 
                         sprintf("%2.0f", rho * 100), "%")
      plot(c(20000,nrow(BAFvals))+20000, c(0, ylim), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
      abline(v = 0, lty = 1, col = "lightgrey")
      abline(h = c(0:ylim), lty = 1, col = "lightgrey")
      segments(x0 = pos_min[[s]], y0 = segment_states_min[[s]], x1 = pos_max[[s]], 
               y1 = segment_states_min[[s]], col = "#2f4f4f", pch = "|", 
               lwd = 6, lend = 1)
      segments(x0 = pos_min[[s]], y0 = segment_states_tot[[s]], x1 = pos_max[[s]], 
               y1 = segment_states_tot[[s]], col = "#E69F00", pch = "|", 
               lwd = 6, lend = 1)
      chrk_tot_len = 0
      for (i in 1:length(chr.segs)) {
        chrk = chr.segs[[i]]
        chrk_hetero = names(BAFvals)[chrk]
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk_hetero)
        vpos = chrk_tot_len
        tpos = (chrk_tot_len + chrk_tot_len_prev)/2
        text(tpos, ylim, chr.names[i], pos = 1, cex = 2)
        abline(v = vpos, lty = 1, col = "lightgrey")
      }
      dev.off()
    }
      
  }
  
  ###############################
  ###Look at consensus amount for spiked/non spiked
  ###############################
  output.gw.figures.prefix = paste0(paste0(samplename, collapse = "_"), spiked_file, "_BattenbergProfile")
  pos_min <- readRDS(paste0(output.gw.figures.prefix, "_pos_min.rds"))
  pos_max <- readRDS(paste0(output.gw.figures.prefix, "_pos_max.rds"))
  segment_states_min <- readRDS(paste0(output.gw.figures.prefix, "_segment_states_min.rds"))
  segment_states_tot <- readRDS(paste0(output.gw.figures.prefix, "_segment_states_tot.rds"))
  
  all_CN <- lapply(1:length(samplename), function(r) {
    CN_df <- tibble(pos_min = pos_min[[r]],
                    pos_max = pos_max[[r]],
                    segment_states_min = segment_states_min[[r]],
                    segment_states_tot = segment_states_tot[[r]])
    return(CN_df)
  })
  
  all_CN_min <- do.call(rbind, lapply(1:length(samplename), function(s) {
    return(unlist(lapply(1:nrow(all_CN[[s]]), function(r) {
      unlist(rep(all_CN[[s]][r,3], 
                 (all_CN[[s]][r,2]-all_CN[[s]][r,1]+2)))
    })))
  }))
  rownames(all_CN_min) <- names(samplename)
  
  all_CN_tot <- do.call(rbind, lapply(1:length(samplename), function(s) {
    return(unlist(lapply(1:nrow(all_CN[[s]]), function(r) {
      unlist(rep(all_CN[[s]][r,4], 
                 (all_CN[[s]][r,2]-all_CN[[s]][r,1]+2)))
    })))
  }))
  rownames(all_CN_tot) <- names(samplename)
  
  consensus_CN_min <- pbapply(all_CN_min, 2, function(p) {
    #Return freq of commonest CN
    return(sort(table(p), decreasing = T)[[1]])
  })
  saveRDS(consensus_CN_min, paste0(output.gw.figures.prefix, spiked_file, "_consensus_CN_min.rds"))
  
  consensus_CN_tot <- pbapply(all_CN_tot, 2, function(p) {
    #Return freq (not CN) of commonest CN
    return(sort(table(p), decreasing = T)[[1]])
  })
  saveRDS(consensus_CN_tot, paste0(output.gw.figures.prefix, spiked_file, "_consensus_CN_tot.rds"))
  
  consensus_CN_min_rle <- rle(consensus_CN_min)
  CN_min_pos <- data.frame(length = consensus_CN_min_rle[[1]], CN = consensus_CN_min_rle[[2]])
  CN_min_pos$sum <- cumsum(CN_min_pos$length)
  CN_min_pos$min <- CN_min_pos$sum - CN_min_pos$length+1
  CN_min_pos$max <- CN_min_pos$min + CN_min_pos$length-1
  
  consensus_CN_tot_rle <- rle(consensus_CN_tot)
  CN_tot_pos <- data.frame(length = consensus_CN_tot_rle[[1]], CN = consensus_CN_tot_rle[[2]])
  CN_tot_pos$sum <- cumsum(CN_tot_pos$length)
  CN_tot_pos$min <- CN_tot_pos$sum - CN_tot_pos$length+1
  CN_tot_pos$max <- CN_tot_pos$min + CN_tot_pos$length-1
  
  png(filename = paste(output.gw.figures.prefix, "_common_CN.png", sep=""), width = 2000, height = 500, res = 200)
  # Plot main frame and title
  par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
  maintitle = paste0("Max number of regions with shared CN")
  plot(c(1,nrow(BAFvals_all[[1]])), c(0,6.5), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
  abline(v=0,lty=1,col="lightgrey")
  # Horizontal lines for y=0 to y=5
  abline(h=c(0:6.5),lty=1,col="lightgrey")
  # Minor allele in gray, total CN in orange
  segments(x0=CN_min_pos$min, y0=CN_min_pos$CN+0.05, x1=CN_min_pos$max, y1=CN_min_pos$CN+0.05, col="#2f4f4f", pch="|", lwd=2, lend=1)
  segments(x0=CN_tot_pos$min, y0=CN_tot_pos$CN-0.05, x1=CN_tot_pos$max, y1=CN_tot_pos$CN-0.05, col="#E69F00", pch="|", lwd=2, lend=1)
  # Plot the vertical lines that show start/end of a chromosome
  chr.segs = lapply(1:length(chr_names), function(ch) { which(BAFvals_all[[1]]$Chromosome==chr_names[ch]) })
  chrk_tot_len = 0
  for (i in 1:length(chr.segs)) {
    chrk = chr.segs[[i]];
    chrk_hetero = names(BAFvals_all[[1]])[chrk]
    chrk_tot_len_prev = chrk_tot_len
    chrk_tot_len = chrk_tot_len + length(chrk_hetero)
    vpos = chrk_tot_len;
    tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
    text(tpos,6.5,chr_names[i], pos = 1, cex = 2)
    abline(v=vpos,lty=1,col="lightgrey")
  }
  dev.off()
  
  png(filename = paste(output.gw.figures.prefix, "_common_CN_pct.png", sep=""), width = 1000, height = 500, res = 200)
  min_table <- data.frame("Region_n" = names(rev(table(consensus_CN_min))), 
                          "Percentage" = round(as.numeric(rev(table(consensus_CN_min)))/length(consensus_CN_min)*100, digits = 2))
  names(min_table) <- c("No. of Regions", "Percentage")
  tot_table <- data.frame("Region_n" = names(rev(table(consensus_CN_tot))), 
                          "Percentage" = round(as.numeric(rev(table(consensus_CN_tot)))/length(consensus_CN_tot)*100, digits = 2))
  names(tot_table) <- c("No. of Regions", "Percentage")
  grid.arrange(textGrob("Minor Allele CN"),
               textGrob("Total CN"),
               tableGrob(min_table, rows = NULL),
               tableGrob(tot_table, rows = NULL), 
               ncol = 2, heights = c(1/5, 4/5), top="Percentage of genome with common copy number")
  dev.off()
}

####################################################################################################################################
### Part 5: Add scDNA to plots
####################################################################################################################################

if (F) {
  #Set % of cells to be considered subclone
  subclone_threshold <- 0.25
  
  #Load in needed files
  spiked_data = T
  if (spiked_data == T) {spiked_file = "_spiked"} else {spiked_file = ""}
  
  samplename = names(TUMOURNAME)
  baf_segmented_file = paste0(if (spiked_data == T) {BB_spiked.dir} else {BB.dir}, samplename, ".BAFsegmented.txt")
  subclones_file = paste0(if (spiked_data == T) {BB_spiked.dir} else {BB.dir}, samplename, "_subclones.txt")
  output.gw.figures.prefix = paste0(paste0(samplename, collapse = "_"), spiked_file, "_BattenbergProfile")
  chr_names = c(1:22,"X")
  
  subclones_all <- lapply(1:length(samplename), function(r) {
    subcloneres = Battenberg::read_table_generic(subclones_file[r])
    subclones = as.data.frame(subcloneres)
    subclones[,2:ncol(subclones)] = sapply(2:ncol(subclones), function(x) { as.numeric(as.character(subclones[,x])) })
    return(subclones)
  })
  BAFvals_all <- lapply(1:length(samplename), function(r) {
    return(as.data.frame(Battenberg:::read_bafsegmented(baf_segmented_file[r])))
  })
  
  chr_probes <- readRDS(paste0(scDNA.dir, "../chr_probes.rds"))
  chr_probes_positions <- lapply(2:nrow(chr_probes), function(r) {
    return((1:chr_probes$total.probes[r] * 500000) - 250000)
  })
  chr_probes_positions <- do.call(rbind, lapply(1:length(chr_probes_positions), function(c) {
    do.call(rbind, lapply(chr_probes_positions[[c]], function (r) {data.frame("CHR" = c, "POS" = r)}))
  }))
  
  #Run plot for each tumour region
  for (i in 1:length(TUMOURNAME)) {
    ###Use freq table for CN 0:10 and rest
    MPNST_scDNA_sc_cn <- readRDS(paste0("MPNST_scDNA_",names(TUMOURNAME)[i],"_sc_cn_freq.rds"))
    MPNST_scDNA_sc_cn <- cbind(c(0:11), MPNST_scDNA_sc_cn) %>% as_tibble
    colnames(MPNST_scDNA_sc_cn) <- c("CN",1:ncol(MPNST_scDNA_sc_cn))
    num_cells <- sum(MPNST_scDNA_sc_cn[,2])
      
    #Remove row for all CN>10
    MPNST_scDNA_sc_cn <- MPNST_scDNA_sc_cn[c(1:11),]
    
    #Retain rows > 25% for position and Keep 2 largest subclones
    MPNST_scDNA_sc_cn_filtered <- lapply(2:ncol(MPNST_scDNA_sc_cn), function(p) {
      sc_tibble <- tibble(CN = unlist(MPNST_scDNA_sc_cn[,1]), Freq = unlist(MPNST_scDNA_sc_cn[,p])) %>% 
        filter(Freq > num_cells*0.25) %>% arrange(desc(Freq))
      if (nrow(sc_tibble) > 1) {
        return(top_n(sc_tibble, 2, Freq))
      } else {
        return(sc_tibble)
      }
    })
    
    #Convert to subclone fractions and add to genomic positions
    MPNST_scDNA_sc_cn_pos <- lapply(MPNST_scDNA_sc_cn_filtered, function (p) {
      ifelse(nrow(p) == 2, 
             sc_profile <- data.frame(CN_1 = pull(p[1,"CN"]),
                                      Frac_1 = prop.table(p["Freq"])[1,1],
                                      CN_2 = pull(p[2,"CN"]),
                                      Frac_2 = prop.table(p["Freq"])[2,1]),
             ifelse(nrow(p) == 1,
                    sc_profile <- data.frame(CN_1 = pull(p[1,"CN"]),
                                             Frac_1 = 1,
                                             CN_2 = NA,
                                             Frac_2 = NA),
                    sc_profile <- data.frame(CN_1 = NA,
                                             Frac_1 = NA,
                                             CN_2 = NA,
                                             Frac_2 = NA)))
      return(sc_profile)
    })
    MPNST_scDNA_sc_cn_pos <- cbind(chr_probes_positions, do.call(rbind, MPNST_scDNA_sc_cn_pos))
    MPNST_scDNA_sc_cn_pos$startpos = MPNST_scDNA_sc_cn_pos$POS-249999
    MPNST_scDNA_sc_cn_pos$endpos = MPNST_scDNA_sc_cn_pos$POS+250000
    
    ###Transfrom into pos_min/pos_max for plotting
      # Map start and end of each segment into the BAF values. The plot uses the index of this BAF table as x-axis
    # sc_pos_min = array(NA, nrow(MPNST_scDNA_sc_cn_pos))
    # sc_pos_max = array(NA, nrow(MPNST_scDNA_sc_cn_pos))
    # for (j in 1:nrow(MPNST_scDNA_sc_cn_pos)) {
    #   segm_chr = MPNST_scDNA_sc_cn_pos$CHR[j] == BAFvals_all[[i]]$Chromosome & 
    #     MPNST_scDNA_sc_cn_pos$startpos[j] < BAFvals_all[[i]]$Position & 
    #     MPNST_scDNA_sc_cn_pos$endpos[j] >= BAFvals_all[[i]]$Position
    #   sc_pos_min[j] = min(which(segm_chr))
    #   sc_pos_max[j] = max(which(segm_chr))
    # }
    
    sc_pos <- pbapply(MPNST_scDNA_sc_cn_pos, 1, function(r) {
      segm_chr = r["CHR"] == BAFvals_all[[i]]$Chromosome & 
        r["startpos"] < BAFvals_all[[i]]$Position & 
        r["endpos"] >= BAFvals_all[[i]]$Position
      sc_pos_min = min(which(segm_chr))
      sc_pos_max = max(which(segm_chr))
      return(c(sc_pos_min, sc_pos_max))
    })
    
    sc_pos_min <- unname(sc_pos[1,])
    sc_pos_max <- unname(sc_pos[2,])

    #Load in data for to battenberg plots
    pos_min <- readRDS(paste0(output.gw.figures.prefix, "_pos_min.rds"))
    pos_max <- readRDS(paste0(output.gw.figures.prefix, "_pos_max.rds"))
    segment_states_min <- readRDS(paste0(output.gw.figures.prefix, "_segment_states_min.rds"))
    segment_states_tot <- readRDS(paste0(output.gw.figures.prefix, "_segment_states_tot.rds"))
    
    ##Make plot and add scDNA subclone lines onto plot
    png(filename = paste0(samplename[i], spiked_file, "_BattenbergProfile", "_average_with_scDNA.png"), width = 2000, height = 500, res = 200)
    ylim = 7
    chr.segs = lapply(1:length(chr_names), function(ch) { which(BAFvals_all[[i]]$Chromosome==chr_names[ch]) })
    # Plot main frame and title
    par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
    maintitle = paste0("Battenberg profile with scDNA subclones")
    plot(c(1,nrow(BAFvals_all[[i]])), c(0,ylim), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    abline(v=0,lty=1,col="lightgrey")
    # Horizontal lines for y=0 to y=5
    abline(h=c(0:ylim),lty=1,col="lightgrey")
    
    #Plot sc CN profiles
    segments(x0=sc_pos_min, y0=MPNST_scDNA_sc_cn_pos$CN_1, x1=sc_pos_max, y1=MPNST_scDNA_sc_cn_pos$CN_1, 
             col="#C77CFF", pch="|", lwd=6*MPNST_scDNA_sc_cn_pos$Frac_1, lend=1)
    segments(x0=sc_pos_min, y0=MPNST_scDNA_sc_cn_pos$CN_2, x1=sc_pos_max, y1=MPNST_scDNA_sc_cn_pos$CN_2, 
             col="#C77CFF", pch="|", lwd=6*MPNST_scDNA_sc_cn_pos$Frac_2, lend=1, )
    
    # Minor allele in gray, total CN in orange
    segments(x0=pos_min[[i]], y0=segment_states_min[[i]], x1=pos_max[[i]], y1=segment_states_min[[i]], col="#2f4f4f", pch="|", lwd=6, lend=1)
    segments(x0=pos_min[[i]], y0=segment_states_tot[[i]], x1=pos_max[[i]], y1=segment_states_tot[[i]], col="#E69F00", pch="|", lwd=6, lend=1)

    # Plot the vertical lines that show start/end of a chromosome
    chrk_tot_len = 0
    for (j in 1:length(chr.segs)) {
      chrk = chr.segs[[j]];
      chrk_hetero = names(BAFvals_all[[i]])[chrk]
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk_hetero)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,ylim,chr_names[j], pos = 1, cex = 2)
      abline(v=vpos,lty=1,col="lightgrey")
    }    
    dev.off()
  }
}


