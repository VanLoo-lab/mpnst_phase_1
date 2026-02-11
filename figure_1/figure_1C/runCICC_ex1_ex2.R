### title: analysis and plotting of mutation clustering
# the code generates Figure 1C in the manuscript (as well as other outputs)

# the code is based on the CICC pipeline for one sample for PCAWG-11
#   written by maxime.tarabichi@crick.ac.uk, 2017-2018
#   more information on: https://github.com/galder-max/CICC

### To run the code:
# (1) Make sure to load SessionInfo() for correct package version
# (2) Set both INPUTDIR and OUTPUTDIR to correct path
# (3) Code necessary to recreate Figure 1c are Parts 0, 1, 2, 3, 4, 7 and 9

# Computational cost is not too high such that the code should be feasible to 
#   run on local machines.

####################################################################################################################################
### Part 0: Code preparation
####################################################################################################################################

# In case plots look weird, the following setting may fix it
#options(bitmapType='cairo')

### Libraries
library(tidyverse)
library(BiocGenerics)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(patchwork)
library(GGally)
library(gridExtra)
library(networkD3)
library(UpSetR)

library(umap)
library(ggrepel)


# number of outlier methods to remove (set to 1)
nbOutlierstoRemove <- 1

# set sample names
samples = c("R1", "R2", "R3", "R4", "R5", "P")
names(samples) = c("VER236A1", "VER236A2", "VER236A3", "VER236A4", "VER236A5", "VER236A6")

# Set code, data input and output directory
CODEDIR <- "~/Documents/GitHub/MPNST-Phase-1/figure_1/figure_1C/"
INPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_1/data/"
OUTPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_1/results/figure_1c/"
setwd(OUTPUTDIR)

# Some additional folder naming for functional code
#cluster.dir <- "/camp/project/proj-vanloo/analyses/hyan/mpnst/bulk/results/DPClust/"
cluster.dir <- paste(INPUTDIR, "bulk/DPClust/", sep="")
#DIRNAME <- "/camp/project/proj-vanloo/analyses/hyan/mpnst/bulk/results/DPClust/CICC_ex1_ex2/"
DIRNAME <- OUTPUTDIR

## number of outlier methods to remove (set to 1)

### Load code from the CICC pipeline
# source CICC functions
source(paste(CODEDIR,"CICC/CICC.PAC.R", sep=""))
# source functions for input and output processing
source(paste(CODEDIR,"CICC/loadData.R", sep=""))
# source dynamic C library for creation of hard-assignment matrix
# the .so file may need to be re-genereated from the .c file depending on the operating system
#   e.g. running "R CMD SHLIB scoringlite.c" in bash
dyn.load(paste(CODEDIR,"CICC/scoringlite.so", sep=""))

####################################################################################################################################
### Part 1: Prepare data loading
####################################################################################################################################

MPNST_load <- function (s, methods) {
  allass <- lapply(methods,function(x) {
    print(x)
    # ex_1 <- unlist(strsplit(methods[1], "_"))[1]
    # ex_2 <- unlist(strsplit(methods[1], "_"))[2]
    t <- read.table(paste0(cluster.dir, "multi_sample_all_ex_", x, "_spiked/MPNST_1_DPoutput_12000iters_2000burnin_seed123/MPNST_1_12000iters_2000burnin_bestConsensusAssignments.bed"), header = T)
    t <-  t[1:(nrow(t)/4),] #DPclust seems to repeat output n times so need to get rid of repeats
    ass <- t[,4]
    names(ass) <- paste(t[,1],t[,3],sep=":")
    return(ass)
  })
  wNull <- which(!sapply(allass,is.null))
  allass <- lapply(wNull,function(x) allass[[x]])
  names(allass) <- methods[wNull]
  nms <- unique(unlist(lapply(allass,names)))
  aa <- lapply(allass,function(x)
  {
    x <- x[nms]
    names(x) <- nms
    x
  })
  passNAcounts <- rowSums(!sapply(aa,is.na))>=round(length(methods)/2+1)
  aa <- lapply(aa,function(x) x[passNAcounts])
  return(aa)  
}

IDS = "MPNST_1"
methods = c(paste0(samples[1], "_", samples[2]), paste0(samples[1], "_", samples[3]), paste0(samples[1], "_", samples[4]),
            paste0(samples[1], "_", samples[5]), paste0(samples[1], "_", samples[6]), paste0(samples[2], "_", samples[3]),
            paste0(samples[2], "_", samples[4]), paste0(samples[2], "_", samples[5]), paste0(samples[2], "_", samples[6]),
            paste0(samples[3], "_", samples[4]), paste0(samples[3], "_", samples[5]), paste0(samples[3], "_", samples[6]),
            paste0(samples[4], "_", samples[5]), paste0(samples[4], "_", samples[6]), paste0(samples[5], "_", samples[6]))

# methods = c(as.vector(samples[1]), as.vector(samples[2]), as.vector(samples[3]), as.vector(samples[4]), as.vector(samples[5]), as.vector(samples[6]), 
#             paste0(samples[1], "_", samples[2]), paste0(samples[1], "_", samples[3]), paste0(samples[1], "_", samples[4]),
#             paste0(samples[1], "_", samples[5]), paste0(samples[1], "_", samples[6]), paste0(samples[2], "_", samples[3]),
#             paste0(samples[2], "_", samples[4]), paste0(samples[2], "_", samples[5]), paste0(samples[2], "_", samples[6]),
#             paste0(samples[3], "_", samples[4]), paste0(samples[3], "_", samples[5]), paste0(samples[3], "_", samples[6]),
#             paste0(samples[4], "_", samples[5]), paste0(samples[4], "_", samples[6]), paste0(samples[5], "_", samples[6]))

####################################################################################################################################
##% Part 2: Loading and saving data for all methods to run on
####################################################################################################################################

print("loading data")
lAA <- lapply(IDS,function(x) try(MPNST_load(x,methods=methods),silent=T))
for (i in 1:length(lAA[[1]])) {
  zero_ccf_clust <- max(lAA[[1]][[i]], na.rm = T)
  lAA[[1]][[i]][lAA[[1]][[i]] == zero_ccf_clust] <- NA
}

save(lAA,file=paste0(DIRNAME,"/", "lAA.nontransformed.",nbOutlierstoRemove,".Rda"))

# remove outliers (here: 1 outlier max and 0 outlier min based on how many an the fraction of mutations they report on
print("removing outliers")
lAA2 <- transformlAA(lAA, nbOutliers=nbOutlierstoRemove)
names(lAA2) <- IDS
save(lAA2,file=paste0(DIRNAME,"/","lAA.transformed.",nbOutlierstoRemove,".Rda"))


####################################################################################################################################
### Part 3: run with different optimal clusters K
####################################################################################################################################

dp_multi_input <- lapply(1:length(samples), function(s) {
  read_tsv(paste0(cluster.dir,"multi_sample_all_ex__spiked/", names(samples[s]),"_allDirichletProcessInfo.txt"), 
           col_types = cols(chr = col_character(),
                            start = col_integer(),
                            end = col_integer(),
                            WT.count = col_integer(),
                            mut.count = col_integer(),
                            subclonal.CN = col_double(),
                            nMaj1 = col_integer(),
                            nMin1 = col_integer(),
                            frac1 = col_double(),
                            nMaj2 = col_integer(),
                            nMin2 = col_integer(),
                            frac2 = col_double(),
                            phase = col_character(),
                            mutation.copy.number = col_double(),
                            subclonal.fraction = col_double(),
                            no.chrs.bearing.mut = col_double()
           ))
})

fastConsensusClustering_MT <- function(allA, nbCluster) {
  nMut <- length(allA[[1]])
  if(nbCluster==1) return(rep(1,nMut))
  nMethods <- length(allA)
  matA <- t(sapply(allA,function(x) x))
  clsA <- NULL
  for(i in 1:nrow(matA))
    clsA <- paste(clsA,matA[i,],sep=if(i==1) "" else "-")
  vuA <- sort(table(clsA),decreasing=T)
  dist1 <- votedist.matrix(names(vuA),nMethods)*max(vuA)*10
  dist2 <- votedist.matrix.nbmut(names(vuA),vuA)
  distF <- as.dist(dist1+dist2)
  hc <- hclust(distF,method="ward.D")
  lClusts <- cutree(hc,k=nbCluster)[clsA]
}

for (k in c(9:20, 25,30,40,50,75,100)) {
  # k = 10
  allRes_k <- fastConsensusClustering_MT(lAA2[["MPNST_1"]], k)
  consensus_clusters_k <- data.frame(chr = as.character(gsub(":.*", "", names(lAA2[["MPNST_1"]][["R1_R2"]]))),
                                     pos = as.integer(gsub(".*:", "", names(lAA2[["MPNST_1"]][["R1_R2"]]))),
                                     cluster = allRes_k)

  saveRDS(consensus_clusters_k, paste0(IDS, "_consensus_clusters_", k, ".rds"))
  # consensus_clusters_k <- readRDS(paste0(IDS, "_consensus_clusters_", k, ".rds"))
  consensus_k <- lapply(1:length(samples), function(s) {
    clusters <- consensus_clusters_k %>% left_join(dp_multi_input[[s]], by= c("chr" = "chr", "pos" = "end")) %>% dplyr::select(chr, pos, cluster, subclonal.fraction)
    return(clusters %>% group_by(cluster) %>% summarise(CCF = round(mean(subclonal.fraction),2), SNVs = n()))
  })
  consensus_CCF_k <- do.call(cbind, consensus_k)[,c(1,2,5,8,11,14,17,18)]
  colnames(consensus_CCF_k) <- c("cluster", samples, "Num SNVs")
  

  png(filename = paste0(IDS, "_consensus_clustering_", k ,".png"), width = 2000, height = 5000, res = 200)
  plot(tableGrob(consensus_CCF_k, rows = NULL))
  dev.off()
  
  #Sankey plot
  if (F) {
    #prev_SNVs <- readRDS("../../DPClust_mutect4.1.2/CICC_ex1_ex2/MPNST_1_consensus_clusters_40.rds")
    prev_SNVs <- readRDS(paste0(INPUTDIR,"CICC_ex1_ex2/MPNST_1_consensus_clusters_40.rds"))
    #Compare to old clusters
    consensus_clusters_k %>% inner_join(prev_SNVs %>% filter(cluster == 7), by = c("chr", "pos")) #R2 subclone cluster
    consensus_clusters_k %>% inner_join(prev_SNVs %>% filter(cluster == 8), by = c("chr", "pos")) #R5 subclone cluster
    
    
    #Incidence matrix
    joint_snv <- consensus_clusters_k %>% rename(rerun = cluster) %>%
      left_join(prev_SNVs, by = c("chr", "pos")) %>% mutate(name = paste0(chr, "_", pos))
    joint_snv$cluster[is.na(joint_snv$cluster)] <- 41
    incidence_mtx <- table(joint_snv[3:4]) %>% as.data.frame.matrix()
    rownames(incidence_mtx) = paste0("rerun_", rownames(incidence_mtx))
    colnames(incidence_mtx) = c("000400", "XXXXXX", "000X00", "0X0000", "X00000", "111111", "040000", "000030", "0000X0", "00X000",
                                "400000", "004000", "X000X0", "00000X", "0XXX00", "222222", "0X0X00", paste0("orig_", 18:40), "New")
    # Transform it to connection data frame with tidyr from the tidyverse:
    links <- incidence_mtx %>% as.data.frame() %>% rownames_to_column(var="source") %>% gather(key="target", value="value", -1) %>% filter(value != 0)
    # From these flows we need to create a node data frame: it lists every entities involved in the flow
    nodes <- data.frame(name=c(as.character(links$source), as.character(links$target)) %>% unique())
    # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
    links$IDsource <- match(links$source, nodes$name)-1 
    links$IDtarget <- match(links$target, nodes$name)-1
    # Make the Network
    sankeyNetwork(Links = links, Nodes = nodes,
                       Source = "IDsource", Target = "IDtarget",
                       Value = "value", NodeID = "name", fontSize = 15,
                       sinksRight=FALSE)

  }
}

####################################################################################################################################
### Part 4: Analysis and plots
####################################################################################################################################

if (T) {
  IDS = "MPNST_1"
  k = 40 #For mutect 4.1.2
  k = 18 #For mutect 4.1.8
  consensus_clusters_k <- readRDS(paste0(IDS, "_consensus_clusters_", k, ".rds"))
  #plot counts per cluster
  consensus_k_counts <- lapply(1:k, function(c) {
    lapply(1:length(samples), function(s) {
      consensus_clusters_k %>% filter(cluster == c) %>% left_join(dp_multi_input[[s]], by= c("chr" = "chr", "pos" = "end")) %>% mutate(region = samples[s])
    })
  })
  
  consensus_k_plot <- do.call(rbind, lapply(consensus_k_counts, function(s) {
    do.call(rbind, s)
  }))
  
  png(filename = paste0(IDS, "_consensus_cluster_", k ,"_mut_counts_violin.png"), width = 2000, height = 1000, res = 200)
  plot(consensus_k_plot %>% ggplot(aes(x = region, y = mut.count, fill = region)) + geom_violin(scale = "width") + 
         ggtitle("Mut read depth per region") + xlab("Region of Origin") + ylab("Read Depth") + ylim(0,150) + theme_minimal()) +
    facet_wrap(~cluster)
  dev.off()
  
  # for (i in 1:k) {
  #   png(filename = paste0(IDS, "_consensus_cluster_", i, "_of_", k ,"_mut_counts.png"), width = 2000, height = 1000, res = 200)
  #   plot(consensus_k_plot %>% filter(cluster == i) %>% ggplot(aes(x = region, y = mut.count, fill = region)) + geom_violin() + 
  #     ggtitle("Mut read depth per region") + xlab("Region of Origin") + ylab("Read Depth") + theme_minimal())
  #   dev.off()
  # }
  
  png(filename = paste0(IDS, "_consensus_cluster_", k ,"_tot_counts_violin.png"), width = 2000, height = 1000, res = 200)
  plot(consensus_k_plot %>% ggplot(aes(x = region, y = WT.count+mut.count, fill = region)) + geom_violin() + 
         ggtitle("Total read depth per region") + xlab("Region of Origin") + ylab("Read Depth") + ylim(0,200) + theme_minimal()) +
    facet_wrap(~cluster)
  dev.off()
  
  ##Try to plot CCFs
  consensus_clusters_CCF_k <- lapply(1:length(samples), function(s) {
    consensus_clusters_k %>% left_join(dp_multi_input[[s]], by= c("chr" = "chr", "pos" = "end")) %>% dplyr::select(chr, pos, cluster, subclonal.fraction)
  })
  
  for (i in 1:length(samples)) {
    plot <- consensus_clusters_CCF_k[[i]] %>% ggplot(aes(x = subclonal.fraction)) + geom_histogram(binwidth = 0.02) + 
      ggtitle(paste0(samples[i], " SNV CCF distribution")) + xlab("CCF") + ylab("Density") + xlim(-0.02,1.52) + theme_minimal() +
      facet_wrap(~cluster, scales = "free")
    
    png(filename = paste0(IDS, "_", samples[i], "_consensus_cluster_", k ,"_CCFs_histo.png"), width = 2000, height = 1000, res = 200)
    print(plot)
    dev.off()
  }
  
  ##Look at CCFs across 2 regions
  for (c in 1:k) {
    cluster_to_look = c
    region_combos <- combn(1:6,2)
    
    png(filename = paste0(IDS, "_consensus_cluster_",cluster_to_look, "of", k ,"_CCFs_all.png"), width = 4000, height = 2000, res = 200)
    plots <- lapply(1:ncol(region_combos), function (p){
      left_join(consensus_clusters_CCF_k[[region_combos[1,p]]] %>% filter(cluster == cluster_to_look),
                consensus_clusters_CCF_k[[region_combos[2,p]]] %>% filter(cluster == cluster_to_look), by = c("chr", "pos", "cluster")) %>%
        ggplot(aes(x=subclonal.fraction.x, y=subclonal.fraction.y)) + geom_point() +
        xlab(samples[region_combos[1,p]]) + ylab(samples[region_combos[2,p]]) + xlim(0,1.5) + ylim(0,1.5) + coord_fixed()
    })
    print(wrap_plots(plots, ncol = 5, nrow = 3, widths = 1, heights = 1) + plot_annotation(title = paste0("CCF across regions for cluster ", cluster_to_look)))
    dev.off()
  }
  
  #Look at those with CCF 0 across regions for cluster 6
  if (F) {
    for (i in 1:length(samples)) {
      cluster_to_look = 6
      cluster_to_look = 1 #for mutect rerun
      region_where_not_zero = i
      region_combos <- combn(1:6,2)
      
      png(filename = paste0(IDS, "_consensus_cluster_",cluster_to_look, "of", k ,"_CCFs_zero_in_",samples[region_where_not_zero],".png"), width = 4000, height = 2000, res = 200)
      plots <- lapply(1:ncol(region_combos), function (p){
        not_zero_in_region <- consensus_clusters_CCF_k[[region_where_not_zero]] %>% filter(subclonal.fraction != 0) %>% dplyr::select(chr, pos, cluster)
        not_zero_in_region %>% inner_join(consensus_clusters_CCF_k[[region_combos[1,p]]], by = c("chr", "pos", "cluster")) %>% filter(cluster == cluster_to_look) %>%
          inner_join(consensus_clusters_CCF_k[[region_combos[2,p]]], by = c("chr", "pos", "cluster")) %>%
          ggplot(aes(x=subclonal.fraction.x, y=subclonal.fraction.y)) + geom_point() +
          xlab(samples[region_combos[1,p]]) + ylab(samples[region_combos[2,p]]) + xlim(0,1.5) + ylim(0,1.5) + coord_fixed()
      })
      print(wrap_plots(plots, ncol = 5, nrow = 3, widths = 1, heights = 1) + plot_annotation(title = paste0("CCF across regions for cluster ", cluster_to_look)))
      dev.off()
    }
  }
  
  #Look at those with CCF 0 across regions
  if (F) {
    regions_clusters_CCFs <- consensus_clusters_CCF_k %>% reduce(left_join, by = c("chr", "pos", "cluster"))
    colnames(regions_clusters_CCFs)[4:9] <- samples
    
    region_subclone_clusters_CCFs <- list()
    for (i in c(1:3,5:6)) { #no need for R4
      cluster_to_look = 6
      cluster_to_look = 1 #for mutect rerun
      region_where_not_zero = i
      other_regions = samples[-i]
      region_combos <- combn(1:6,2)
      
      filter_rule <- paste0(other_regions, " == 0")

      filtered_clusters_CCFs <- regions_clusters_CCFs %>% filter(eval(parse(text = filter_rule[1])), 
                                                                eval(parse(text = filter_rule[2])), 
                                                                eval(parse(text = filter_rule[3])), 
                                                                eval(parse(text = filter_rule[4])), 
                                                                eval(parse(text = filter_rule[5])), 
                                                                eval(parse(text = paste0(samples[i], " > 0"))))

      png(filename = paste0(IDS, "_consensus_cluster_",cluster_to_look, "of", k ,"_CCFs_not_zero_in_",samples[region_where_not_zero],".png"), width = 4000, height = 2000, res = 200)
      plots <- lapply(1:ncol(region_combos), function (p) {
        not_zero_in_region <- filtered_clusters_CCFs %>% dplyr::select(chr, pos, cluster)
        not_zero_in_region %>% inner_join(consensus_clusters_CCF_k[[region_combos[1,p]]], by = c("chr", "pos", "cluster")) %>% filter(cluster == cluster_to_look) %>%
          inner_join(consensus_clusters_CCF_k[[region_combos[2,p]]], by = c("chr", "pos", "cluster")) %>%
          ggplot(aes(x=subclonal.fraction.x, y=subclonal.fraction.y)) + geom_point() +
          xlab(samples[region_combos[1,p]]) + ylab(samples[region_combos[2,p]]) + xlim(0,1.5) + ylim(0,1.5) + coord_fixed()
      })
      print(wrap_plots(plots, ncol = 5, nrow = 3, widths = 1, heights = 1) + plot_annotation(title = paste0("CCF across regions for cluster ", cluster_to_look)))
      dev.off()
      region_subclone_clusters_CCFs[[i]] <- filtered_clusters_CCFs %>% filter(cluster == cluster_to_look)
    }
    saveRDS(region_subclone_clusters_CCFs, paste0(IDS, "_region_subclone_",k,"_clusters_CCFs.rds"))

    #Sankey plot with manually retrieved SNVs
    prev_SNVs <- readRDS("../../DPClust_mutect4.1.2/CICC_ex1_ex2/MPNST_1_consensus_clusters_40.rds")
    region_subclone_clusters_mod <- lapply(which(lapply(region_subclone_clusters_CCFs, length) > 0), function(s) {
      region_subclone_clusters_CCFs[[s]]$cluster <- paste0("Subclone_", samples[s])
      return(region_subclone_clusters_CCFs[[s]][,1:3])
    })
    consensus_clusters_k_mod <- consensus_clusters_k %>% left_join(do.call(rbind, region_subclone_clusters_mod), by = c("chr", "pos")) %>%
      mutate(cluster = ifelse(is.na(cluster.y), cluster.x, cluster.y)) %>% select(-cluster.x, -cluster.y)
    saveRDS(consensus_clusters_k_mod, paste0(IDS, "_consensus_clusters_",k,"_mod.rds"))
    
    #Incidence matrix
    joint_snv <- consensus_clusters_k_mod %>% dplyr::rename(rerun = cluster) %>%
      left_join(prev_SNVs, by = c("chr", "pos")) %>% mutate(name = paste0(chr, "_", pos))
    joint_snv$cluster[is.na(joint_snv$cluster)] <- 41
    incidence_mtx <- table(joint_snv[3:4]) %>% as.data.frame.matrix()
    rownames(incidence_mtx) = paste0("rerun_", rownames(incidence_mtx))
    colnames(incidence_mtx) = c("000400", "XXXXXX", "000X00", "0X0000", "X00000", "111111", "040000", "000030", "0000X0", "00X000",
                                "400000", "004000", "X000X0", "00000X", "0XXX00", "222222", "0X0X00", paste0("orig_", 18:40), "New")
    # Transform it to connection data frame with tidyr from the tidyverse:
    links <- incidence_mtx %>% as.data.frame() %>% rownames_to_column(var="source") %>% gather(key="target", value="value", -1) %>% filter(value != 0)
    # From these flows we need to create a node data frame: it lists every entities involved in the flow
    nodes <- data.frame(name=c(as.character(links$source), as.character(links$target)) %>% unique())
    # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
    links$IDsource <- match(links$source, nodes$name)-1 
    links$IDtarget <- match(links$target, nodes$name)-1
    # Make the Network
    sankeyNetwork(Links = links, Nodes = nodes,
                  Source = "IDsource", Target = "IDtarget",
                  Value = "value", NodeID = "name", fontSize = 15,
                  sinksRight=FALSE) %>% saveNetwork(file = paste0(IDS, "_consensus_clustering_", k ,"_mod.html"))

    
    #Output CCF
    k = 40
    consensus_clusters_k_mod <- readRDS(paste0(IDS, "_consensus_clusters_",k,"_mod.rds"))
    consensus_k_mod <- lapply(1:length(samples), function(s) {
      clusters <- consensus_clusters_k_mod %>% left_join(dp_multi_input[[s]], by= c("chr" = "chr", "pos" = "end")) %>% dplyr::select(chr, pos, cluster, subclonal.fraction)
      return(clusters %>% group_by(cluster) %>% summarise(CCF = round(mean(subclonal.fraction),2), SNVs = n()))
    })
    consensus_CCF_k_mod <- do.call(cbind, consensus_k_mod)[,c(1,2,5,8,11,14,17,18)]
    colnames(consensus_CCF_k_mod) <- c("cluster", samples, "Num SNVs")
    
    
    png(filename = paste0(IDS, "_consensus_clustering_", k ,"_mod.png"), width = 2000, height = 5000, res = 200)
    plot(tableGrob(consensus_CCF_k_mod %>% arrange(desc(`Num SNVs`)), rows = NULL))
    dev.off()
    
    #Tidy up clusters
    good_clusters <- c(3, 10, "Subclone_P", 9, 6, "Subclone_R1", 7, "Subclone_R5", 11, 13, 5, "Subclone_R2", 4, 2, 8, "Subclone_R3")
    renamed_clusters <- c("XXXXXX", "00000X", "000003", "X000X0", "X00000", "300000", "0000X0",  "000030",
                          "0XXX00", "0X0X00", "0X0000", "030000", "000X00", "000400", "00X000", "003000")
    consensus_clusters_k_mod_clean <- consensus_clusters_k_mod %>% 
      filter(cluster %in% good_clusters) %>% 
      left_join(data.frame(cluster = good_clusters,
                           Cluster = renamed_clusters)) %>%
      select(-cluster) %>% rename(cluster = Cluster) 
    saveRDS(consensus_clusters_k_mod_clean, paste0(IDS, "_consensus_clusters_",k,"_mod_clean.rds"))

    consensus_CCF_k_mod_clean <- data.frame(Cluster = good_clusters) %>% 
      left_join(consensus_CCF_k_mod, by = c("Cluster" = "cluster")) %>% rename("No. SNVs" = "Num SNVs")
    consensus_CCF_k_mod_clean$Cluster = renamed_clusters
    pdf(file = paste0(IDS, "_consensus_clustering_", k ,"_mod.pdf"), width = 5, height = 5)
    plot(tableGrob(consensus_CCF_k_mod_clean, rows = NULL))
    dev.off()
  }
  # for (c in 1:ncol(region_combos)) {
  #   png(filename = paste0(IDS, "_", samples[i], "_consensus_cluster_",cluster_to_look, "of", k ,"_CCFs_",samples[region_combos[1,c]],"v",samples[region_combos[2,c]],".png"), width = 1000, height = 1000, res = 200)
  #   left_join(consensus_clusters_CCF_k[[region_combos[1,c]]] %>% filter(cluster == cluster_to_look),
  #             consensus_clusters_CCF_k[[region_combos[2,c]]] %>% filter(cluster == cluster_to_look), by = c("chr", "pos", "cluster")) %>%
  #     ggplot(aes(x=subclonal.fraction.x, y=subclonal.fraction.y)) + geom_point() + ggtitle(paste0("CCF across regions for cluster ", cluster_to_look)) +
  #     xlab(samples[region_combos[1,c]]) + ylab(samples[region_combos[2,c]]) + xlim(0,1.5) + ylim(0,1.5)
  #   dev.off()
  # }
  
  #Check CCF across all regions (parallel coordinates)
  consensus_clusters_SNV_CCF_k <- do.call(cbind, consensus_clusters_CCF_k)[c(1:4,8,12,16,20,24)]
  colnames(consensus_clusters_SNV_CCF_k) <- c("chr", "pos", "cluster", samples)
  consensus_clusters_SNV_CCF_k$cluster <- as.factor(consensus_clusters_SNV_CCF_k$cluster)
  
  png(filename = paste0(IDS, "_consensus_cluster_", k ,"_para_coord_CCFs_all.png"), width = 4000, height = 2000, res = 200)
  ggparcoord(consensus_clusters_SNV_CCF_k, columns = 4:9, groupColumn = 3, scale = "globalminmax",
             title = "CCF of all clusters of SNVs across regions") + ylim(0,1.5) + theme_minimal()
  dev.off()
  
  for (c in 1:k) {
    png(filename = paste0(IDS, "_consensus_cluster_",c, "of", k ,"_para_coord_CCFs.png"), width = 4000, height = 2000, res = 200)
    print(ggparcoord(consensus_clusters_SNV_CCF_k %>% filter(cluster == c), columns = 4:9, groupColumn = 3, scale = "globalminmax",
                     title = paste0("CCF of cluster ",c , " across regions")) + 
            geom_point(data = pivot_longer(consensus_CCF_k[c,], cols = 2:7, names_to = "region", values_to = "CCF"), aes(x = region, y = CCF), inherit.aes = F) + 
            ylim(0,1.5) + theme_minimal())
    dev.off()
  }
  
  #Look at CCFs of clonal cluster
  clusters_to_plot <- c(2)
  consensus_clusters_CCF_k_plot <- do.call(rbind, lapply(1:length(samples), function(r) {
    consensus_clusters_CCF_k[[r]] %>% mutate(region = samples[r]) %>% 
      arrange(factor(cluster, levels = clusters_to_plot)) %>% rowid_to_column(var="ID")
  })) %>% filter(cluster %in% clusters_to_plot)
  
  consensus_clusters_CCF_k_wide <- consensus_clusters_CCF_k_plot %>% mutate(subclonal.fraction = ifelse(subclonal.fraction>0,1,0)) %>% pivot_wider(names_from = region, values_from = subclonal.fraction)
  consensus_clusters_CCF_k_upset <- consensus_clusters_CCF_k_wide %>% filter(R1 == 0 | R2 == 0 | R3 == 0 | R4 == 0 | R5 == 0 | P == 0)
  consensus_clusters_CCF_k_upset_list <- lapply(samples, function(s) {
    return(consensus_clusters_CCF_k_upset %>% dplyr::select(ID,chr,pos,s) %>% filter(consensus_clusters_CCF_k_upset[,s] > 0) %>% mutate(SNV_name = paste0(chr, "_", pos)) %>% pull(SNV_name))
  })
  names(consensus_clusters_CCF_k_upset_list) <- samples
  
  png(filename = paste0("Truncal_SNVs_upset.png"), width = 2000, height = 2000, res = 200)
  upset(fromList(consensus_clusters_CCF_k_upset_list), 
        nintersects = 40, 
        nsets = length(consensus_clusters_CCF_k_upset_list), 
        sets = rev(sort(names(consensus_clusters_CCF_k_upset_list))), #sets primary first
        order.by = "freq", 
        decreasing = T, 
        mb.ratio = c(0.7, 0.3),
        number.angles = 0,
        text.scale = 1.5,
        point.size = 2.2, 
        line.size = 0.7,
        keep.order = TRUE)
  dev.off()
  
  #Extract overlap groups
  overlapGroups <- function (listInput, sort = TRUE) {
    # listInput could look like this:
    # $one
    # [1] "a" "b" "c" "e" "g" "h" "k" "l" "m"
    # $two
    # [1] "a" "b" "d" "e" "j"
    # $three
    # [1] "a" "e" "f" "g" "h" "i" "j" "l" "m"
    listInputmat    <- fromList(listInput) == 1
    #     one   two three
    # a  TRUE  TRUE  TRUE
    # b  TRUE  TRUE FALSE
    #...
    # condensing matrix to unique combinations elements
    listInputunique <- unique(listInputmat)
    grouplist <- list()
    print(paste0("Going through ",nrow(listInputunique)," combinations"))
    # going through all unique combinations and collect elements for each in a list
    for (i in 1:nrow(listInputunique)) {
      if (i%%10 == 0) {print(i)}
      currentRow <- listInputunique[i,]
      myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
      attr(myelements, "groups") <- currentRow
      grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
      myelements
      # attr(,"groups")
      #   one   two three 
      # FALSE FALSE  TRUE 
      #  f  i 
      # 12 13 
    }
    if (sort) {
      grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
    }
    attr(grouplist, "elements") <- unique(unlist(listInput))
    return(grouplist)
    # save element list to facilitate access using an index in case rownames are not named
  }
  
  consensus_clusters_CCF_overlap_groups <- overlapGroups(consensus_clusters_CCF_k_upset_list)
  saveRDS(consensus_clusters_CCF_overlap_groups, "MPNST_bulk_SNVs_overlap_groups.rds")
  #To get SNV names:
  attr(consensus_clusters_CCF_overlap_groups, "elements")[consensus_clusters_CCF_overlap_groups[["R1:R2:R3:R4:R5"]]]
  saveRDS(attr(consensus_clusters_CCF_overlap_groups, "elements")[consensus_clusters_CCF_overlap_groups[["R1:R2:R3:R4:R5"]]], "MPNST_recurrence_SNVs.rds")
  
  
  Recurrence_CCF <- consensus_clusters_CCF_k_plot %>% pivot_wider(names_from = region, values_from = subclonal.fraction) %>%
    mutate(name = paste0(chr,"_",pos)) %>% filter(name %in% attr(consensus_clusters_CCF_overlap_groups, "elements")[consensus_clusters_CCF_overlap_groups[["R1:R2:R3:R4:R5"]]])
  
  png(filename = paste0("Recurrence_SNVs_CCF.png"), width = 2000, height = 500, res = 200)
  Recurrence_CCF_plot <- do.call(cbind, lapply(unname(samples), function(s) {
    round(mean(Recurrence_CCF %>% pull(s)),2)
  }))
  colnames(Recurrence_CCF_plot) <- samples
  plot(tableGrob(cbind(Recurrence_CCF_plot, data.frame("Num SNVs" = nrow(Recurrence_CCF))), rows = NULL))
  dev.off()
}
####################################################################################################################################
### Part 5: Trinucleotide context of clusters of mutation
####################################################################################################################################

if (F) {
  library(MutationalPatterns)
  library(BSgenome)
  ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
  library(BSgenome.Hsapiens.UCSC.hg38)
  vcf.input.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/bulk/results/snv_mnv_indel_tumouronly/", names(samples), "/")
  vcf_files <- paste0(vcf.input.dir,names(samples),"_PASS_snvs_indels.vcf.gz")
  
  #Vignette code
  if (F) {
    vcfs <- read_vcfs_as_granges(vcf_files, samples, ref_genome)
    # muts = mutations_from_vcf(vcfs[[1]]) #base substitutions
    # types = mut_type(vcfs[[1]]) #convert to 6 types of base substitutions
    # context = mut_context(vcfs[[1]], ref_genome) #Trinucleotide context
    # type_context = type_context(vcfs[[1]], ref_genome) #Get type and context
    
    # type_occurrences <- mut_type_occurrences(vcfs, ref_genome) #Get type and context freq (useful one)
    # plot_spectrum(type_occurrences, by = samples) #6 types of mutations plot
    mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome) #Make 96 matrix
    plot_96_profile(mut_mat, ymax = 0.1, condensed = T) #Make 96 plot
  }
  
  vcfs <- read_vcfs_as_granges(vcf_files, samples, ref_genome)
  combined_vcf <- unlist(vcfs) 
  combined_vcf <- combined_vcf[!duplicated(combined_vcf)]
  
  for (c in 1:k) {
    cluster_to_look = 17
    cluster_snvs <- GRanges(seqnames = paste0("chr", consensus_clusters_k %>% filter(cluster == cluster_to_look) %>% pull(chr)),
                            ranges = IRanges(start = consensus_clusters_k %>% filter(cluster == cluster_to_look) %>% pull(pos),
                                             end = consensus_clusters_k %>% filter(cluster == cluster_to_look) %>% pull(pos)))
    cluster_vcfs <- combined_vcf[combined_vcf %in% cluster_snvs]
    print(paste0(length(cluster_vcfs), " SNVs in cluster"))
    mut_mat <- mut_matrix(vcf_list = list(cluster_vcfs), ref_genome = ref_genome) #Make 96 matrix
    png(filename = paste0(IDS, "_", samples[i], "_consensus_cluster_",cluster_to_look, "of", k ,"_context.png"), width = 4000, height = 2000, res = 200)
    print(plot_96_profile(mut_mat, ymax = 0.1, condensed = T)) #Make 96 plot
    dev.off()
  }  
  
  #Combine cluster 18:40
  cluster_to_look = 18:40
  cluster_snvs <- GRanges(seqnames = paste0("chr", consensus_clusters_k %>% filter(cluster %in% cluster_to_look) %>% pull(chr)),
                          ranges = IRanges(start = consensus_clusters_k %>% filter(cluster %in% cluster_to_look) %>% pull(pos),
                                           end = consensus_clusters_k %>% filter(cluster %in% cluster_to_look) %>% pull(pos)))
  cluster_vcfs <- combined_vcf[combined_vcf %in% cluster_snvs]
  print(paste0(length(cluster_vcfs), " SNVs in cluster"))
  mut_mat <- mut_matrix(vcf_list = list(cluster_vcfs), ref_genome = ref_genome) #Make 96 matrix
  png(filename = paste0(IDS, "_", samples[i], "_consensus_cluster_18_to_40_context.png"), width = 4000, height = 2000, res = 200)
  print(plot_96_profile(mut_mat, ymax = 0.1, condensed = T)) #Make 96 plot
  dev.off()
}

####################################################################################################################################
### Part 6: Quick try to see UMAP projections of clusters
####################################################################################################################################

if (F) {
  MPNST_umap_input <- cbind(do.call(cbind, lAA2[[1]]), cluster = consensus_clusters_k[,"cluster"])
  MPNST_umap_data <- as.data.frame(do.call(cbind, lAA2[[1]]))
  MPNST_umap_data[is.na(MPNST_umap_data)] <- 0#Replace NAs with Os as can't calculate correlation with NAs
  MPNST_umap_labels <- consensus_clusters_k[,"cluster"]
  
  ##UMAP of correlation matrix (feed in data, let it calculate correlation distance)
  umap.settings <- umap.defaults
  umap.settings$input = "data"
  umap.settings$metric = "euclidean"
  umap.settings$verbose = T
  
  MPNST_umap <-umap(MPNST_umap_data, config=umap.settings)
  
  tibble(SNVs = rownames(MPNST_umap[["layout"]]), x = MPNST_umap$layout[,1], y = MPNST_umap$layout[,2], ) %>%
    cbind(cluster = as.factor(MPNST_umap_labels)) -> MPNST_umap_plot
  # MPNST_obs_umap_plot$x <- -(MPNST_obs_umap_plot$x)
  # MPNST_obs_umap_plot$y <- -(MPNST_obs_umap_plot$y)
  
  gg_default_colours <- hcl(h = seq(15, 375, length = k + 1), l = 65, c = 100)[1:k]
  gg_colours <- gg_default_colours
  barplot(1:k, col = gg_colours)
  
  png(filename = paste0("MPNST_CICC_cluster_UMAP.png"), width = 2000, height = 2000, res = 200)
  ggplot(MPNST_umap_plot, aes(x, y)) + geom_point(size = 0.5) + 
    theme_classic()
  dev.off()
  
  png(filename = paste0("MPNST_CICC_cluster_UMAP_lab.png"), width = 2000, height = 2000, res = 200)
  ggplot(MPNST_umap_plot, aes(x, y, colour = cluster)) + geom_point(size = 0.5) + 
    scale_color_manual(values = gg_colours) + theme_classic()
  dev.off()
  
  MPNST_umap_plot_labels <- MPNST_umap_plot %>% group_by(cluster) %>% summarise(x=median(x), y=median(y))
  png(filename = paste0("MPNST_CICC_cluster_UMAP_lab_text.png"), width = 2000, height = 2000, res = 200)
  ggplot(MPNST_umap_plot, aes(x, y, colour = cluster)) + geom_point(size = 0.5) + 
    scale_color_manual(values = gg_colours) + geom_label_repel(data = MPNST_umap_plot_labels, aes(label = cluster)) + theme_classic()
  dev.off()
}

####################################################################################################################################
### Part 7: Prep inputs (needed for plots below)
####################################################################################################################################

if (T) {
  #Load in CN profiles and CCFs
  dp_multi_input <- lapply(1:length(samples), function(s) {
    read_tsv(paste0(cluster.dir, "multi_sample_all_ex__spiked/", names(samples[s]),"_allDirichletProcessInfo.txt"), 
             col_types = cols(chr = col_character(),
                              start = col_integer(),
                              end = col_integer(),
                              WT.count = col_integer(),
                              mut.count = col_integer(),
                              subclonal.CN = col_double(),
                              nMaj1 = col_integer(),
                              nMin1 = col_integer(),
                              frac1 = col_double(),
                              nMaj2 = col_integer(),
                              nMin2 = col_integer(),
                              frac2 = col_double(),
                              phase = col_character(),
                              mutation.copy.number = col_double(),
                              subclonal.fraction = col_double(),
                              no.chrs.bearing.mut = col_double()
             ))
  })
  IDS = "MPNST_1"
  k = 40
  consensus_clusters_k <- readRDS(paste0(cluster.dir, "CICC_ex1_ex2/" ,IDS, "_consensus_clusters_", k, ".rds"))
  consensus_clusters_CCF_k <- lapply(1:length(samples), function(s) {
    consensus_clusters_k %>% left_join(dp_multi_input[[s]], by= c("chr" = "chr", "pos" = "end")) %>% dplyr::select(chr, pos, cluster, subclonal.fraction)
  })
  consensus_clusters_k_mod <- readRDS(paste0(cluster.dir, "CICC_ex1_ex2/", IDS, "_consensus_clusters_",k,"_mod.rds"))
  
  cluster_names <- c("000400", "XXXXXX", "000X00", "0X0000", "X00000", "111111", "040000", "000030", "0000X0", "00X000",
                     "400000", "004000", "X000X0", "00000X", "0XXX00", "222222", "0X0X00")
  cluster_names <- c("111111", "000400", "XXXXXX", "000X00", "0X0000", "X00000", "0000X0", "00X000", "X000X0", "00000X",
                     "0XXX00", "222222", "0X0X00") #Rerun
  
}

####################################################################################################################################
### Part 8: Plot SNVs along genome
####################################################################################################################################

if (F) {
  plotCNProfile <- function(cn,ccfs,region,cluster) {
    col1 <- rgb(.7,.3,.2,.5)
    col2 <- rgb(.5,.5,.5,.5)
    lwd <- 2
    chrs <- c(1:22,"X")
    cs <- cumsum(sapply(chrs,function(x) {
      maxend <- max(cn[as.character(cn$chr)==x,"endpos"],na.rm=T)
      as.numeric(maxend)/10000
    }))
    maxX <- max(cs)
    maxY <- max(cn$nMaj1_A+cn$nMin1_A,na.rm=T)+1.5
    if(maxY>6) maxY <- 6
    xlims <- c(0,maxX)
    ylims <- c(0,maxY)
    par(mar=c(4,4,3,5))
    plot(0,0,xlab="", ylab="Copy Number States", main=paste0(region,": Cluster ", cluster_names[cluster], " CCFs"),
         col=rgb(0,0,0,0), xlim=xlims, ylim=ylims)
    ats <- seq(0,maxY+0.01,maxY/5)
    axis(side=4,ylab="CCFs",at=ats,as.character(signif(ats/maxY*2,2)))
    ccfs <- ccfs[[cluster]]
    addCCFs(ccfs,c(0,cs), maxX, maxY)
    allends <- unlist(lapply(1:nrow(cn),function(x) {
      wc <- which(chrs==cn$chr[x])
      wc <- wc[!is.na(wc)]
      wc <- wc[1]
      cn$endpos[x]/10000+c(0,cs)[wc]
    }))
    allstarts <- allends-(cn$endpos-cn$startpos)/10000
    # Determine whether it's the major or the minor allele that is represented by two states
    is_subclonal_maj = abs(cn$nMaj1_A - cn$nMaj2_A) > 0
    is_subclonal_min = abs(cn$nMin1_A - cn$nMin2_A) > 0
    is_subclonal_maj[is.na(is_subclonal_maj)] = F
    is_subclonal_min[is.na(is_subclonal_min)] = F
    #Calculate subclonal states
    segment_states_min = cn$nMin1_A * ifelse(is_subclonal_min, cn$frac1_A, 1)  + ifelse(is_subclonal_min, cn$nMin2_A, 0) * ifelse(is_subclonal_min, cn$frac2_A, 0) 
    segment_states_maj = cn$nMaj1_A * ifelse(is_subclonal_maj, cn$frac1_A, 1)  + ifelse(is_subclonal_maj, cn$nMaj2_A, 0) * ifelse(is_subclonal_maj, cn$frac2_A, 0) 
    segment_states_tot = segment_states_maj + segment_states_min
    #Plot segments
    segments(allstarts, segment_states_tot, allends, segment_states_tot,col=col1,lwd=lwd)
    segments(allstarts, segment_states_min, allends, segment_states_min,col=col2,lwd=lwd)
    abline(v=cs,lty=2,col=rgb(.5,.5,.5,.5))
    return(list(cs=c(0,cs),maxX=maxX,maxY=maxY))
  }
  
  addCCFs <- function(ccfs,cs,maxX,maxY) {
    chrs <- c(1:22,"X")
    xx <- sapply(1:nrow(ccfs),function(x)
      ccfs$pos[x]/10000+cs[which(chrs==ccfs$chr[x])])
    points(xx,ccfs$subclonal.fraction*maxY/2,col=rgb(.5,.5,.5,.5),pch=19,cex=.4)
  }
  
  for (r in 1:length(samples)) {
    consensus_CCF_by_clusters <- lapply(1:k, function(c) {
      cluster <- consensus_clusters_k %>% filter(cluster == c) %>% 
        left_join(consensus_clusters_CCF_k[[r]], by = c("chr", "pos")) %>% dplyr::select(chr, pos, subclonal.fraction)
    })
    names(consensus_CCF_by_clusters) <- cluster_names
    CN_profile <- read.delim(paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/bulk/results/BB_spiked/All_with_Ext/", samples[r], "_subclones.txt"))
    for (c in 1:length(cluster_names)) {
      png(filename = paste0("MPNST_",samples[r],"_cluster_",cluster_names[c],"_CCFvCN.png"), width = 4000, height = 2000, res = 200)
      plotCNProfile(CN_profile, consensus_CCF_by_clusters, samples[r], c)
      dev.off()
    }
  }
}

####################################################################################################################################
### Part 9: Plot CCF for samples
####################################################################################################################################

if (T) {
  if (T) {
    clusters_to_plot <- c(2,6,16,14, 13,5,11,9,8, 15,17,4,7,3,1,10,12)
    cluster_colours <- c("deepskyblue", "black", "#619CFF", "#00BA38", "#B79F00", "grey30", "chartreuse", "pink", "#F564E3", "#00BFC4",
                         "gold", "cyan", "coral", "#F8766D", "darkslategray4", "grey70", "darkgreen")
    consensus_clusters_CCF_k_plot <- do.call(rbind, lapply(1:length(samples), function(r) {
      consensus_clusters_CCF_k[[r]] %>% mutate(region = samples[r]) %>% 
        arrange(factor(cluster, levels = clusters_to_plot)) %>% rowid_to_column(var="ID")
    })) %>% filter(cluster %in% clusters_to_plot)
  } else {
    clusters_to_plot <- c(3, 1, 12, 10, "Subclone_P", 9, 6, "Subclone_R1", 7, "Subclone_R5", 11, 13, 5, "Subclone_R2", 4, 2, 8, "Subclone_R3")
    cluster_colours <- c("black", "grey70", "grey30", "#F8766D", "pink", "coral", "#B79F00", "gold", "#F564E3", "mediumpurple1", 
                         "darkslategray4", "darkgreen", "#00BA38", "chartreuse", "#619CFF", "deepskyblue", "#00BFC4", "cyan")
    names(cluster_colours) <- clusters_to_plot
    consensus_clusters_CCF_k_plot <- do.call(rbind, lapply(1:length(samples), function(r) {
      consensus_clusters_CCF_k[[r]] %>% mutate(region = samples[r], cluster = consensus_clusters_k_mod$cluster) %>% 
        arrange(factor(cluster, levels = clusters_to_plot)) %>% rowid_to_column(var="ID")
    })) %>% filter(cluster %in% clusters_to_plot)
  }
  
  png(filename = paste0("MPNST_all_region_all_clusters_CCF.png"), width = 12000, height = 6000, res = 200)
  consensus_clusters_CCF_k_plot %>% arrange(region) %>% 
    ggplot(aes(x = ID, y = subclonal.fraction, fill = as.factor(cluster))) + geom_bar(stat = "identity", width = 1) + scale_fill_manual(values = cluster_colours) +
    scale_y_continuous(breaks = c(0,0.5,1)) + coord_cartesian(ylim=c(0, 1)) + xlab("Variant") + ylab("CCF") + facet_wrap(~region, ncol = 1, strip.position = "right") +
    theme(text=element_text(size=60), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.y.left = element_line(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background = element_blank(), panel.spacing = unit(2, "lines"))
  dev.off() 
  
  if (T) {
    clusters_to_plot <- c(2,14, 13,5,11,9,8, 15,17,4,7,3,1,10,12) #without artefact clusters
    cluster_colours <- c("deepskyblue", "black", "#619CFF", "#00BA38", "#B79F00", "chartreuse", "pink", "#F564E3", "#00BFC4",
                         "gold", "cyan", "coral", "#F8766D", "darkslategray4", "darkgreen")
    consensus_clusters_CCF_k_plot <- do.call(rbind, lapply(1:length(samples), function(r) {
      consensus_clusters_CCF_k[[r]] %>% mutate(region = samples[r]) %>% 
        arrange(factor(cluster, levels = clusters_to_plot)) %>% rowid_to_column(var="ID")
    })) %>% filter(cluster %in% clusters_to_plot)
  } else {
    clusters_to_plot <- c(3, 10, "Subclone_P", 9, 6, "Subclone_R1", 7, "Subclone_R5", 11, 13, 5, "Subclone_R2", 4, 2, 8, "Subclone_R3")
    cluster_colours <- c("black", "#F8766D", "pink", "coral", "#B79F00", "gold", "#F564E3", "mediumpurple1", 
                         "darkslategray4", "darkgreen", "#00BA38", "chartreuse", "#619CFF", "deepskyblue", "#00BFC4", "cyan")
    names(cluster_colours) <- clusters_to_plot
    consensus_clusters_CCF_k_plot <- do.call(rbind, lapply(1:length(samples), function(r) {
      consensus_clusters_CCF_k[[r]] %>% mutate(region = samples[r], cluster = consensus_clusters_k_mod$cluster) %>% 
        arrange(factor(cluster, levels = clusters_to_plot)) %>% rowid_to_column(var="ID")
    })) %>% filter(cluster %in% clusters_to_plot)
  }
  
  png(filename = paste0("MPNST_all_region_good_clusters_CCF.png"), width = 12000, height = 6000, res = 200)
  consensus_clusters_CCF_k_plot %>% arrange(region) %>% 
    ggplot(aes(x = ID, y = subclonal.fraction, fill = as.factor(cluster))) + geom_bar(stat = "identity", width = 1) + scale_fill_manual(values = cluster_colours) +
    scale_y_continuous(breaks = c(0,0.5,1)) + coord_cartesian(ylim=c(0, 1)) + xlab("Variant") + ylab("CCF") + facet_wrap(~region, ncol = 1, strip.position = "right") +
    theme(text=element_text(size=60), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.y.left = element_line(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background = element_blank(), panel.spacing = unit(2, "lines"))
  dev.off() 
  
  pdf(file = paste0("MPNST_all_region_good_clusters_CCF.pdf"), width = 14, height = 7)
  consensus_clusters_CCF_k_plot %>% arrange(region) %>% 
    ggplot(aes(x = ID, y = subclonal.fraction, fill = as.factor(cluster))) + geom_bar(stat = "identity", width = 1) + scale_fill_manual(values = cluster_colours) +
    scale_y_continuous(breaks = c(0,1)) + coord_cartesian(ylim=c(0, 1)) + xlab("Variant") + ylab("CCF") + facet_wrap(~region, ncol = 1, strip.position = "right") +
    theme(text=element_text(size=28), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.y.left = element_line(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background = element_blank(), panel.spacing = unit(2, "lines"))
  dev.off() 
}

####################################################################################################################################
### Part 10: Reassign cluster based on next closest cluster?
####################################################################################################################################
if (F) {
  cluster_to_reassign <- 11
  
  nMut <- length(lAA2[["MPNST_1"]][[1]])
  # if(nbCluster==1) return(rep(1,nMut))
  nMethods <- length(lAA2[["MPNST_1"]])
  matA <- t(sapply(lAA2[["MPNST_1"]],function(x) x))
  clsA <- NULL
  for(i in 1:nrow(matA))
    clsA <- paste(clsA,matA[i,],sep=if(i==1) "" else "-") #generate combinations of votes
  vuA <- sort(table(clsA),decreasing=T) #frequency table of votes
  dist1 <- votedist.matrix(names(vuA),nMethods)*max(vuA)*10 #matrix of combos. each value is sum of votes?
  dist2 <- votedist.matrix.nbmut(names(vuA),vuA)
  distT <- dist1+dist2
  distF <- as.dist(dist1+dist2)
  hc <- hclust(distF,method="ward.D")
  lClusts <- cutree(hc,k=11)[clsA]
  
  clsA_rm <- names(lClusts[lClusts == 11])
  clsA_rm_Clusts <- pbmclapply(1:length(clsA_rm), function(x) {
    vuA_sub <- sort(table(c(clsA[!clsA %in% clsA_rm], clsA_rm[x])),decreasing=T)
    dist1_sub <- votedist.matrix(names(vuA_sub),nMethods)*max(vuA)*10 #matrix of combos. each value is sum of votes?
    dist2_sub <- votedist.matrix.nbmut(names(vuA_sub),vuA_sub)
    # distT_sub <- dist1+dist2
    distF_sub <- as.dist(dist1_sub+dist2_sub)
    hc_sub <- hclust(distF_sub,method="ward.D")
    lClusts_sub <- cutree(hc_sub,k=k-1)[c(clsA[!clsA %in% clsA_rm], clsA_rm[x])]
    lClusts_sub[clsA_rm[x]]
  }, mc.cores = 20)

  #Extract reassigned cluster numbers
  lClusts_reassigned <- left_join(enframe(lClusts) %>% mutate(ID = names(lAA2[["MPNST_1"]][[1]])), 
                                  enframe(unlist(clsA_rm_Clusts)) %>% mutate(ID = names(lAA2[["MPNST_1"]][[1]])[lClusts == 11]), 
                                  by = c("ID", "name"))
  lClusts_reassigned <- lClusts_reassigned %>% mutate(cluster = ifelse(is.na(value.y), value.x, value.y)) %>% pull(cluster) %>% as.integer()
  
  consensus_clusters_k_rm <- data.frame(chr = as.character(gsub(":.*", "", names(lAA2[["MPNST_1"]][["R1_R2"]]))),
                                     pos = as.integer(gsub(".*:", "", names(lAA2[["MPNST_1"]][["R1_R2"]]))),
                                     cluster = lClusts_reassigned)
  saveRDS(consensus_clusters_k_rm, paste0(IDS, "_consensus_clusters_", k, "_removed_cluster_",cluster_to_reassign, ".rds"))
  
  consensus_k_rm <- lapply(1:length(samples), function(s) {
    clusters <- consensus_clusters_k_rm %>% left_join(dp_multi_input[[s]], by= c("chr" = "chr", "pos" = "end")) %>% select(chr, pos, cluster, subclonal.fraction)
    return(clusters %>% group_by(cluster) %>% summarise(CCF = round(mean(subclonal.fraction),3), SNVs = n()))
  })
  consensus_CCF_k_rm <- do.call(cbind, consensus_k_rm)[,c(1,2,5,8,11,14,17,18)]
  colnames(consensus_CCF_k_rm) <- c("cluster", samples, "Num SNVs")
  
  png(filename = paste0(IDS, "_consensus_clustering_", k , "_removed_cluster_",cluster_to_reassign, ".png"), width = 2000, height = 2000, res = 200)
  plot(tableGrob(consensus_CCF_k_rm, rows = NULL))
  dev.off()
}

fastConsensusClustering_MT <- function(allA, nbCluster) {
  nMut <- length(allA[[1]])
  if(nbCluster==1) return(rep(1,nMut))
  nMethods <- length(allA)
  matA <- t(sapply(allA,function(x) x))
  clsA <- NULL
  for(i in 1:nrow(matA)) {
    clsA <- paste(clsA,matA[i,],sep=if(i==1) "" else "-")
  }
    
  vuA <- sort(table(clsA),decreasing=T)
  dist1 <- votedist.matrix(names(vuA),nMethods)*max(vuA)*10
  dist2 <- votedist.matrix.nbmut(names(vuA),vuA)
  distF <- as.dist(dist1+dist2)
  hc <- hclust(distF,method="ward.D")
  lClusts <- cutree(hc,k=nbCluster)[clsA]
}

