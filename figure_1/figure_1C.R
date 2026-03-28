### title: plotting clustered subclonal architecture based on CCF values

# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"

# Define input directory with data and output directory to save the figure
input.dir <- paste0(zenodo.dir, "data/bulk/CCF_clustering/")
output.dir <- paste0(zenodo.dir, "results/figure_1/")
setwd(output.dir)

####################################################################################################################################
### Part 0: Load libraries and data
####################################################################################################################################

### Load libraries
library(tidyverse)

### Load data
consensus_clusters_k_mod <- readRDS(paste0(input.dir, "consensus_clusters_k_mod.rds"))
consensus_clusters_CCF_k<- readRDS(paste0(input.dir, "consensus_clusters_CCF_k.rds"))

####################################################################################################################################
### Part 1: Make the figure
####################################################################################################################################

### Set sample names
samples = c("R1", "R2", "R3", "R4", "R5", "P")
names(samples) = c("VER236A1", "VER236A2", "VER236A3", "VER236A4", "VER236A5", "VER236A6")

### Define cluster order manually
clusters_to_plot <- c(3, 10, "Subclone_P", 9, 6, "Subclone_R1", 7, "Subclone_R5", 11, 13, 5, "Subclone_R2", 4, 2, 8, "Subclone_R3")

### Define color palette
cluster_colours <- c("black", "#F8766D", "pink", "coral", "#B79F00", "gold", "#F564E3", "mediumpurple1", 
                     "darkslategray4", "darkgreen", "#00BA38", "chartreuse", "#619CFF", "deepskyblue", "#00BFC4", "cyan")
names(cluster_colours) <- clusters_to_plot

### Create a data frame for plotting with ggplot2
consensus_clusters_CCF_k_plot <- do.call(rbind, lapply(1:length(samples), function(r) {
  consensus_clusters_CCF_k[[r]] %>% mutate(region = samples[r], cluster = consensus_clusters_k_mod$cluster) %>% 
    arrange(factor(cluster, levels = clusters_to_plot)) %>% rowid_to_column(var="ID")
})) %>% filter(cluster %in% clusters_to_plot)

### Plotting
pdf(file = paste0("figure_1C.pdf"), width = 14, height = 7)
consensus_clusters_CCF_k_plot %>% arrange(region) %>%
  ggplot(aes(x = ID, y = subclonal.fraction, fill = as.factor(cluster))) + geom_bar(stat = "identity", width = 1) + scale_fill_manual(values = cluster_colours) +
  scale_y_continuous(breaks = c(0,1)) + coord_cartesian(ylim=c(0, 1)) + xlab("Variant") + ylab("CCF") + facet_wrap(~region, ncol = 1, strip.position = "right") +
  theme(text=element_text(size=28), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line.y.left = element_line(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background = element_blank(), panel.spacing = unit(2, "lines"))
dev.off()
