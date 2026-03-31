### title: Cell type deconvolution of spatial transcriptomics data

# Define the zenodo repository containing input and output folders
zenodo.dir <- "~/Documents/GitHub/MPNST-Zenodo/"

# Define input directory with data and output directory to save the figure
input.dir <- paste0(zenodo.dir, "data/spRNA/samples/")
output.dir <- paste0(zenodo.dir, "results/figure_4/")
setwd(output.dir)

####################################################################################################################################
### Part 0: Load libraries and data
####################################################################################################################################

### Load libraries
library(tidyverse)
library(CARD)
library(SingleCellExperiment)

### Define sample names
samples = c("Pa", "R1", "R3", "Pb", "R2", "R4", "R5a", "R5b")
names(samples) = c("VER683A1", "VER683A2", "VER683A3", "VER683A4", "VER683A5", "VER683A6", "VER683A7", "VER683A8")

### Load data immediately in for loop below
# see sc_count, sc_meta, sp_count, sp_location for each sample in the loop below

####################################################################################################################################
### Part 1: Load scDNA and spRNA data, run CARD and plot cell type deconcolution
####################################################################################################################################

# Note that the code is unable to overwrite files and errors if files already exist

for (sample in samples[1:8]) {
  ### Load sc and sp data
  sc_count <- readRDS(paste0(input.dir, sample, "_sc_count.rds"))
  sc_meta <- readRDS(paste0(input.dir, sample, "_sc_meta.rds"))
  sp_count <- readRDS(paste0(input.dir, sample, "_sp_count.rds"))
  sp_location <- readRDS(paste0(input.dir, sample, "_sp_location.rds"))
  
  ### Create CARD object
  CARD_obj = createCARDObject(
    sc_count = sc_count,
    sc_meta = sc_meta,
    spatial_count = sp_count,
    spatial_location = sp_location,
    ct.varname = "cellType",
    ct.select = unique(sc_meta$cellType),
    sample.varname = "sampleInfo",
    minCountGene = 100,
    minCountSpot = 5) 
  
  ### Run CARD
  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
  
  ### Single cell resolution mapping
  # Note that here the shapeSpot is the user defined variable which indicates the capturing area of single cells. Details see above.
  scMapping = CARD_SCMapping(CARD_obj,shapeSpot="Circle",numCell=10,ncore=40)
  
  ### Spatial location info and expression count of the single cell resolution data
  MapCellCords = as.data.frame(colData(scMapping))
  count_SC = assays(scMapping)$counts
  
  ### Define colors according to cell types
  malig_colours <- c("lightcoral", "darkgoldenrod2", "dodgerblue3")[1:length(grep("Malig", unique(MapCellCords$CT), value = T))]
  names(malig_colours) <- grep("Malig", unique(MapCellCords$CT), value = T)
  cell_colours <- c("Endothelial" = "mediumorchid3", 
                    "Macrophage" = "skyblue2",
                    "Skeletal Muscle" = "darkslateblue",
                    "T cell" = "seagreen4", malig_colours)
  
  ### Make the spatial plot
  pdf(file = paste0("figure_4A_", sample, ".pdf"), height = 7, width = 10)
  print(MapCellCords %>% ggplot(aes(x = -x, y = -y, colour = CT)) + geom_point(size = 1.0) + 
          scale_color_manual(values = cell_colours) +
          theme(plot.margin = margin(0, 0, 0, 0, "cm"),
                panel.background = element_rect(colour = NA, fill="transparent"),
                plot.background = element_rect(colour = NA, fill="transparent"),
                legend.position="right",
                axis.text =element_blank(),
                axis.ticks =element_blank(),
                axis.title =element_blank(),
                legend.title=element_text(size = 16,face="bold"),
                legend.text=element_text(size = 16),
                legend.key = element_rect(colour = NA, fill = "transparent"),
                legend.key.size = unit(0.45, 'cm'),
                strip.text = element_text(size = 15,face="bold"))+
          guides(color=guide_legend(title="Cell Type", override.aes = list(size=8), ncol = 1)))
  dev.off()
}


