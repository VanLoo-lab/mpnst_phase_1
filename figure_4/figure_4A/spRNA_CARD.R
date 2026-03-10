#ml rgdal/1.5-16-foss-2020a-R-4.0.0
#ml R/4.0.0-foss-2020a
#R

# options(bitmapType='cairo')

library(tidyverse)
library(CARD) 
# Had to install CARD in home folder due to V8
# Note also had to install older version of MuSiC (used https://github.com/omnideconv/MuSiC which is older version forked from xuranw/MuSiC)

# library(tidyverse, lib.loc = "/camp/lab/vanloop/working/yanh/R/library_4.0.0/")
# library(Seurat, lib.loc = "/camp/lab/vanloop/working/yanh/R/library_4.0.0/")
samples = c("Pa", "R1", "R3", "Pb", "R2", "R4", "R5a", "R5b")
names(samples) = c("VER683A1", "VER683A2", "VER683A3", "VER683A4", "VER683A5", "VER683A6", "VER683A7", "VER683A8")

OUTPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_4/results/figure_4A/"
INPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_4/data/"

output.dir <- OUTPUTDIR
setwd(output.dir)

####################################################################################################################################
### Part 2: 26/05/23 Run on MPNST
####################################################################################################################################

# The code is unable to overwrite files and errors if files already exist

for (sample in samples[1:8]) {
  #Load sc and sp data
  sc_count <- readRDS(paste0(INPUTDIR, sample, "_sc_count.rds"))
  sc_meta <- readRDS(paste0(INPUTDIR, sample, "_sc_meta.rds"))
  sp_count <- readRDS(paste0(INPUTDIR, sample, "_sp_count.rds"))
  sp_location <- readRDS(paste0(INPUTDIR, sample, "_sp_location.rds"))
  
  #Create CARD object
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
  
  #Run CARD
  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
  #saveRDS(CARD_obj, paste0(sample,"_CARD_obj.rds"))
  
  # Alex Stein: While this object existed in the data of Haixi, it was inconsistent
  # with the code above creating CARD_obj from the other data.
  #CARD_obj <- readRDS(paste0(INPUTDIR, sample,"_CARD_obj.rds"))
  
  ##Visualise cell type
  ## set the colors. Here, I just use the colors in the manuscript, if the color is not provided, the function will use default color in the package. 
  # colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
  #            "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
  #            "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
  p1 <- CARD.visualize.pie(proportion = CARD_obj@Proportion_CARD,spatial_location = CARD_obj@spatial_location)#, colors = colors)
  pdf(file = paste0(sample, "_CARD_Proportion.pdf"), height = 14, width = 14)
  # png(filename = paste0(sample, "_CARD_Proportion.png"), height = 2000, width = 2000, res = 200)
  print(p1)
  dev.off()
  
  ## select the cell type that we are interested
  ct.visualize = as.character(unique(sc_meta$cellType))
  ## visualize the spatial distribution of the cell type proportion
  p2 <- CARD.visualize.prop(
    proportion = CARD_obj@Proportion_CARD,        
    spatial_location = CARD_obj@spatial_location, 
    ct.visualize = ct.visualize[],                 ### selected cell types to visualize
    # colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
    NumCols = 4)                                 ### number of columns in the figure panel
  pdf(file = paste0(sample, "_CARD_cell_type.pdf"), height = 14, width = 28)
  print(p2)
  dev.off()
  
  ##Check against genotyping information
  # tumour_normal_prop =  data.frame(Tumour = rowSums(CARD_obj@Proportion_CARD[,grepl("Malignant", colnames(CARD_obj@Proportion_CARD)), drop = F]),
  #                                  Normal = rowSums(CARD_obj@Proportion_CARD[,!grepl("Malignant", colnames(CARD_obj@Proportion_CARD)), drop = F])) %>%
  #   rownames_to_column("Barcode") %>% mutate(Barcode = gsub("-1", "", paste0(sample, "_", Barcode)))
  #
  # both_haplo_count_combined_by_cell <- readRDS("../Genotype_SNPs/both_haplo_count_combined_by_cell.rds") %>% 
  #   filter(str_detect(Barcode, sample)) %>%
  #   mutate(Total_Count = Haplo_1_Count + Haplo_2_Count, Haplo_2_pct = Haplo_2_Count/Total_Count) %>% left_join(tumour_normal_prop, by = "Barcode")
  # 
  # png(filename = paste0(sample, "_haplo_counts_log_vs_prop.png"), width = 4000, height = 4000, res = 200)
  # print(both_haplo_count_combined_by_cell %>% filter(Total_Count >=50) %>%
  #         ggplot(aes(x=Haplo_2_pct, y=Normal)) + geom_point(color = "#F8766D") + geom_smooth(method=lm, level = 0.99) +
  #         xlim(0,0.5) + ylim(0,1) +
  #         xlab("Haplotype 2 Proportion") + ylab("Normal Proportion") +
  #         ggtitle(sample) + theme(plot.title = element_text(size = 40), text=element_text(size=30), aspect.ratio=1))
  # dev.off()
  # 
  # png(filename = paste0(sample, "_haplo_counts_log_prop.png"), width = 4000, height = 4000, res = 200)
  # print(both_haplo_count_combined_by_cell %>% 
  #         ggplot(aes(x=Haplo_1_Count, y=Haplo_2_Count, colour = Tumour)) + geom_point() + 
  #         scale_x_log10(limits = c(1,1000)) + scale_y_log10(limits = c(1,1000)) +
  #         xlab("Total number of counts of haplotype 1") + ylab("Total number of counts of haplotype 2") +
  #         ggtitle(sample) + theme(plot.title = element_text(size = 40), text=element_text(size=30), aspect.ratio=1))
  # dev.off()
  
  ## visualize the spatial distribution of two cell types on the same plot
  if (length(grep("Malignant", ct.visualize, value = T)) < 2) {
    ct2.visualize = c(grep("Malignant", ct.visualize, value = T), "Macrophage")
  } else {
    ct2.visualize = grep("Malignant", ct.visualize, value = T)[1:2]
  }
  
  p3 = CARD.visualize.prop.2CT(
    proportion = CARD_obj@Proportion_CARD,                             ### Cell type proportion estimated by CARD
    spatial_location = CARD_obj@spatial_location,                      ### spatial location information
    ct2.visualize = ct2.visualize)             ### two cell types you want to visualize
  # colors = list(c("lightblue","lightyellow","red"),c("lightblue","lightyellow","black")))       ### two color scales                             
  pdf(file = paste0(sample, "_CARD_2_cell_type.pdf"), height = 7, width = 7)
  print(p3)
  dev.off()
  
  ## visualise the cell proportion correlation
  p4 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
  pdf(file = paste0(sample, "_CARD_cell_prop_corr.pdf"), height = 7, width = 7)
  print(p4)
  dev.off()
  
  ## Single cell resolution mapping
  ## Note that here the shapeSpot is the user defined variable which indicates the capturing area of single cells. Details see above.
  scMapping = CARD_SCMapping(CARD_obj,shapeSpot="Circle",numCell=10,ncore=40)
  print(scMapping)
  saveRDS(scMapping, paste0(sample, "_CARD_scMapping.rds"))
  
  scMapping <- readRDS(paste0(sample, "_CARD_scMapping.rds"))
  ### spatial location info and expression count of the single cell resolution data
  library(SingleCellExperiment)
  MapCellCords = as.data.frame(colData(scMapping))
  count_SC = assays(scMapping)$counts
  
  ## visualise single cell resolution mapping
  # colors = c("#8DD3C7","#CFECBB","#F4F4B9","#CFCCCF","#D1A7B9","#E9D3DE","#F4867C","#C0979F",
  #            "#D5CFD6","#86B1CD","#CEB28B","#EDBC63","#C59CC5","#C09CBF","#C2D567","#C9DAC3","#E1EBA0",
  #            "#FFED6F","#CDD796","#F8CDDE")
  p10 = ggplot(MapCellCords, aes(x = x, y = y, colour = CT)) + 
    geom_point(size = 1.0) +
    # scale_colour_manual(values =  colors) +
    #facet_wrap(~Method,ncol = 2,nrow = 3) + 
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          panel.background = element_rect(colour = "white", fill="white"),
          plot.background = element_rect(colour = "white", fill="white"),
          legend.position="bottom",
          panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
          axis.text =element_blank(),
          axis.ticks =element_blank(),
          axis.title =element_blank(),
          legend.title=element_text(size = 13,face="bold"),
          legend.text=element_text(size = 12),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.key.size = unit(0.45, 'cm'),
          strip.text = element_text(size = 15,face="bold"))+
    guides(color=guide_legend(title="Cell Type"))
  pdf(file = paste0(sample, "_CARD_single_cell.pdf"), height = 7, width = 7)
  print(p10)
  dev.off()
  
  malig_colours <- c("lightcoral", "darkgoldenrod2", "dodgerblue3")[1:length(grep("Malig", unique(MapCellCords$CT), value = T))]
  names(malig_colours) <- grep("Malig", unique(MapCellCords$CT), value = T)
  cell_colours <- c("Endothelial" = "mediumorchid3", 
                    "Macrophage" = "skyblue2",
                    "Skeletal Muscle" = "darkslateblue",
                    "T cell" = "seagreen4", malig_colours)
  
  pdf(file = paste0(sample, "_CARD_single_cell_col.pdf"), height = 7, width = 10)
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


