options(bitmapType='cairo') #to solve plotting issue on CAMP
library(tidyverse, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(ggpubr, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(gridExtra, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(VariantAnnotation, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(liftOver, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(copynumber, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(ComplexHeatmap, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(pbapply, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(pbmcapply, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(dendextend, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(grid, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(gridExtra, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(rdist, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")
library(RColorBrewer, lib.loc = "/camp/lab/vanloop/working/yanh/R/library/")


#Load metadata
source(file = "/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/code/10X_DLP_metadata.R")

output.dir = paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/results/snv_mnv_indel_tumouronly/SNV_analysis/")
system(paste0("mkdir -p ", output.dir))
setwd(output.dir)

####################################################################################################################################
#Choose which clustering run to use
####################################################################################################################################
filter_n_cells <- 2 #SNVs with less than or equal to filter_n_cells cell removed 
sub_folder <- paste0(run, "_filter_n_", filter_n_cells)
system(paste0("mkdir -p ", output.dir, sub_folder))

#Use no cells as proxy of coverage instead
scDNA_CN_cluster_ids <- scDNA_CN_cluster_ids[str_sort(names(scDNA_CN_cluster_ids))]
# scDNA_CN_clusters_coverage <- unlist(lapply(scDNA_CN_cluster_ids, function(k) {
#   return(length(k)/15) #15 is just a scaling factor
# }))

# #Calculate coverage for each cluster of cells (region pseudobulk coverage * no cells in cluster / total cells from region)
# region_10X_coverage <- readRDS("../../../../10X_DNA/results/Genotype_de_novo_SNVs/MPNST_region_coverage.rds")
# region_DLP_coverage <- readRDS("../../../../DLP_plus/results/Genotype_SNVs_CICC/MPNST_DLP_region_coverage.rds")
# 
# scDNA_CN_clusters_coverage <- unlist(lapply(1:length(scDNA_CN_cluster_ids), function(k) {
#   tenX_coverage <- lapply(samples, function(s) {
#     num_cells <- sum(grepl(paste0(s, "_10X"), scDNA_CN_cluster_ids[[k]]))
#     cells_region <- sum(grepl(paste0(s, "_10X"), unlist(scDNA_CN_cluster_ids)))
#     return(region_10X_coverage[[names(samples[samples == s])]] * (num_cells / cells_region))
#   })
#   DLP_coverage <- lapply(samples, function(s) {
#     num_cells <- sum(grepl(paste0(s, "_DLP"), scDNA_CN_cluster_ids[[k]]))
#     cells_region <- sum(grepl(paste0(s, "_DLP"), unlist(scDNA_CN_cluster_ids)))
#     return(region_DLP_coverage[[s]] * (num_cells / cells_region))
#   })
#   return(sum(unlist(tenX_coverage)) + sum(unlist(DLP_coverage)))
# }))

####################################################################################################################################
#Functions for clustering SNVs and plotting
####################################################################################################################################
if (T) {
  cluster_plot_SNVs <- function(SNV_set = "de_novo", G1000 = F, upper_limit = NULL, ex_clonal = F, J = 36, ann_cluster = T, plot = T, reorder = NULL) {
    if (G1000 == T) {G1000F = "G1000filtered_"} else {G1000F = NULL}
    if (SNV_set == "DLP") {
      scDNA_CN_clusters_SNVs <- readRDS(paste0(sub_folder, "/MPNST_",SNV_set,"_scDNA_",G1000F,"CN_clusters_SNVs.rds"))
      sub_suffix <- NULL
    }
    if (SNV_set == "de_novo") {
      scDNA_CN_clusters_SNVs <- readRDS(paste0(sub_folder, "/MPNST_",SNV_set,"_scDNA_",G1000F,"CN_clusters_SNVs.rds"))
      sub_suffix <- NULL
    }
    if (SNV_set == "all" | SNV_set == "bulk") {
      scDNA_CN_clusters_de_novo_SNVs <- readRDS(paste0(sub_folder, "/MPNST_de_novo_scDNA_",G1000F,"CN_clusters_SNVs.rds"))
      scDNA_CN_clusters_bulk_SNVs <- readRDS(paste0(sub_folder, "/MPNST_scDNA_CN_clusters_bulk_SNVs.rds"))
      CICC_output <- do.call(rbind, readRDS("../../../../10X_DNA/results/Genotype_SNVs/MPNST_SNVs_clusters.rds"))
      if (ex_clonal) {
        sub_suffix <- "sub_"
        remove_clusters <- "XXXXXX"
        sub_suffix <- paste0(sub_suffix, "clonal_")
        # CICC_output <- readRDS("../../../../bulk/results/DPClust/CICC_ex1_ex2/MPNST_1_consensus_clusters_40.rds")
        # CICC_output$chr <- paste0("chr", CICC_output$chr)
        SNVs_to_keep <- CICC_output %>% filter(!cluster %in% remove_clusters) %>% mutate(ID = paste0(chr, "_", pos)) %>% pull(ID) #remove rogue clusters
      } else {
        SNVs_to_keep <- rownames(scDNA_CN_clusters_bulk_SNVs)
        sub_suffix <- NULL
      }
      if (SNV_set == "bulk") {
        scDNA_CN_clusters_SNVs <- rbind(scDNA_CN_clusters_bulk_SNVs[rownames(scDNA_CN_clusters_bulk_SNVs) %in% SNVs_to_keep,])
      } else {
        scDNA_CN_clusters_bulk_SNVs_keep <- scDNA_CN_clusters_bulk_SNVs[rownames(scDNA_CN_clusters_bulk_SNVs) %in% SNVs_to_keep,]
        scDNA_CN_clusters_SNVs <- rbind(scDNA_CN_clusters_bulk_SNVs_keep, 
                                        scDNA_CN_clusters_de_novo_SNVs[which(!rownames(scDNA_CN_clusters_de_novo_SNVs) %in% rownames(scDNA_CN_clusters_bulk_SNVs_keep)),])
      }
    }
    
    if (!is.null(upper_limit)) {
      scDNA_CN_clusters_SNVs[scDNA_CN_clusters_SNVs > upper_limit] <- upper_limit
      #Need to removes where values are identical
      scDNA_CN_clusters_SNVs <- scDNA_CN_clusters_SNVs[rowSums(scDNA_CN_clusters_SNVs) != ncol(scDNA_CN_clusters_SNVs) * upper_limit,]
      upper_suffix <- paste0("_upper_", upper_limit)
    } else {upper_suffix <- NULL}
    
    if (!file.exists(paste0(sub_folder, "/MPNST_",SNV_set,"_",sub_suffix,"scDNA_SNV_",G1000F,"per_CN_cluster_hclust_cor_ward_ex_WGD_N",upper_suffix,".rds"))) {
      print("Clustering")
      saveRDS(hclust(as.dist(1-cor(t(scDNA_CN_clusters_SNVs))), method = "ward.D2"),
              paste0(sub_folder, "/MPNST_",SNV_set,"_",sub_suffix,"scDNA_SNV_",G1000F,"per_CN_cluster_hclust_cor_ward_ex_WGD_N",upper_suffix,".rds"))#Use correlation distance
      MPNST_CN_clusters_SNVs_hclust_cor_ward <- readRDS(paste0(sub_folder, "/MPNST_",SNV_set,"_",sub_suffix,"scDNA_SNV_",G1000F,"per_CN_cluster_hclust_cor_ward_ex_WGD_N",upper_suffix,".rds"))
      png(filename = paste0(sub_folder, "/MPNST_",SNV_set,"_",sub_suffix,"SNV_",G1000F,"per_CN_cluster_hclust_cor_ward_ex_WGD_N_dendrogram",upper_suffix,".png"), width = 4000, height = 2000, res = 200)
      plot(rev(MPNST_CN_clusters_SNVs_hclust_cor_ward), main = "Dendrogram - Correlation Ward.D2", xlab = "", sub = "", labels = F)
      dev.off()
    } 
    MPNST_CN_clusters_SNVs_hclust_cor_ward <- readRDS(paste0(sub_folder, "/MPNST_",SNV_set,"_",sub_suffix,"scDNA_SNV_",G1000F,"per_CN_cluster_hclust_cor_ward_ex_WGD_N",upper_suffix,".rds"))
    
    if (plot) {
      print("Plotting")
      #Reorder in size of coverage and plot to proportion of cells
      CNA_cluster_resize <- unlist(lapply(1:length(scDNA_CN_cluster_ids), function(x) { #For each cluster
        return(rep(x, ifelse(round(length(scDNA_CN_cluster_ids[[x]])/30) < 1, 1, round(length(scDNA_CN_cluster_ids[[x]])/30))))
      }))
      
      if (SNV_set == "de_novo" | SNV_set =="DLP" | SNV_set == "bulk") {
        if (is.null(upper_limit)) {
          # scDNA_mtx_breaks <- c(-0.5, filter_n_cells+0.5, filter_n_cells+1.5, filter_n_cells+2.5, filter_n_cells+3.5)
          SNV_counts <- c(0:filter_n_cells, sort(unique(as.numeric(scDNA_CN_clusters_SNVs)))[sort(unique(as.numeric(scDNA_CN_clusters_SNVs))) > filter_n_cells])
          color_scale <- c(rep("gray80", filter_n_cells+1), colorRampPalette(c("gray80", "darkred"))(4)[-1], rep("darkred",length(SNV_counts)-filter_n_cells-4))
          names(color_scale) <- SNV_counts
        }
      }
      
      if (SNV_set == "all") {
        # scDNA_mtx_breaks <- c(-filter_n_cells-3.5, -filter_n_cells-2.5, -filter_n_cells-1.5, -filter_n_cells-.5, 1.5, 2.5, 3.5, 4.5)
        # color_scale <- colorRampPalette(c("darkblue", "gray80", "darkred"))(length(scDNA_mtx_breaks)-1)
        if (is.null(upper_limit)) {
          SNV_counts <- c(0:filter_n_cells, sort(unique(as.numeric(scDNA_CN_clusters_SNVs)))[sort(unique(as.numeric(scDNA_CN_clusters_SNVs))) > filter_n_cells])
          color_scale_bulk <- c(rep("gray80", filter_n_cells), colorRampPalette(c("gray80", "darkred"))(4)[-1], rep("darkred",length(SNV_counts)-filter_n_cells-4))
          color_scale_denovo <- c(rep("darkblue",length(SNV_counts)-filter_n_cells-4), colorRampPalette(c("darkblue", "gray80"))(4)[-1], rep("gray80", filter_n_cells))
          color_scale <- c(color_scale_denovo, "gray80", color_scale_bulk)
          names(color_scale) <- c(-c(sort(unique(as.numeric(scDNA_CN_clusters_SNVs)), decreasing = T)[sort(unique(as.numeric(scDNA_CN_clusters_SNVs)), decreasing = T) > filter_n_cells], filter_n_cells:0),
                                  c(1:filter_n_cells, sort(unique(as.numeric(scDNA_CN_clusters_SNVs)))[sort(unique(as.numeric(scDNA_CN_clusters_SNVs))) > filter_n_cells]))
        }
        #After clustering convert de novo SNVs into negative for visibility on plot
        scDNA_CN_clusters_SNVs <- rbind(scDNA_CN_clusters_bulk_SNVs, 
                                        -scDNA_CN_clusters_de_novo_SNVs[which(!rownames(scDNA_CN_clusters_de_novo_SNVs) %in% rownames(scDNA_CN_clusters_bulk_SNVs)),])
      }
      
      if (!is.null(upper_limit)) {
        # scDNA_mtx_breaks <- c(-2.5, -1.5, 1.5, 2.5)
        # color_scale <- colorRampPalette(c("darkblue", "gray80", "darkred"))(length(scDNA_mtx_breaks)-1)
        color_scale <- c("darkblue", "gray80", "darkred") #binary colour scheme
        names(color_scale) <- c(-2,0,2)
        #Need to removes where values are identical
        scDNA_CN_clusters_SNVs <- scDNA_CN_clusters_SNVs[MPNST_CN_clusters_SNVs_hclust_cor_ward[["labels"]],]
        scDNA_CN_clusters_SNVs[scDNA_CN_clusters_SNVs > upper_limit] <- upper_limit #apply upper limit
        scDNA_CN_clusters_SNVs[scDNA_CN_clusters_SNVs < -upper_limit] <- -upper_limit
      }
      
      scDNA_CN_clusters_SNVs_resize <- scDNA_CN_clusters_SNVs[,CNA_cluster_resize]
      colnames(scDNA_CN_clusters_SNVs_resize) <- colnames(scDNA_CN_clusters_SNVs)[CNA_cluster_resize]
      scDNA_CN_clusters_SNVs_resize <- t(scDNA_CN_clusters_SNVs_resize) #plot horizontal
      saveRDS(scDNA_CN_clusters_SNVs_resize, paste0(sub_folder, "/MPNST_",SNV_set,"_",sub_suffix,"SNV_",G1000F,"clusters_ex_WGD_N_mtx.rds"))
      
      #Annotate
      # annotate_row <- data.frame("CNA_Region" = colnames(scDNA_CN_clusters_SNVs))
      #                            #"CNA_Cluster" = as.factor(scDNA_CN_cluster_orig_names))
      # rownames(annotate_row) <- colnames(scDNA_CN_clusters_SNVs)
      # region_colours <- c("#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3", "#F8766D")
      # names(region_colours) <- samples
      # CNA_Region_col = region_colours[gsub("_.*", "", colnames(scDNA_CN_clusters_SNVs))]
      # names(CNA_Region_col) = annotate_row$CNA_Region
      # CNA_Cluster_col = hcl(h = seq(15, 375, length = as.numeric(gsub("_.*", "", sub_folder)) + 1), l = 65, c = 100)[scDNA_CN_cluster_orig_names]
      # names(CNA_Cluster_col) = annotate_row$CNA_Cluster
      # if (SNV_set == "bulk") {
      #   annotate_col <- data.frame("SNV_name" = MPNST_CN_clusters_SNVs_hclust_cor_ward[["labels"]]) %>% 
      #     left_join(CICC_output %>% mutate(SNV_name = paste0(chr, "_", pos), SNV_cluster = as.character(cluster))) %>% dplyr::select(SNV_cluster)
      #   rownames(annotate_col) <- MPNST_CN_clusters_SNVs_hclust_cor_ward[["labels"]]
      #   SNV_Cluster_col <- hcl(h = seq(15, 375, length = length(unique(annotate_col$SNV_cluster)) + 1), l = 65, c = 100)[-1]
      #   names(SNV_Cluster_col) <- unique(annotate_col$SNV_cluster)
      #   ann_colors = list(CNA_Region = CNA_Region_col, CNA_Cluster = CNA_Cluster_col, SNV_cluster = SNV_Cluster_col)
      # } else{
      #   annotate_col = NA
      #   ann_colors = list(CNA_Region = CNA_Region_col, CNA_Cluster = CNA_Cluster_col)
      # }
      ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_clusters_SNVs_resize))),
                             col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D")),
                             annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 20), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), ncol = 1)), show_annotation_name = F)
      if (ann_cluster) {
        ha_col = columnAnnotation(df = data.frame("SNV_name" = MPNST_CN_clusters_SNVs_hclust_cor_ward[["labels"]]) %>% 
                                    left_join(CICC_output %>% mutate(SNV_name = paste0(chr, "_", pos), SNV_cluster = as.character(cluster))) %>% dplyr::mutate(SNV_cluster = replace_na(SNV_cluster, "De novo")) %>% dplyr::select(SNV_cluster),
                                  col = list(SNV_cluster = c("XXXXXX" = "black", "00000X" = "#F8766D", "000003" = "pink", "X000X0" = "coral", "X00000" = "#B79F00", "300000" = "gold", "0000X0" = "#F564E3", "000030" = "mediumpurple1",
                                                             "0XXX00" = "darkslategray4", "0X0X00" = "darkgreen", "0X0000" = "#00BA38", "030000" = "chartreuse", "000X00" = "#619CFF", "000400" = "deepskyblue", "00X000" = "#00BFC4", "003000" = "cyan", "De novo" = "Gray50")),
                                  annotation_legend_param = list(SNV_cluster = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 20), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), ncol = 1)), show_annotation_name = F)
      } else {
        ha_col = NULL
      }
      
      png(filename = paste0(sub_folder, "/MPNST_",SNV_set,"_",sub_suffix,"SNV_",G1000F,"per_CN_cluster_ex_WGD_N_heatmap_resize",upper_suffix,".png"), width = 4000, height = 4000, res = 200)
      print(Heatmap(scDNA_CN_clusters_SNVs_resize, col = color_scale, cluster_rows = F, row_split = names(scDNA_CN_cluster_ids)[CNA_cluster_resize], cluster_columns = MPNST_CN_clusters_SNVs_hclust_cor_ward, 
                    row_title_rot = 0, show_row_names = F, show_column_names = F, row_title_gp = grid::gpar(fontsize = 24),
                    left_annotation = ha_row, top_annotation = ha_col, show_heatmap_legend = F))
      # pheatmap(scDNA_CN_clusters_SNVs_resize, color = color_scale, breaks = scDNA_mtx_breaks, cluster_rows = F, cluster_cols = MPNST_CN_clusters_SNVs_hclust_cor_ward,
      #          annotation_col = annotate_col, annotation_row = annotate_row, annotation_colors = ann_colors, show_rownames = F, show_colnames = F)
      dev.off()
      
      #Get SNV clusters after cutting and plot
      if (!is.null(J)) {
        if (!file.exists(paste0(sub_folder, "/MPNST_",SNV_set,"_",sub_suffix,"SNV_",G1000F,"clusters_ex_WGD_N_ids_",J,".rds"))) {
          CN_based_SNV_clusters <- cutree(MPNST_CN_clusters_SNVs_hclust_cor_ward, k=J, order_clusters_as_data = F)
          as.vector(table(CN_based_SNV_clusters)) #nice just to see how big clusters are
          CN_based_SNV_clusters_ids <- list()
          for (i in 1:J) {
            CN_based_SNV_clusters_ids[[i]] <- names(CN_based_SNV_clusters[CN_based_SNV_clusters==i])
          }
          saveRDS(CN_based_SNV_clusters_ids, paste0(sub_folder, "/MPNST_",SNV_set,"_",sub_suffix,"SNV_",G1000F,"clusters_ex_WGD_N_ids_",J,".rds"))
        }
        CN_based_SNV_clusters <- cutree(MPNST_CN_clusters_SNVs_hclust_cor_ward, k=J, order_clusters_as_data = F)
        CN_based_SNV_clusters_ids <- readRDS(paste0(sub_folder, "/MPNST_",SNV_set,"_",sub_suffix,"SNV_",G1000F,"clusters_ex_WGD_N_ids_",J,".rds"))
        
        #Plot cut
        png(filename = paste0(sub_folder, "/MPNST_",SNV_set,"_",sub_suffix,"SNV_",G1000F,"per_CN_cluster_ex_WGD_N_heatmap_resize",upper_suffix,"_", J, ".png"), width = 4000, height = 4000, res = 200)
        print(Heatmap(scDNA_CN_clusters_SNVs_resize, col = color_scale, cluster_rows = F, row_split = names(scDNA_CN_cluster_ids)[CNA_cluster_resize], 
                      cluster_columns = MPNST_CN_clusters_SNVs_hclust_cor_ward, column_split = J, column_title = NULL, show_column_dend = F,
                      cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                      row_title_rot = 0, show_row_names = F, show_column_names = F, row_title_gp = grid::gpar(fontsize = 24), left_annotation = ha_row, top_annotation = ha_col, show_heatmap_legend = F))
        # pheatmap(scDNA_CN_clusters_SNVs_resize, color = color_scale, breaks = scDNA_mtx_breaks, cutree_rows = J, cluster_rows = F, cluster_cols = MPNST_CN_clusters_SNVs_hclust_cor_ward, 
        #          annotation_col = annotate_col, annotation_colors = ann_colors, show_rownames = F, show_colnames = F)
        dev.off()
        
        #Plot reordered
        if (!is.null(reorder)) {
          #Reorder SNV clusters in order provided
          CN_based_SNV_clusters_ids_reorder <- CN_based_SNV_clusters_ids[reorder]
          SNV_cluster_reorder <- sapply(colnames(scDNA_CN_clusters_SNVs_resize), function(snv) grep(snv, CN_based_SNV_clusters_ids_reorder)) %>% unname()
          # SNV_cluster_reorder <- CN_based_SNV_clusters_ids[reorder]
          
          # pdf(file = paste0(sub_folder, "/MPNST_",SNV_set,"_",sub_suffix,"SNV_",G1000F,"per_CN_cluster_ex_WGD_N_heatmap_resize",upper_suffix,"_", J, "_reordered.pdf"), width = 14, height = 14)
          png(filename = paste0(sub_folder, "/MPNST_",SNV_set,"_",sub_suffix,"SNV_",G1000F,"per_CN_cluster_ex_WGD_N_heatmap_resize",upper_suffix,"_", J, "_reordered.png"), width = 4000, height = 4000, res = 200)
          print(Heatmap(scDNA_CN_clusters_SNVs_resize, col = color_scale, cluster_rows = F, row_split = names(scDNA_CN_cluster_ids)[CNA_cluster_resize], 
                        cluster_columns = F, column_split = SNV_cluster_reorder, show_column_dend = F, column_title = NULL, 
                        cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                        row_title_rot = 0, show_row_names = F, show_column_names = F, row_title_gp = grid::gpar(fontsize = 24), left_annotation = ha_row, top_annotation = ha_col, show_heatmap_legend = F))
          # pheatmap(scDNA_CN_clusters_SNVs_resize[,unlist(SNV_cluster_reorder)], color = color_scale, breaks = scDNA_mtx_breaks, cluster_rows = F, cluster_cols = F, 
          #          gaps_row = cumsum(unlist(lapply(SNV_cluster_reorder, length))), annotation_col = annotate_col, annotation_colors = ann_colors, show_rownames = F, show_colnames = F)
          dev.off()
        }
      }
    }
  }
}

####################################################################################################################################
###Try to cluster SNVs of single cells based on CN (WITHOUT WGD AND NORMAL CELLS)
####################################################################################################################################
#Create SNV matrix for de novo calls
if (T) {
  #Load in SNVs of all cells and subset to pass filter cells
  de_novo_mt_10X_counts_SNV_vec_cell <- readRDS("../../../../10X_DNA/results/Genotype_de_novo_SNVs/MPNST_scDNA_mt_counts_SNV_vec_cell.rds")
  names(de_novo_mt_10X_counts_SNV_vec_cell) <- unlist(barcodes)
  
  de_novo_mt_DLP_counts_SNV_vec_cell <- readRDS("../../../../DLP_plus/results/snv_mnv_indel_tumouronly/SNV_analysis/MPNST_10X_scDNA_mt_counts_SNV_vec_cell.rds")
  names(de_novo_mt_DLP_counts_SNV_vec_cell) <- all_DLP_barcodes
  
  de_novo_mt_counts_SNV_vec_cell <- lapply(c(de_novo_mt_10X_counts_SNV_vec_cell, de_novo_mt_DLP_counts_SNV_vec_cell), function(x) {
    x %>% arrange(chr, pos) #make sure in same order
  })

  #For each CN cluster get number of SNVs
  scDNA_CN_clusters_SNVs <- do.call(cbind, lapply(1:length(scDNA_CN_cluster_ids), function(x) {
    SNV_vec_cells <- de_novo_mt_counts_SNV_vec_cell[scDNA_CN_cluster_ids[[x]]]
    SNV_vec_cells_sum <- do.call(cbind, lapply(SNV_vec_cells, function(c) {
      return(c$Present)
    })) %>% rowSums()
    return(SNV_vec_cells_sum)
  }))
  rownames(scDNA_CN_clusters_SNVs) <- paste0(de_novo_mt_counts_SNV_vec_cell[[1]][,"chr"], "_", de_novo_mt_counts_SNV_vec_cell[[1]][,"pos"])
  
  #Apply a filter to remove SNVs in just n cells
  scDNA_CN_clusters_SNVs[scDNA_CN_clusters_SNVs <= filter_n_cells] <- 0
  
  #Have to remove SNVs with no calls in all clusters as std dev is 0 and correlation distance doesn't work (1366 SNVs left)
  scDNA_CN_clusters_SNVs <- scDNA_CN_clusters_SNVs[rowSums(scDNA_CN_clusters_SNVs) != 0,]
  
  #Annotate CNA names
  colnames(scDNA_CN_clusters_SNVs) <- names(scDNA_CN_cluster_ids)
  
  saveRDS(scDNA_CN_clusters_SNVs, paste0(sub_folder, "/MPNST_de_novo_scDNA_CN_clusters_SNVs.rds"))
  # 
  # #Keep just filtered SNVs
  # mt_counts_SNV_mtx_filter <- readRDS("../../../../10X_DNA/results/snv_mnv_indel_tumouronly/SNV_analysis/MPNST_mt_counts_SNV_mtx_filter.rds")
  # scDNA_CN_clusters_SNVs <- scDNA_CN_clusters_SNVs[rownames(scDNA_CN_clusters_SNVs) %in% colnames(mt_counts_SNV_mtx_filter),]
  # 
  # 
  # #Remove just those with >2 read in bulk normal
  # if (!file.exists("DLP_PN_vcf_SNVs_af.rds")) {
  #   #Apply filter to remove de novo SNVs with >1 reads in normal
  #   if (F) {
  #     NBCORES = 8
  #     MIN_BASE_QUAL = 20
  #     MIN_MAP_QUAL = 0
  #     
  #     getAlleleCounts <- function (bam.file, output.file, g1000.loci, min.base.qual = 20,
  #                                  min.map.qual = 35, allelecounter.exe = "alleleCounter") {
  #       cmd = paste(allelecounter.exe, "-b", bam.file, "-l", g1000.loci,
  #                   "-o", output.file, "-m", min.base.qual, "-q", min.map.qual, "-f 0 -F 0")
  #       system(cmd, wait = T)
  #     }
  #     # system("module load alleleCount") doesn't work in R studio
  #     ALLELECOUNTER = "/camp/apps/eb/software/alleleCount/4.0.0-foss-2016b/bin/alleleCounter"
  #     
  #     getAlleleCounts(bam.file=paste0("../../../../bulk/data/merged_bam/VER236A7_merged_rmdup_recal.bam"),
  #                     output.file=paste0("scDNA_Allele_Freq/bulk_normal_alleleFrequencies_all.txt"),
  #                     g1000.loci=paste0("SNV_loci/DLP_SNV_loci.txt"),
  #                     min.base.qual=MIN_BASE_QUAL,
  #                     min.map.qual=MIN_MAP_QUAL,
  #                     allelecounter.exe=ALLELECOUNTER)
  #   }
  #   
  #   bulk_normal_allele_freq_all <- read_tsv(paste0("scDNA_Allele_Freq/bulk_normal_alleleFrequencies_all.txt")) %>% mutate(name = paste0(`#CHR`, "_", POS))
  #   DLP_PN_vcf_SNVs <- do.call(rbind, readRDS("DLP_PN_vcf_SNVs.rds")) %>% 
  #     mutate(ref = gsub("/.", "", gsub(".*_", "", vcf_name)), alt = gsub(".*/", "", vcf_name)) %>% 
  #     left_join(bulk_normal_allele_freq_all, by = "name") 
  #   DLP_PN_vcf_SNVs_af <- do.call(rbind, pbmclapply(1:nrow(DLP_PN_vcf_SNVs), function(v) {
  #     return(DLP_PN_vcf_SNVs[v,] %>% mutate(ref_count = get(paste0("Count_", ref)), alt_count = get(paste0("Count_", alt))))
  #   }, mc.cores = 60))
  #   saveRDS(DLP_PN_vcf_SNVs_af, "DLP_PN_vcf_SNVs_af.rds")
  # } else {
  #   DLP_PN_vcf_SNVs_af <- readRDS("DLP_PN_vcf_SNVs_af.rds")
  # }
  # SNVs_in_bulk_normal <- DLP_PN_vcf_SNVs_af %>% filter(alt_count > 2) %>% pull(name)
  # scDNA_CN_clusters_SNVs <- scDNA_CN_clusters_SNVs[!rownames(scDNA_CN_clusters_SNVs) %in% SNVs_in_bulk_normal, ]
  # 
  # #Remove bulk calls (only for DLP runs)
  # CICC_output <- do.call(rbind, readRDS("../../Genotype_SNVs_CICC/MPNST_SNVs_clusters.rds"))
  # scDNA_CN_clusters_SNVs <- scDNA_CN_clusters_SNVs[!rownames(scDNA_CN_clusters_SNVs) %in% paste0(CICC_output$chr, "_", CICC_output$pos),]
  # 
  # #Apply a filter to remove SNVs in just n cells in CN cluster
  # scDNA_CN_clusters_SNVs[scDNA_CN_clusters_SNVs <= filter_n_cells] <- 0
  # 
  # #Have to remove SNVs with no calls in all clusters as std dev is 0 and correlation distance doesn't work (from 31924 SNVs to just 9974 SNVs left)
  # scDNA_CN_clusters_SNVs <- scDNA_CN_clusters_SNVs[rowSums(scDNA_CN_clusters_SNVs) != 0,]
  # 
  # #Annotate CNA names
  # colnames(scDNA_CN_clusters_SNVs) <- names(scDNA_CN_cluster_ids)
  # 
  # saveRDS(scDNA_CN_clusters_SNVs, paste0(sub_folder, "/MPNST_DLP_scDNA_CN_clusters_SNVs.rds"))
}
cluster_plot_SNVs(SNV_set = "de_novo", G1000 = F, upper_limit = NULL, J = 36, plot = T)

#Create SNV matrix for de novo calls and apply filter for 1000G positions
if (T) {
  G1000PREFIX_AC = "../../../../../ref_files/Battenberg/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_loci_chrstring_" #modified maxime's file so doesn't have chr
  G1000_loci <- do.call(rbind, lapply(chrom_names, function(c) {
    print(c)
    read.table(file = paste0(G1000PREFIX_AC, "chr", c, ".txt"), header = F)
  }))
  G1000_loci$loci <- paste0(G1000_loci$V1, "_", G1000_loci$V2)
  
  mt_counts_SNV_raw <- de_novo_mt_counts_SNV_vec_cell[[1]][,1:2]
  mt_counts_SNV_raw$loci <- paste0(mt_counts_SNV_raw[,1], "_", mt_counts_SNV_raw[,2])
  mt_counts_SNV_filtered <- which(!mt_counts_SNV_raw$loci %in% G1000_loci$loci) #find positions which are in 1000G
  
  scDNA_CN_clusters_SNVs <- readRDS(paste0(sub_folder, "/MPNST_de_novo_scDNA_CN_clusters_SNVs.rds"))
  
  scDNA_CN_clusters_SNVs <- scDNA_CN_clusters_SNVs[!rownames(scDNA_CN_clusters_SNVs) %in% G1000_loci$loci,] #Keep only those no in 1000G loci
  # # 
  # # mt_counts_SNV_vec_cell_filtered <- pblapply(mt_counts_SNV_vec_cell, function(v) {
  # #   return(v[mt_counts_SNV_filtered,]) #remove those in 1000 genomes
  # # })
  # # 
  # # #For each CN cluster get number of SNVs
  # # scDNA_CN_clusters_SNVs <- do.call(cbind, lapply(1:length(scDNA_CN_cluster_ids), function(x) {
  # #   SNV_vec_cells <- mt_counts_SNV_vec_cell_filtered[scDNA_CN_cluster_ids[[x]]]
  # #   SNV_vec_cells_sum <- do.call(cbind, lapply(SNV_vec_cells, function(c) {
  # #     return(c$Present)
  # #   })) %>% rowSums()
  # #   return(SNV_vec_cells_sum)
  # # }))
  # # rownames(scDNA_CN_clusters_SNVs) <- paste0(mt_counts_SNV_vec_cell_filtered[[1]][,"chr"], "_", mt_counts_SNV_vec_cell_filtered[[1]][,"pos"])
  # # 
  # # #Keep just filtered SNVs
  # # mt_counts_SNV_mtx_filter <- readRDS("MPNST_mt_counts_SNV_mtx_filter.rds")
  # # scDNA_CN_clusters_SNVs <- scDNA_CN_clusters_SNVs[rownames(scDNA_CN_clusters_SNVs) %in% colnames(mt_counts_SNV_mtx_filter),]
  # # 
  # # #Remove bulk calls (only for DLP runs)
  # # CICC_output <- do.call(rbind, readRDS("../../Genotype_SNVs_CICC/MPNST_SNVs_clusters.rds"))
  # # scDNA_CN_clusters_SNVs <- scDNA_CN_clusters_SNVs[!rownames(scDNA_CN_clusters_SNVs) %in% paste0(CICC_output$chr, "_", CICC_output$pos),]
  # # 
  # # #Apply a filter to remove SNVs in just n cells
  # # scDNA_CN_clusters_SNVs[scDNA_CN_clusters_SNVs <= filter_n_cells] <- 0
  # # 
  # #Have to remove SNVs with no calls in all clusters as std dev is 0 and correlation distance doesn't work (6980 SNVs left)
  # scDNA_CN_clusters_SNVs <- scDNA_CN_clusters_SNVs[rowSums(scDNA_CN_clusters_SNVs) != 0,]
  # 
  # #Annotate CNA names
  # colnames(scDNA_CN_clusters_SNVs) <- names(scDNA_CN_cluster_ids)
  # 
  saveRDS(scDNA_CN_clusters_SNVs, paste0(sub_folder, "/MPNST_de_novo_scDNA_G1000filtered_CN_clusters_SNVs.rds"))
}
cluster_plot_SNVs(SNV_set = "de_novo", G1000 = T, upper_limit = NULL, J = 36, plot = T)

####################################################################################################################################
###Plot with bulk SNV calls
####################################################################################################################################
#Create SNV matrix for bulk calls
if (T) {
  #Get bulk calls for run
  #Load in SNVs of all cells and subset to pass filter cells
  bulk_mt_10X_counts_SNV_vec_cell <- readRDS("../../../../10X_DNA/results/Genotype_SNVs/MPNST_scDNA_mt_counts_SNV_vec_cell.rds")
  names(bulk_mt_10X_counts_SNV_vec_cell) <- unlist(barcodes)

  bulk_mt_DLP_counts_SNV_vec_cell <- readRDS("../../../../DLP_plus/results/Genotype_SNVs_CICC/MPNST_scDNA_mt_counts_SNV_vec_cell.rds")
  names(bulk_mt_DLP_counts_SNV_vec_cell) <- all_DLP_barcodes

  bulk_mt_counts_SNV_vec_cell <- lapply(c(bulk_mt_10X_counts_SNV_vec_cell, bulk_mt_DLP_counts_SNV_vec_cell), function(x) {
    x %>% arrange(chr, pos) #make sure in same order
  })
  
  #For each CN cluster get number of SNVs
  scDNA_CN_clusters_bulk_SNVs <- do.call(cbind, lapply(1:length(scDNA_CN_cluster_ids), function(x) {
    SNV_vec_cells <- bulk_mt_counts_SNV_vec_cell[scDNA_CN_cluster_ids[[x]]]
    SNV_vec_cells_sum <- do.call(cbind, lapply(SNV_vec_cells, function(c) {
      return(c$Present)
    })) %>% rowSums()
    return(SNV_vec_cells_sum)
  }))
  rownames(scDNA_CN_clusters_bulk_SNVs) <- paste0(bulk_mt_counts_SNV_vec_cell[[1]][,"chr"], "_", bulk_mt_counts_SNV_vec_cell[[1]][,"pos"])
  
  #Apply a filter to remove just one read
  scDNA_CN_clusters_bulk_SNVs[scDNA_CN_clusters_bulk_SNVs == 1] <- 0
  
  #Have to remove SNVs with no calls in all clusters as std dev is 0 and correlation distance doesn't work
  scDNA_CN_clusters_bulk_SNVs <- scDNA_CN_clusters_bulk_SNVs[rowSums(scDNA_CN_clusters_bulk_SNVs) != 0,]
  
  colnames(scDNA_CN_clusters_bulk_SNVs) <- names(scDNA_CN_cluster_ids)
  saveRDS(scDNA_CN_clusters_bulk_SNVs, paste0(sub_folder, "/MPNST_scDNA_CN_clusters_bulk_SNVs.rds"))
}
cluster_plot_SNVs(SNV_set = "all", G1000 = F, upper_limit = NULL, J = 36, plot = T)

#Introduce upper limit for clustering (worth doing!)
cluster_plot_SNVs(SNV_set = "all", G1000 = F, upper_limit = 2, J = 36, plot = T)

#Plot just bulk calls and colour by cluster
cluster_plot_SNVs(SNV_set = "bulk", G1000 = F, upper_limit = 2, J = NULL, plot = T)

cluster_plot_SNVs(SNV_set = "bulk", G1000 = F, upper_limit = 2, J = 30, plot = T)

####################################################################################################################################
###Plot with bulk SNV calls (try for G1000 filtered) and remove clonal/artefact clusters
####################################################################################################################################
#Run with G1000 filter
cluster_plot_SNVs(SNV_set = "all", G1000 = T, upper_limit = NULL, J = 36, plot = T)

#Introduce upper limit for clustering (worth it!)
cluster_plot_SNVs(SNV_set = "all", G1000 = T, upper_limit = 2, J = 36, plot = T)

#Remove clonal clusters
cluster_plot_SNVs(SNV_set = "all", G1000 = T, upper_limit = 2, ex_clonal = T, J = 28, plot = T)

#Reorder cluster
if (run == "22_kmeans_10X_DLP") {
  cluster_plot_SNVs(SNV_set = "all", G1000 = T, upper_limit = 2, ex_clonal = T, J = 32, ann_cluster = F, plot = T, 
                    reorder = c(7,8,9,#Clonal
                                29,30,#P
                                28,#R1extra
                                13,14,15,#R1/R5
                                32,31,12,#R1
                                17,18,16,19,10,6,5,11,#R5
                                23,#R2,3,4
                                22,24,20,25,21,#R2
                                26,#R3
                                4,3,1,2,27#R4
                    ))
}
