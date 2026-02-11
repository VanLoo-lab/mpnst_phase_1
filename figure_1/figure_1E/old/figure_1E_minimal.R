####################################################################################################################################
### Part 6
####################################################################################################################################

gamma = 5
as_res <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/all_ASCAT.sc_ascn_mpcf_", gamma,".rds"))
asCN_probes <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/MPNST_all_asCN_probes_mpcf_", gamma,".rds"))
asCN_chr_probes <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/MPNST_all_asCN_chr_probes_mpcf_", gamma,".rds"))


#Generate asCN heatmap
scDNA_asCN_mtx <- exist_load_save(file = "MPNST_all_asfixed_asCN_mtx.rds",
                                  do.call(rbind, lapply(which(lengths(as_res[["allProfiles_AS"]])!=0), function(b){
                                    prof <- as_res$allProfiles_AS[[b]]$nprof.fixed[!is.na(as_res$allProfiles_AS[[b]]$nprof.fixed$total_copy_number),] #Using fixed
                                    prof$nA[is.na(prof$nA)] <- 0
                                    prof$nB[is.na(prof$nB)] <- prof$total_copy_number[which(is.na(prof$nB))]
                                    asCN_df <- data.frame(allele1 = prof$nA,
                                                          allele2 = prof$nB) %>% #flip alleles so A is larger
                                      mutate(alleleA = ifelse(allele1>=allele2, allele1, allele2),
                                             alleleB = ifelse(allele1>=allele2, allele2, allele1)) %>%
                                      mutate(alleleA = ifelse(alleleA<0, 0, alleleA), #code to fix CN of -1
                                             alleleB = ifelse(alleleB<0, 0, alleleB)) %>%
                                      mutate(alleleT = alleleA + alleleB) %>% #get total CN
                                      mutate(plotT = alleleT*10+alleleB) %>% select(-c(1:2))#Convert to matrix format
                                    rep(na.omit(as.numeric(asCN_df$plotT)), asCN_probes$n.probes) #NAs are omitted...
                                  })) %>% set_rownames(names(as_res[["allProfiles_AS"]])))


MPNST_asCN_hclust_ward <- readRDS(paste0(INPUTDIR, "10x_DLP/ASCAT.sc/MPNST_all_asCN_hclust_ward.rds"))

#Plot CN profiles
ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_asCN_mtx)),
                                       Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_asCN_mtx))),
                       col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
                                  Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")),
                       annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 22), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1),
                                                      Tech = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 22), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
png(filename = paste0("MPNST_all_asfixed_asCN_heatmaptest.png"), width = 4000, height = 4000, res = 200)
print(sc_asCN_heatmap(CN_mtx = scDNA_asCN_mtx, hclust = MPNST_asCN_hclust_ward, row_ann = ha_row, probes = asCN_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X"))))
dev.off()

as_CN_colour <- c("royalblue3", "skyblue2", "grey80", "white", "gold1", "khaki1", "darkorange3", "darkorange1", "orange", "red4", "red", "orangered2", "purple4")
as_CN_breaks = c(0,10,20,21,30,31,40,41,42,50,51,52,60,1000)-0.1
png(filename = paste0("MPNST_all_asCN_colourkey.png"), width = 2000, height = 500, res = 200)
barplot <- barplot(rep(1,length(as_CN_colour)), col = as_CN_colour)
text(barplot, .2, c("0+0", "1+0", "2+0", "1+1", "3+0", "2+1", "4+0", "3+1", "2+2", "5+0", "4+1", "3+2", ">6"), 0, cex =1, pos=3)
dev.off()

#Plot CN profiles of clusters

k_means_clusters <- readRDS(paste0("MPNST_all_k_means_K",K,"_clusters_hclust.rds")) #using hclust on each kmeans

ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_asCN_mtx[names(named_clusters),])),
                                       Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_asCN_mtx[names(named_clusters),]))),
                       col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
                                  Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")),
                       annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 22), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1),
                                                      Tech = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 22), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
png(filename = paste0("MPNST_all_asfixed_asCN_K", K, "_named.png"), width = 4000, height = 4000, res = 200)
named_clusters <- unlist(lapply(1:length(scDNA_CN_cluster_ids), function(i) return(rep(names(scDNA_CN_cluster_ids)[i], length(scDNA_CN_cluster_ids[[i]]))))) %>% set_names(unlist(scDNA_CN_cluster_ids))
print(sc_asCN_heatmap(CN_mtx = scDNA_asCN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, row_ann = ha_row, column_title = NA, probes = asCN_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X"))))
dev.off()