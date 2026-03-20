############################################################
#Plot totCN heatmap function
############################################################
sc_totCN_heatmap <- function(CN_mtx, hclust = F, km = NULL, row_split = NULL, probes, row_ann = NULL, max_CN = 10, column_title = "Single Cell Copy Number Heatmap", title = "Total copy number state", 
                             colour_scheme = c("#00008B", "#7B7BC3", "#FFFFFF", "#FFCCCC", "#FF8080", "#FF3333", "#E60000", "#B31000", "#8B0000", "#660000", "#330000")) {
  show_raw_dend = T
  if (!is.null(km)) {
    hclust = T
    show_raw_dend = F
  } #Cluster each kmeans cluster but don't show dendrogram
  CN_mtx[CN_mtx<0] <- 0 #set min CN value
  CN_mtx[CN_mtx>max_CN] <- max_CN #set max CN value
  names(colour_scheme) <- 0:max_CN #Colours
  CN_labels <- paste0(0:max_CN, " ")[which(names(colour_scheme) %in% sort(unique(as.numeric(CN_mtx))))] #Get labels in legend
  #Annotation
  ha_column = HeatmapAnnotation(odd = anno_empty(border = F),
                                even = anno_empty(border = F)) #Chromosome annotation at bottom
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(matrix = CN_mtx, cluster_rows = hclust, row_km = km, row_split = row_split, cluster_row_slices = FALSE, show_row_dend = show_raw_dend, 
                                               row_title_rot = 0, cluster_columns = F, show_row_names = F, show_column_names = F, row_title_gp = grid::gpar(fontsize = 24),
                                               heatmap_legend_param = list(labels = CN_labels, title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 20), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), border = "black", nrow=1),
                                               bottom_annotation = ha_column, left_annotation = row_ann, col = colour_scheme, column_title = column_title, name = title),
                       heatmap_legend_side="bottom", annotation_legend_side="bottom")
  #Add chr lines
  probes <- c(0, probes) #Adds posiiton zero for start of first chromosome
  decorate_heatmap_body(heatmap = title, {
    for (k in 1:(length(probes))) {
      grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1), gp=gpar(col="black", lty = 1, lwd = 1.5))
    }
  }) #Chr lines
  if (!is.null(row_split)) {
    for (i in 2:length(unique(row_split))) {
      decorate_heatmap_body(heatmap = title, row_slice = i, {
        for (k in 1:(length(probes))) {
          grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
        }
      }) #Chr lines
    }
  }
  if (!is.null(km)) {
    for (i in 2:K) {
      decorate_heatmap_body(heatmap = title, row_slice = i, {
        for (k in 1:(length(probes))) {
          grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
        }
      }) #Chr lines
    }
  }
  decorate_annotation(annotation = "odd", {
    for (k in 2:(length(probes))) {
      if (k%%2 == 0) {grid.text(names(probes)[k], x=((probes[k-1]+((probes[k]-probes[k-1])/2))/ncol(CN_mtx)) + 0.001, y=0.5, just = "left", gp=gpar(cex = 3))}
    }
  }) #Chr numbers
  decorate_annotation(annotation = "even", {
    for (k in 2:(length(probes))) {
      if (k%%2 == 1) {grid.text(names(probes)[k], x=((probes[k-1]+((probes[k]-probes[k-1])/2))/ncol(CN_mtx)) + 0.001, y=0.5, just = "left", gp=gpar(cex = 3))}
    }
  }) #Chr numbers
}

# test_mtx <- matrix(c(rep(0,10),rep(6,10),rep(2,10),rep(1,10),rep(5,10),rep(8,10)), ncol = 10)
# sc_totCN_heatmap(test_mtx, probes = c(1:10))

### Adapted function form Alex Stein
sc_totCN_heatmap2 <- function(CN_mtx, hclust = F, km = NULL, row_split = NULL, probes, row_ann = NULL, max_CN = 10, column_title = "Single Cell Copy Number Heatmap", title = "Total copy number state", 
                              colour_scheme = c("#00008B", "#7B7BC3", "#FFFFFF", "#FFCCCC", "#FF8080", "#FF3333", "#E60000", "#B31000", "#8B0000", "#660000", "#330000")) {
  show_raw_dend = T
  if (!is.null(km)) {
    hclust = T
    show_raw_dend = F
  } #Cluster each kmeans cluster but don't show dendrogram
  CN_mtx[CN_mtx<0] <- 0 #set min CN value
  CN_mtx[CN_mtx>max_CN] <- max_CN #set max CN value
  names(colour_scheme) <- 0:max_CN #Colours
  #CN_labels <- paste0(0:max_CN, " ")[which(names(colour_scheme) %in% sort(unique(as.numeric(CN_mtx))))] #Get labels in legend
  CN_at <- intersect(0:max_CN, sort(unique(as.numeric(CN_mtx))))
  CN_labels <- paste0(CN_at, " ")
  #Annotation
  ha_column = HeatmapAnnotation(odd = anno_empty(border = F),
                                even = anno_empty(border = F)) #Chromosome annotation at bottom
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(matrix = CN_mtx, cluster_rows = hclust, row_km = km, row_split = row_split, cluster_row_slices = FALSE, show_row_dend = show_raw_dend, 
                                               row_title_rot = 0, cluster_columns = F, show_row_names = F, show_column_names = F, row_title_gp = grid::gpar(fontsize = 24),
                                               heatmap_legend_param = list(at = CN_at, labels = CN_labels, title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 20), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), border = "black", nrow=1),
                                               use_raster = FALSE,
                                               bottom_annotation = ha_column, left_annotation = row_ann, col = colour_scheme, column_title = column_title, name = title),
                       heatmap_legend_side="bottom", annotation_legend_side="bottom")
  #Add chr lines
  probes <- c(0, probes) #Adds posiiton zero for start of first chromosome
  decorate_heatmap_body(heatmap = title, {
    for (k in 1:(length(probes))) {
      grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1), gp=gpar(col="black", lty = 1, lwd = 1.5))
    }
  }) #Chr lines
  if (!is.null(row_split)) {
    for (i in 2:length(unique(row_split))) {
      decorate_heatmap_body(heatmap = title, row_slice = i, {
        for (k in 1:(length(probes))) {
          grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
        }
      }) #Chr lines
    }
  }
  if (!is.null(km)) {
    for (i in 2:K) {
      decorate_heatmap_body(heatmap = title, row_slice = i, {
        for (k in 1:(length(probes))) {
          grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
        }
      }) #Chr lines
    }
  }
  decorate_annotation(annotation = "odd", {
    for (k in 2:(length(probes))) {
      if (k%%2 == 0) {grid.text(names(probes)[k], x=((probes[k-1]+((probes[k]-probes[k-1])/2))/ncol(CN_mtx)) + 0.001, y=0.5, just = "left", gp=gpar(cex = 3))}
    }
  }) #Chr numbers
  decorate_annotation(annotation = "even", {
    for (k in 2:(length(probes))) {
      if (k%%2 == 1) {grid.text(names(probes)[k], x=((probes[k-1]+((probes[k]-probes[k-1])/2))/ncol(CN_mtx)) + 0.001, y=0.5, just = "left", gp=gpar(cex = 3))}
    }
  }) #Chr numbers
}

############################################################
#Plot asCN heatmap function 
############################################################
sc_asCN_heatmap <- function(CN_mtx, hclust = F, km = NULL, row_split = NULL, probes, row_ann = NULL, column_title = "Single Cell Copy Number Heatmap", title = "Allele-specific copy number state", 
                            colour_scheme = c("royalblue3", "skyblue2", "grey80", "white", "gold1", "khaki1", "darkorange3", "darkorange1", "orange", "red4", "red", "orangered2", "purple4")) {
  show_raw_dend = T
  if (!is.null(km)) {
    hclust = T
    show_raw_dend = F
  } #Cluster each kmeans cluster but don't show dendrogram
  CN_mtx[CN_mtx<0] <- 0 #set min CN value
  CN_mtx[CN_mtx>60] <- 60 #set max CN value
  names(colour_scheme) <- c(0, 10, 20, 21, 30, 31, 40, 41, 42, 50, 51, 52, 60)#Colours
  CN_labels <- c("0+0 ", "1+0 ", "2+0 ", "1+1 ", "3+0 ", "2+1 ", "4+0 ", "3+1 ", "2+2 ", "5+0 ", "4+1 ", "3+2 ", ">6 ")[which(names(colour_scheme) %in% sort(unique(as.numeric(CN_mtx))))] #Get labels in legend
  ha_column = HeatmapAnnotation(odd = anno_empty(border = F),
                                even = anno_empty(border = F)) #Chromosome annotation at bottom
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(matrix = CN_mtx, cluster_rows = hclust, row_km = km, row_split = row_split, cluster_row_slices = FALSE, show_row_dend = show_raw_dend, 
                                               row_title_rot = 0, cluster_columns = F, show_row_names = F, show_column_names = F, row_title_gp = grid::gpar(fontsize = 24),
                                               heatmap_legend_param = list(labels = CN_labels, title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 20), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), border = "black", nrow=1),
                                               bottom_annotation = ha_column, left_annotation = row_ann, col = colour_scheme, column_title = column_title, name = title),
                       heatmap_legend_side="bottom", annotation_legend_side="bottom")
  probes <- c(0, probes) #Adds posiiton zero for start of first chromosome
  
  decorate_heatmap_body(heatmap = title, {
    for (k in 1:(length(probes))) {
      grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1), gp=gpar(col="black", lty = 1, lwd = 1.5))
    }
  }) #Chr lines
  if (!is.null(row_split)) {
    for (i in 2:length(unique(row_split))) {
      decorate_heatmap_body(heatmap = title, row_slice = i, {
        for (k in 1:(length(probes))) {
          grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
        }
      }) #Chr lines
    }
  }
  if (!is.null(km)) {
    for (i in 2:K) {
      decorate_heatmap_body(heatmap = title, row_slice = i, {
        for (k in 1:(length(probes))) {
          grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
        }
      }) #Chr lines
    }
  }
  decorate_annotation(annotation = "odd", {
    for (k in 2:(length(probes))) {
      if (k%%2 == 0) {grid.text(names(probes)[k], x=((probes[k-1]+((probes[k]-probes[k-1])/2))/ncol(CN_mtx)) + 0.001, y=0.5, just = "left", gp=gpar(cex = 3))}
    }
  }) #Chr numbers
  decorate_annotation(annotation = "even", {
    for (k in 2:(length(probes))) {
      if (k%%2 == 1) {grid.text(names(probes)[k], x=((probes[k-1]+((probes[k]-probes[k-1])/2))/ncol(CN_mtx)) + 0.001, y=0.5, just = "left", gp=gpar(cex = 3))}
    }
  }) #Chr numbers
}

# test_mtx <- matrix(c(rep(10,10),rep(60,10),rep(51,10),rep(20,10),rep(41,10),rep(40,10)), ncol = 10)
# sc_asCN_heatmap(test_mtx, probes = c(1:10))

############################################################
#Plot infercnv heatmap function
############################################################
infercnv_heatmap <- function(CN_mtx, hclust = F, km = NULL, row_split = NULL, probes, row_ann = NULL, column_title = "Inferred Copy Number Heatmap", title = "InferCNV signal intensity ", 
                             colour_scheme = c("darkblue", "white", "darkred")) {
  show_raw_dend = T
  if (!is.null(km)) {
    hclust = T
    show_raw_dend = F
  } #Cluster each kmeans cluster but don't show dendrogram
  colour_scheme <- circlize::colorRamp2(c(0.85, 1, 1.15), colour_scheme) #Colours
  #Annotation
  ha_column = HeatmapAnnotation(odd = anno_empty(border = F),
                                even = anno_empty(border = F)) #Chromosome annotation at bottom
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(matrix = CN_mtx, cluster_rows = hclust, row_km = km, row_split = row_split, cluster_row_slices = FALSE, show_row_dend = show_raw_dend, 
                                row_title_rot = 0, cluster_columns = F, show_row_names = F, show_column_names = F, row_title_gp = grid::gpar(fontsize = 24),
                                heatmap_legend_param = list(direction = "horizontal", at = seq(0.85,1.15,0.05), labels = c("", 0.9, 0.95, 1, 1.05, 1.1, ""), title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 20), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), legend_width = unit(10, "cm"), nrow=1),
                                bottom_annotation = ha_column, left_annotation = row_ann, col = colour_scheme, column_title = column_title, name = title),
                       heatmap_legend_side="bottom", annotation_legend_side="bottom")
  #Add chr lines
  probes <- c(0, probes) #Adds posiiton zero for start of first chromosome
  decorate_heatmap_body(heatmap = title, {
    for (k in 1:(length(probes))) {
      grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1), gp=gpar(col="black", lty = 1, lwd = 1.5))
    }
  }) #Chr lines
  if (!is.null(row_split)) {
    if (length(row_split) > 1) {
      for (i in 2:length(unique(row_split))) {
        decorate_heatmap_body(heatmap = title, row_slice = i, {
          for (k in 1:(length(probes))) {
            grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
          }
        }) #Chr lines
      }
    } else {
      for (i in 2:row_split) {
        decorate_heatmap_body(heatmap = title, row_slice = i, {
          for (k in 1:(length(probes))) {
            grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
          }
        }) #Chr lines
      }
    }
  }
  if (!is.null(km)) {
    for (i in 2:km) {
      decorate_heatmap_body(heatmap = title, row_slice = i, {
        for (k in 1:(length(probes))) {
          grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
        }
      }) #Chr lines
    }
  }
  decorate_annotation(annotation = "odd", {
    for (k in 2:(length(probes))) {
      if (k%%2 == 0) {grid.text(names(probes)[k], x=((probes[k-1]+((probes[k]-probes[k-1])/2))/ncol(CN_mtx)) + 0.001, y=0.5, just = "left", gp=gpar(cex = 3))}
    }
  }) #Chr numbers
  decorate_annotation(annotation = "even", {
    for (k in 2:(length(probes))) {
      if (k%%2 == 1) {grid.text(names(probes)[k], x=((probes[k-1]+((probes[k]-probes[k-1])/2))/ncol(CN_mtx)) + 0.001, y=0.5, just = "left", gp=gpar(cex = 3))}
    }
  }) #Chr numbers
}

# test_mtx <- matrix(c(rep(1,10),rep(0.9,10),rep(1.1,10),rep(1.15,10),rep(-0.85,10),rep(1.05,10)), ncol = 10)
# infercnv_heatmap(test_mtx, probes = c(1:10))

############################################################
#Plot infercnv heatmap function
############################################################
Visium_heatmap <- function(CN_mtx, hclust = F, km = NULL, row_split = NULL, probes, row_ann = NULL, column_title = "MPNST Visium InferCNV Relative Copy Number Heatmap", title = "InferCNV Intensity ", 
                             colour_scheme = c("darkblue", "royalblue1", "grey95", "firebrick3", "darkred")) {
  show_raw_dend = T
  if (!is.null(km)) {
    hclust = T
    show_raw_dend = F
  } #Cluster each kmeans cluster but don't show dendrogram
  names(colour_scheme) <- c(-2, -1, 0, 1, 2)#Colours
  CN_labels <- c("-2", "-1", "0", "1", "2")[which(names(colour_scheme) %in% sort(unique(as.numeric(CN_mtx))))] #Get labels in legend
  ha_column = HeatmapAnnotation(odd = anno_empty(border = F),
                                even = anno_empty(border = F)) #Chromosome annotation at bottom
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(matrix = CN_mtx, cluster_rows = hclust, row_km = km, row_split = row_split, cluster_row_slices = FALSE, show_row_dend = show_raw_dend, 
                                               row_title_rot = 0, cluster_columns = F, show_row_names = T, show_column_names = F, row_title_gp = grid::gpar(fontsize = 24),
                                               heatmap_legend_param = list(labels = CN_labels, title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 20), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), border = "black", nrow=1),
                                               bottom_annotation = ha_column, left_annotation = row_ann, col = colour_scheme, column_title = column_title, name = title),
                       heatmap_legend_side="bottom", annotation_legend_side="bottom")
  probes <- c(0, probes) #Adds posiiton zero for start of first chromosome
  
  decorate_heatmap_body(heatmap = title, {
    for (k in 1:(length(probes))) {
      grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1), gp=gpar(col="black", lty = 1, lwd = 1.5))
    }
  }) #Chr lines
  if (!is.null(row_split)) {
    for (i in 2:length(unique(row_split))) {
      decorate_heatmap_body(heatmap = title, row_slice = i, {
        for (k in 1:(length(probes))) {
          grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
        }
      }) #Chr lines
    }
  }
  if (!is.null(km)) {
    for (i in 2:K) {
      decorate_heatmap_body(heatmap = title, row_slice = i, {
        for (k in 1:(length(probes))) {
          grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
        }
      }) #Chr lines
    }
  }
  decorate_annotation(annotation = "odd", {
    for (k in 2:(length(probes))) {
      if (k%%2 == 0) {grid.text(names(probes)[k], x=((probes[k-1]+((probes[k]-probes[k-1])/2))/ncol(CN_mtx)) + 0.001, y=0.5, just = "left", gp=gpar(cex = 3))}
    }
  }) #Chr numbers
  decorate_annotation(annotation = "even", {
    for (k in 2:(length(probes))) {
      if (k%%2 == 1) {grid.text(names(probes)[k], x=((probes[k-1]+((probes[k]-probes[k-1])/2))/ncol(CN_mtx)) + 0.001, y=0.5, just = "left", gp=gpar(cex = 3))}
    }
  }) #Chr numbers
}

############################################################
#Do hclust on matrix object and save or load if already exists
############################################################
hclust_save_load <- function(object, sample = "", distance = "manhattan", method = "ward.D2") {
  file = paste0(sample, "_hclust_", distance, "_", method, ".rds")
  if (!file.exists(file)) {
    print(paste0("Clustering and saving hclust to ", file))
    hclust_obj <- fastcluster::hclust(dist(object, method = distance), method = method)
    saveRDS(hclust_obj, file)
    return(hclust_obj)
  } else {
    print(paste0("Loading hclust object from ", file))
    return(readRDS(file))
  }
}

############################################################
#Function for louvain clustering (takes in a cell x CN matrix)
############################################################
run_louvain <- function(sample = "sample", cell_matrix, save_intermediate = F, resolution = 1.5) {
  print("Generating KNN matrix")
  knn.info <- RANN::nn2(cell_matrix, k = 20)
  knn <- knn.info$nn.idx
  adj <- matrix(0, nrow(cell_matrix), nrow(cell_matrix)) #Adjacency matrix
  rownames(adj) <- colnames(adj) <- rownames(cell_matrix)
  for(i in seq_len(nrow(cell_matrix))) {
    adj[i,rownames(cell_matrix)[knn[i,]]] <- 1
  }
  #Convert as graph
  print("Converting to KNN graph")
  knn_graph <- igraph::graph.adjacency(adj, mode="undirected")
  knn_graph <- simplify(knn_graph) ## remove self loops
  if (save_intermediate) {
    saveRDS(knn_graph, paste0(sample, "_knn_graph.rds"))
  }
  #Do louvain clustering
  print("Running Louvain clustering")
  knn_louvain <- igraph::cluster_louvain(knn_graph, weights = NULL, resolution = resolution)
  saveRDS(knn_louvain, paste0(sample, "_knn_louvain.rds"))
  return(knn_louvain)
}

############################################################
#If doesn't exist save
############################################################
exist_load_save <- function(file, exp) {
  if (!file.exists(file)) {
    saveRDS(eval(exp), file)
    eval(exp)
  } else {
    readRDS(file)
  }
}
