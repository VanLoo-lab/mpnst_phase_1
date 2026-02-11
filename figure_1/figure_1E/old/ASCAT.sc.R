### title: Plot the copy number profiles of single cells. 
# The following code is used to generate Figure 1E (and other figures)

### To run the code:
# (1) Make sure to load SessionInfo() for correct package version 
# (2) Set both INPUTDIR and OUTPUTDIR to correct path
# (3) Part 4 contains essential code to recreate Figure 1B.
#     Parts can be ignored by setting False in the "if" commands, i.e. "if(F)".

# Computational cost is not too high such that the code should be feasible to 
#   run on local machines.

####################################################################################################################################
### Part 0: Preperation 
####################################################################################################################################

# In case plots look weird, the following setting may fix it
#options(bitmapType='cairo')

### Load libraries
library(tidyverse)
library(ASCAT.sc)
library(parallel)
library(pbmcapply)
library(fastcluster)
library(ComplexHeatmap)
library(grid)
library(gridExtra)
library(magrittr)


INPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_1/data/"
OUTPUTDIR <- "~/Documents/GitHub/MPNST-Zenodo/figure_1/results/figure_1e/"
CODEDIR <- "~/Documents/GitHub/MPNST-Phase-1/figure_1/figure_1E/"

### Load metadata
source( paste0(CODEDIR, "metadata/MPNST_common_functions.R") )
source( paste0(CODEDIR, "metadata/10X_DLP_metadata_minimal.R") )

### Prepare data input
#FASTA <- "/camp/project/proj-vanloo/reference_files/human/references/alignment/hs38DH/hs38DH.fa"

scDNA_CN_colour <- c("#00008B", "#7B7BC3", "#FFFFFF", "#FFCCCC", "#FF8080", "#FF3333", "#E60000", "#B31000", "#8B0000", "#660000", "#330000")
names(scDNA_CN_colour) <- 0:10
scDNA_CN_breaks = seq(-0.5,10.5, by = 1)

####################################################################################################################################
### Part 1: Get probes for plotting non-mpcf run
####################################################################################################################################

if (T) {
  data("lSe_filtered_30000.hg38", package = "ASCAT.sc")
  binsize = 500000
  nlSe <- treatlSe(lSe.hg38.filtered[1:23], window = ceiling(binsize/30000))
  
  chr_probe_positions <- do.call(rbind, lapply(1:length(nlSe), function(c){
    return(data.frame(chr = names(nlSe)[c], probe.pos = nlSe[[c]][["starts"]]))
  }))
  chr_probes <- rbind(data.frame(chr = 0, probes = 0), chr_probe_positions %>% 
                        mutate(chr = as.numeric(gsub("X", 23, gsub("chr", "", chr)))) %>% 
                        group_by(chr) %>% 
                        summarise(probes = n())) %>% mutate(cum.probes = cumsum(probes))
}

# ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx))),
#                        col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D")), show_annotation_name = F)

output.dir = OUTPUTDIR
setwd(output.dir)

####################################################################################################################################
### Part 2: 25/05/2023
####################################################################################################################################

#First run (no multipcf) to remove bad quality cells from each region for 10X
#For DLP already run
if (F) {
  for (s in samples[5]) {
    print(s)
    system(paste0("mkdir -p ", output.dir, s, "_10X"))
    BAMS=paste0(input.dir[which(names(input.dir) == s)], "possorted_bam.bam")
    
    # barcodes_orig <- gsub("-1", "", as.character(read.csv(barcode.dir[which(names(input.dir) == s)])[,"barcode"]))
    barcodes_orig <- as.character(read.csv(barcode.dir[which(names(input.dir) == s)])[,"barcode"])
    #Run on individual regions
    if (!file.exists(paste0(s, "_10X_ASCAT.sc.rds"))) {
      res=run_sc_sequencing(tumour_bams=BAMS,
                            res = NULL,
                            allchr=paste0("chr", c(1:22,"X")),
                            sex=rep('male',length(BAMS)),
                            chrstring_bam='chr',
                            purs = c(0.5,1),
                            ploidies = seq(1.7,5,0.01),
                            maxtumourpsi=5,
                            binsize=500000,
                            segmentation_alpha = 0.01,
                            predict_refit = TRUE,
                            print_results = TRUE,
                            build="hg38",
                            MC.CORES=30,
                            barcodes_10x = barcodes_orig, 
                            outdir = paste0("./", s, "_10X/"),  #Subfolder to save files
                            probs_filters = 0.1, 
                            path_to_phases = NULL, 
                            list_ac_counts_paths = NULL,
                            sc_filters = TRUE, 
                            projectname = paste0("MPNST_10X_",s),
                            steps = NULL,
                            smooth_sc = FALSE, 
                            multipcf = FALSE)
      saveRDS(res,file=paste0(s, "_10X_ASCAT.sc.rds"))
    } else {
      res <- readRDS(paste0(s, "_10X_ASCAT.sc.rds"))
      if (F) {
        getFilters_hxy <- function (res, thresholdNrec = NULL, probs = 0.1, outdir = "./", projectname = "") {
          filters <- NULL
          try({
            allT <- res$allTracks.processed
            allS <- res$allSolutions
            getloess <- function(qu, nr) {
              nms <- paste0("n", 1:length(qu))
              names(qu) <- nms
              quo <- qu[order(nr, decreasing = F)]
              fitted <- stats::runmed(quo, k = 31, endrule = "keep")
              names(fitted) <- names(quo)
              list(fitted = fitted[nms], residuals = quo[nms] - 
                     fitted[nms])
            }
            getQuality.SD <- function(allT) {
              sapply(allT, function(x) {
                median(abs(diff(unlist(lapply(x$lCTS, function(y) y$smoothed)))))
              })
            }
            nrecords <- sapply(allT, function(x) sum(unlist(lapply(x$lCTS, function(y) y$records))))
            thresholdNrec <- ifelse(is.null(thresholdNrec), quantile(nrecords, probs = probs), thresholdNrec) #Use 10% quantile unless specified
            ambiguous <- sapply(allS, function(x) x$ambiguous)
            doublet <- sapply(allS, function(x) if (!is.null(x$bestfit)) 
              !x$bestfit$ambiguous
              else F)
            qualities <- getQuality.SD(allT)
            thresholdQual <- quantile(qualities[nrecords > thresholdNrec], probs = 1 - probs) #Quantile only on those that pass thresholdNrec
            keep <- qualities <= thresholdQual & nrecords >= thresholdNrec & 
              !ambiguous & !doublet
            keep2 <- !(qualities < thresholdQual & nrecords < thresholdNrec)
            keep2 <- keep2 & !ambiguous
            ll <- getloess(qualities[keep2], log2(nrecords)[keep2])
            ord <- order(log2(nrecords)[keep2], decreasing = F)
            filters <- (1:length(nrecords)) %in% (which(keep2)[ll$residuals <=  0.02]) & keep
            #Generate PDF
            pdf(paste0(outdir, "/", projectname, "_sc_filters.pdf"))
            plot(qualities, log2(nrecords), xlab = "Noise logr", ylab = "Total number of reads", 
                 pch = ifelse(doublet, 15, 19), cex = ifelse(doublet, 0.3, 0.1), #Doublets larger
                 col = ifelse(ambiguous, rgb(1, 0, 0.5, 0.5), rgb(0, 0, 0, 0.5))) #Ambiguous in pink
            points(ll$fitted[ord], log2(nrecords)[keep2][ord], type = "l", col = rgb(1, 0, 0, 0.5), lwd = 1.5) #loess line
            abline(v = thresholdQual, h = log2(thresholdNrec))
            points(qualities[keep2], log2(nrecords)[keep2], col = ifelse(abs(ll$residuals) > 0.02, rgb(1, 0, 0, 1), rgb(0, 0, 0, 0)), cex = 0.15) #Noisier cells
            #Density plot
            try({
              plot(density(log2(nrecords[!doublet])))
              polygon(density(log2(nrecords[doublet])), border = "red")
            }, silent = T)
            #Filter plot
            plot(qualities, log2(nrecords), xlab = "Noise logr", ylab = "Total number of reads", 
                 pch = ifelse(filters, 15, 19), cex = ifelse(filters, 0.3, 0.1), 
                 col = ifelse(filters, rgb(1, 0, 0.5, 0.5), rgb(0, 0, 0, 0.5))) #cells passing are pink and square
            dev.off()
          })
          res$filters <- filters
          res
        }
        res <- getFilters_hxy(res, thresholdNrec = 3e5, probs = 0.1, outdir = paste0("./", s, "/"), 
                              projectname = paste0("MPNST_DLP_",s))
        
        
        # res=run_sc_sequencing(tumour_bams=names(res[["allTracks"]]),
        #                       res = res,
        #                       allchr=paste0("chr", c(1:22,"X")),
        #                       sex=rep('male',length(res[["allTracks"]])),
        #                       chrstring_bam='chr',
        #                       purs = c(0.5,1),
        #                       ploidies = seq(1.7,5,0.01),
        #                       maxtumourpsi=5,
        #                       binsize=500000,
        #                       segmentation_alpha = 0.01,
        #                       predict_refit = FALSE,
        #                       print_results = FALSE, #No need to reprint
        #                       build="hg38",
        #                       MC.CORES=60,
        #                       barcodes_10x = NULL, 
        #                       outdir = paste0("./", s, "/"),  #Subfolder to save files
        #                       probs_filters = 0.2, #R1 0.1, R2 0.1, R3 0.12, R4 0.1, R5 0.2, P 0.25
        #                       path_to_phases = NULL, 
        #                       list_ac_counts_paths = NULL,
        #                       sc_filters = TRUE, 
        #                       projectname = paste0("MPNST_DLP_",s),
        #                       steps = NULL,
        #                       smooth_sc = FALSE, 
        #                       multipcf = FALSE) #Rerun with different filters
        saveRDS(res,file=paste0(s, "_filtered_ASCAT.sc.rds"))
      } #didn't run this filter for 10X
    }
    
    #Only needed if first run didn't have refit
    if (F) {
      res_mod <- res
      nms <- names(res_mod$allTracks)
      res_mod$allTracks <- lapply(nms, function(x) {
        res_mod$allTracks[[x]] <- res_mod$allTracks[[x]]$lCTS.tumour
        return(res_mod$allTracks[[x]])
      })
      names(res_mod$allTracks) <- nms
      
      res=run_sc_sequencing(tumour_bams=BAMS,
                            res = res_mod,
                            allchr=paste0("chr", c(1:22,"X")),
                            sex=rep('male',length(BAMS)),
                            chrstring_bam='chr',
                            purs = c(0.5,1),
                            ploidies = seq(1.7,5,0.01),
                            maxtumourpsi=5,
                            binsize=500000,
                            segmentation_alpha = 0.01,
                            predict_refit = TRUE,
                            print_results = TRUE,
                            build="hg38",
                            MC.CORES=60,
                            barcodes_10x = barcodes_orig, 
                            outdir = paste0("./", s, "_10X/"),  #Subfolder to save files
                            probs_filters = 0.1, 
                            path_to_phases = NULL, 
                            list_ac_counts_paths = NULL,
                            sc_filters = TRUE, 
                            projectname = paste0("MPNST_10X_",s),
                            steps = NULL,
                            smooth_sc = FALSE, 
                            multipcf = FALSE)
      saveRDS(res_mod,file=paste0(s, "_10X_ASCAT.sc.rds"))
    }
    
    # #Generate CN mtx
    # allProfiles <- res$allProfiles #NOTE using refitted here
    # if (!file.exists(paste0("MPNST_10X_",s,"_CN_mtx.rds"))) {
    #   scDNA_CN_mtx <- do.call(rbind, pbmclapply(1:length(allProfiles), function(b){
    #     unlist(lapply(1:nrow(chr_probe_positions), function(r) {
    #       probe_CN <- allProfiles[[b]] %>% as_tibble %>% filter(chromosome == chr_probe_positions[r,1], as.numeric(start) < chr_probe_positions[r,2], as.numeric(end) > chr_probe_positions[r,2]) %>%
    #         select(total_copy_number) %>% pull()
    #       return(ifelse(is.null(probe_CN), NA, as.numeric(probe_CN)))
    #     }))
    #   }, mc.cores = 30))
    #   saveRDS(scDNA_CN_mtx, paste0("MPNST_10X_",s,"_CN_mtx.rds"))
    # } else {
    #   scDNA_CN_mtx <- readRDS(paste0("MPNST_10X_",s,"_CN_mtx.rds"))
    # }
    
    #Generate CN mtx
    allProfiles <- res$allProfiles.refitted.auto #NOTE using refitted here
    if (!file.exists(paste0("MPNST_10X_refit_",s,"_CN_mtx.rds"))) {
      scDNA_CN_mtx <- do.call(rbind, pbmclapply((1:length(allProfiles))[unlist(lapply(allProfiles, function(x) !is.null(dim(x))))], function(b){
        unlist(lapply(1:nrow(chr_probe_positions), function(r) {
          probe_CN <- allProfiles[[b]] %>% as_tibble %>% filter(chromosome == chr_probe_positions[r,1], as.numeric(start) < chr_probe_positions[r,2], as.numeric(end) > chr_probe_positions[r,2]) %>%
            select(total_copy_number) %>% pull()
          return(ifelse(is.null(probe_CN), NA, as.numeric(probe_CN)))
        }))
      }, mc.cores = 30))
      saveRDS(scDNA_CN_mtx, paste0("MPNST_10X_refit_",s,"_CN_mtx.rds"))
    } else {
      scDNA_CN_mtx <- readRDS(paste0("MPNST_10X_refit_",s,"_CN_mtx.rds"))
    }

    #Generate raw heatmap
    MPNST_CN_hclust_man_ward <- hclust_save_load(scDNA_CN_mtx, sample = paste0("MPNST_10X_refit_CN_",s,"_raw"), distance = "manhattan", method = "ward.D2" )
    png(filename = paste0("MPNST_10X_refit_",s,"_scDNA_CN_raw.png"), width = 4000, height = 4000, res = 200)
    print(sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx, hclust = MPNST_CN_hclust_man_ward, probes = chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), title = "MPNST 10X \nTotal CN "))
    dev.off()

    #Generate pass heatmap
    MPNST_CN_hclust_man_ward <- hclust_save_load(scDNA_CN_mtx[res$filters[unlist(lapply(allProfiles, function(x) !is.null(dim(x))))],], sample = paste0("MPNST_10X_refit_CN_",s,"_pass"), distance = "manhattan", method = "ward.D2" )
    png(filename = paste0("MPNST_10X_refit_",s,"_scDNA_CN_pass.png"), width = 4000, height = 4000, res = 200)
    print(sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[res$filters[unlist(lapply(allProfiles, function(x) !is.null(dim(x))))],], hclust = MPNST_CN_hclust_man_ward, probes = chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), title = "MPNST 10X \nTotal CN "))
    dev.off()
  }
}

####################################################################################################################################
### Part 3: Run multi-PCF on 10X cells (didn't complete, mpcf crashing)
####################################################################################################################################
if (F) {
  gamma = 5
  system(paste0("mkdir -p ", output.dir, "10X_mpcf_", gamma))
  res_list <- c(lapply(samples[1:6], function(s) readRDS(paste0(s, "_10X_ASCAT.sc.rds"))))
  # res_list <- c(lapply(samples[1:2], function(s) readRDS(paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/DLP_plus/results/ASCAT.sc/", s, "_ASCAT.sc.rds"))),
  #               lapply(samples[3:6], function(s) readRDS(paste0(s, "_ASCAT.sc.rds"))))
  names(res_list) <- NULL
  
  res_comb <- list(allTracks = do.call(c, lapply(1:length(samples), function(s) {
    res <- res_list[[s]]$allTracks
    names(res) <- paste0(samples[s], "_", names(res))
    return(res[res_list[[s]]$filters])
  })))
  timetoread_tumours = 1
  
  # BAMS <- unlist(lapply(1:length(samples), function(s) {
  #   paste0(DLP.singlebam.dir[s], "/", names(res_list[[s]][["filters"]])[res_list[[s]][["filters"]]])
  # }))
  
  all_res=run_sc_sequencing(tumour_bams=names(res_comb[["allTracks"]]),
                            res = res_comb,
                            allchr=paste0("chr", c(1:22,"X")),
                            sex=rep('male',length(res_comb[["allTracks"]])),
                            chrstring_bam='chr',
                            purs = c(0.5,1),
                            ploidies = seq(1.7,5,0.01),
                            maxtumourpsi=5,
                            binsize=500000,
                            segmentation_alpha = (1/gamma),
                            predict_refit = TRUE,
                            print_results = TRUE,
                            build="hg38",
                            MC.CORES=40,
                            barcodes_10x = NULL, 
                            outdir = paste0("./10X_mpcf_",gamma,"/"),  #Subfolder to save files
                            probs_filters = 0.1, 
                            path_to_phases = NULL, 
                            list_ac_counts_paths = NULL,
                            sc_filters = TRUE, 
                            projectname = paste0("MPNST_10X_All"),
                            steps = NULL,
                            smooth_sc = FALSE, 
                            multipcf = TRUE)
  saveRDS(all_res,file=paste0("all_10X_ASCAT.sc_mpcf_", gamma,".rds"))
  
  if (!file.exists(paste0("MPNST_10X_probes_mpcf_", gamma,".rds"))) {
    tenX_probes <- do.call(rbind, lapply(all_res$allTracks.processed[[1]][["lSegs"]], function(s) return(s[["output"]]))) %>% 
      select(chrom, num.mark, loc.start, loc.end) %>% dplyr::rename(n.probes = num.mark, startpos = loc.start, endpos = loc.end) %>%
      filter(n.probes>1) #Remove segments of just 1/2 probes(these are NAs in ascat.sc outputs)
    
    tenX_chr_probes <- rbind(data.frame(chrom = 0, total.probes = 0), tenX_probes[,1:2] %>% group_by(chrom) %>% summarise(total.probes = sum(n.probes))) %>% 
      mutate(cum.probes = cumsum(total.probes))
    saveRDS(tenX_probes, paste0("MPNST_10X_probes_mpcf_", gamma,".rds"))
    saveRDS(tenX_chr_probes, paste0("MPNST_10X_chr_probes_mpcf_", gamma,".rds"))
  } else {
    tenX_probes <- readRDS(paste0("MPNST_10X_probes_mpcf_", gamma,".rds"))
    tenX_chr_probes <- readRDS(paste0("MPNST_10X_chr_probes_mpcf_", gamma,".rds"))
  }
  
  if (!file.exists(paste0("MPNST_10X_CN_mtx_mpcf_", gamma,".rds"))) {
    scDNA_CN_mtx <- do.call(rbind, lapply(1:length(all_res$allProfiles), function(b){
      rep(na.omit(as.numeric(all_res$allProfiles[[b]][,"total_copy_number"])), tenX_probes$n.probes)
    }))
    # names(all_barcodes) <- unlist(lapply(barcodes, function(s) return(as.character(s$sample_id)))) %>% unname()
    # rownames(scDNA_CN_mtx) <- all_barcodes[gsub(".bam", "", names(all_res$allProfiles))]
    saveRDS(scDNA_CN_mtx, paste0("MPNST_10X_CN_mtx_mpcf_", gamma,".rds"))
  } else {
    scDNA_CN_mtx <- readRDS(paste0("MPNST_10X_CN_mtx_mpcf_", gamma,".rds"))
  }
  
  #Generate raw heatmap
  MPNST_CN_hclust_man_ward <- hclust_save_load(scDNA_CN_mtx, sample = paste0("MPNST_10X_CN_all_mpcf_", gamma,"_raw"), distance = "manhattan", method = "ward.D2" )
  ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx))),
                         col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D")), show_annotation_name = F)
  png(filename = paste0("MPNST_10X_all_scDNA_CN_mpcf_", gamma,"_raw.png"), width = 4000, height = 4000, res = 200)
  print(sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx, hclust = MPNST_CN_hclust_man_ward, row_ann = ha_row, probes = tenX_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), title = "MPNST 10X \nTotal CN "))
  dev.off()
  
  #Generate pass heatmap
  MPNST_CN_hclust_man_ward <- hclust_save_load(scDNA_CN_mtx[all_res$filters,], sample = paste0("MPNST_10X_CN_all_mpcf_", gamma,"_pass"), distance = "manhattan", method = "ward.D2" )
  ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[all_res$filters,]))),
                         col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D")), show_annotation_name = F)
  png(filename = paste0("MPNST_10X_all_scDNA_CN_mpcf_", gamma,"_pass.png"), width = 4000, height = 4000, res = 200)
  print(sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[all_res$filters,], hclust = MPNST_CN_hclust_man_ward, row_ann = ha_row, probes = tenX_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), title = "MPNST 10X \nTotal CN "))
  dev.off()
  
  #Probably don't need second filter so used raw
  MPNST_CN_hclust_man_ward <- hclust_save_load(scDNA_CN_mtx, sample = paste0("MPNST_DLP_CN_all_mpcf_", gamma,"_raw"), distance = "manhattan", method = "ward.D2" )
  nonwgd_cells <- names(cutree(MPNST_CN_hclust_man_ward, k = 3)[cutree(MPNST_CN_hclust_man_ward, k = 3) == 1]) #Cluster 2 and 3 is WGD cells
  K = 22
  
  png(filename = paste0("MPNST_DLP_all_scDNA_CN_mpcf_", gamma,"_K", K, ".png"), width = 4000, height = 4000, res = 200)
  #Get clusters
  if (!file.exists(paste0("MPNST_DLP_k_means_K",K,"_clusters.rds"))) {
    set.seed(123)
    hm <- ComplexHeatmap::Heatmap(matrix = scDNA_CN_mtx[nonwgd_cells,], cluster_rows = T, cluster_row_slices = FALSE, row_km = K, row_km_repeats = 100)
    clusters <- row_order(hm)
    k_means_clusters <- lapply(1:length(clusters), function(j) {
      return(rownames(scDNA_CN_mtx[nonwgd_cells,])[clusters[[j]]])
    })
    saveRDS(k_means_clusters, paste0("MPNST_DLP_k_means_K",K,"_clusters.rds"))
  } else {
    k_means_clusters <- readRDS(paste0("MPNST_DLP_k_means_K",K,"_clusters.rds"))
  }
  named_clusters <- unlist(lapply(1:length(k_means_clusters), function(i) return(rep(i, length(k_means_clusters[[i]]))))) %>% set_names(unlist(k_means_clusters))
  ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[names(named_clusters),]))),
                         col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D")), show_annotation_name = F)
  sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, probes = DLP_chr_probes$cum.probes[-1] %>% set_names(nm = 1:23), row_ann = ha_row, title = "MPNST DLP \nTotal CN ")
  dev.off()
}

####################################################################################################################################
### Part 4: Get allele specifc calls for 10X (didn't run)
####################################################################################################################################

if (F) {
  gamma = 5
  system(paste0("mkdir -p ", output.dir, "ascn_mpcf_", gamma))
  
  path_to_phases <- list(paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/bulk/results/BB_spiked/All_with_Ext/",
                                "normal_multisample_haplotypes_chr",c(1:22,"23"),".vcf"))
  ac_path <- "/camp/project/proj-vanloo/analyses/hyan/mpnst/DLP_plus/results/ASCAT.sc/asCN/Allele_Counts/"
  
  list_ac_counts_paths <- lapply(unlist(lapply(samples, function(s) as.character(barcodes[[s]]$sample_id))), function(b){
    paste0(ac_path, b, "_alleleFrequencies_chr", c(1:22,"23"), ".txt")
  }) %>% set_names(paste0(unlist(lapply(samples, function(s) as.character(barcodes[[s]]$sample_id))), ".bam"))
  
  #Turn allele freq into per cell
  if (F) {
    for (s in samples[5:6]) {
      for (c in 1:23) {
        print(c)
        allele_counts <- read.delim(paste0("../../Genotype_SNPs/Allele_Freq/",s,"/", s,"_alleleFrequencies_chr",c,".txt"), stringsAsFactors = F)
        colnames(allele_counts)[1] <- "#CHR"
        pbmclapply(1:nrow(barcodes[[s]]), function(b) {
          write.table(allele_counts %>% filter(Barcode == paste0(barcodes[[s]]$BARCODE[b], "-0")), file = paste0(ac_path, barcodes[[s]]$sample_id[b], "_alleleFrequencies_chr", c, ".txt"), quote = F, row.names = F, sep = "\t")
        }, mc.cores = 10)
      }
    }
  }
  
  if (!file.exists(paste0("all_ASCAT.sc_ascn_mpcf_", gamma,".rds"))) {
    all_res <- readRDS(file=paste0("all_ASCAT.sc_mpcf_", gamma,".rds"))
    
    as_res <- getAS_CNA(all_res, path_to_phases = path_to_phases, 
                        list_ac_counts_paths = list_ac_counts_paths[names(all_res[["allTracks"]])], purs = c(0.5,1), 
                        ploidies = seq(1.7,5,0.01), outdir = aste0("./ascn_mpcf_",gamma,"/"), projectname = paste0("MPNST_DLP_All"), 
                        steps = NULL, mc.cores = 30)
    names(as_res$allProfiles_AS) <- names(as_res$allProfiles)
    
    # as_res <- run_sc_sequencing(tumour_bams=names(all_res[["allTracks"]]),
    #                             res = all_res,
    #                             allchr=paste0("chr", c(1:22,"X")),
    #                             sex=rep('male',length(all_res[["allTracks"]])),
    #                             chrstring_bam='chr',
    #                             purs = c(0.5,1),
    #                             ploidies = seq(1.7,5,0.01),
    #                             maxtumourpsi=5,
    #                             binsize=500000,
    #                             segmentation_alpha = (1/gamma),
    #                             predict_refit = FALSE,
    #                             print_results = FALSE,
    #                             build="hg38",
    #                             MC.CORES=40,
    #                             barcodes_10x = NULL, 
    #                             outdir = paste0("./ascn_mpcf_",gamma,"/"),  #Subfolder to save files
    #                             probs_filters = 0.1, 
    #                             path_to_phases = path_to_phases, 
    #                             list_ac_counts_paths = list_ac_counts_paths[names(all_res[["allTracks"]])],
    #                             sc_filters = TRUE, 
    #                             projectname = paste0("MPNST_DLP_All"),
    #                             steps = NULL,
    #                             smooth_sc = FALSE, 
    #                             multipcf = TRUE)
    
    saveRDS(as_res,file=paste0("all_ASCAT.sc_ascn_mpcf_", gamma,".rds"))
    as_res <- getAS_CNA_smoothed(as_res, mc.cores = 30)
    
  }
  else {
    as_res <- readRDS(paste0("all_ASCAT.sc_ascn_mpcf_", gamma,".rds"))
  }
  
  if (!file.exists(paste0("MPNST_DLP_asCN_probes_mpcf_", gamma,".rds"))) {
    asCN_probes <- do.call(rbind, lapply(as_res$allTracks.processed[[1]][["lSegs"]], function(s) return(s[["output"]]))) %>% 
      select(chrom, num.mark, loc.start, loc.end) %>% dplyr::rename(n.probes = num.mark, startpos = loc.start, endpos = loc.end) %>%
      filter(n.probes>1) #Remove segments of just 1/2 probes(these are NAs in ascat.sc outputs)
    
    asCN_chr_probes <- rbind(data.frame(chrom = 0, total.probes = 0), asCN_probes[,1:2] %>% group_by(chrom) %>% summarise(total.probes = sum(n.probes))) %>% 
      mutate(cum.probes = cumsum(total.probes))
    saveRDS(asCN_probes, paste0("MPNST_DLP_asCN_probes_mpcf_", gamma,".rds"))
    saveRDS(asCN_chr_probes, paste0("MPNST_DLP_asCN_chr_probes_mpcf_", gamma,".rds"))
  } else {
    asCN_probes <- readRDS(paste0("MPNST_DLP_asCN_probes_mpcf_", gamma,".rds"))
    asCN_chr_probes <- readRDS(paste0("MPNST_DLP_asCN_chr_probes_mpcf_", gamma,".rds"))
  }
  
  #Generate asCN heatmap
  scDNA_asCN_mtx <- exist_load_save(file = "MPNST_DLP_asfree_asCN_mtx.rds",
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
  
  #Function to calculate allele specific distance
  ascn_dist <- function(cna1, cna2) {
    # ind1 <- is.na(cna1[,"nA"]) | is.na(cna2[,"nA"])
    # ind2 <- !ind1 
    ind1 <- (is.na(cna1[,"nA"]) | is.na(cna2[,"nB"])) #Edited as one segment has NA
    ind2 <- !(is.na(cna1[,"nA"]) | is.na(cna2[,"nB"]))
    sizes <- cna1[,"endpos"]-cna1[,"startpos"]
    sizes <- sizes/500000/sum(sizes/500000)
    dists <- sum(pmin(2,abs(cna1[ind2,"nA"]-cna2[ind2,"nA"])+abs(cna1[ind2,"nB"]-cna2[ind2,"nB"]))*sizes[ind2])
    # dists <- dists+sum(pmin(2,abs(cna1[ind1,"fitted"]-cna2[ind1,"fitted"])*sizes[ind1]))
    dists
  }
  
  alldists <- matrix(NA,length(as_res[["allProfiles_AS"]]),length(as_res[["allProfiles_AS"]]))
  for(i in 1:(length(as_res[["allProfiles_AS"]])-1)) {
    cat(".")
    for(j in (i+1):length(as_res[["allProfiles_AS"]]))
      alldists[i,j] <- ascn_dist(as_res[["allProfiles_AS"]][[i]][["nprof.fixed"]],as_res[["allProfiles_AS"]][[j]][["nprof.fixed"]])
  }
  rownames(alldists) <- names(as_res[["allProfiles_AS"]])
  colnames(alldists) <- names(as_res[["allProfiles_AS"]])
  alldists.dist <- as.dist(t(alldists))
  saveRDS(alldists.dist, paste("MPNST_DLP_asCN_dist.rds"))
  
  #Manhattan clustering
  saveRDS(hclust(alldists.dist, method = "ward.D2"), 
          paste0("MPNST_DLP_asCN_hclust_ward.rds"))#maxime's pairwise distance calc
  MPNST_asCN_hclust_ward <- readRDS(paste0("MPNST_DLP_asCN_hclust_ward.rds"))
  
  as_CN_colour <- c("royalblue3", "skyblue2", "grey80", "white", "gold1", "khaki1", "darkorange3", "darkorange1", "orange", "red4", "red", "orangered2", "purple4")
  as_CN_breaks = c(0,10,20,21,30,31,40,41,42,50,51,52,60,1000)-0.1
  
  #Plot CN profiles
  ha_row = rowAnnotation(df = data.frame(Region = samples[unlist(mclapply(rownames(scDNA_asCN_mtx), function(b) grep(b, lapply(barcodes, function(s) as.character(s$sample_id))), mc.cores = 60))]),
                         col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF")), show_annotation_name = F)
  png(filename = paste0("MPNST_DLP_asfixed_asCN_heatmap.png"), width = 4000, height = 4000, res = 200)
  print(sc_asCN_heatmap(CN_mtx = scDNA_asCN_mtx, hclust = MPNST_asCN_hclust_ward, probes = asCN_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), title = "MPNST DLP \nAS CN "))
  dev.off()
  
  png(filename = paste0("MPNST_DLP_asCN_colourkey.png"), width = 2000, height = 500, res = 200)
  barplot <- barplot(rep(1,length(as_CN_colour)), col = as_CN_colour)
  text(barplot, .2, c("0+0", "1+0", "2+0", "1+1", "3+0", "2+1", "4+0", "3+1", "2+2", "5+0", "4+1", "3+2", ">6"), 0, cex =1, pos=3)
  dev.off()
  
}

####################################################################################################################################
### Part 5: Run multi-PCF on all 10X and DLP cells
####################################################################################################################################

if (T) {
  gamma = 5
  
  if (!file.exists(paste0("all_ASCAT.sc_mpcf_", gamma,".rds"))) {
    system(paste0("mkdir -p ", output.dir, "All_mpcf_", gamma))
    res_list <- c(lapply(samples[1:6], function(s) readRDS(paste0(s, "_10X_ASCAT.sc.rds"))),
                  lapply(samples[1:6], function(s) readRDS(paste0("../../../DLP_plus/results/ASCAT.sc/", s, "_filtered_ASCAT.sc.rds"))))
    # res_list <- c(lapply(samples[1:2], function(s) readRDS(paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/DLP_plus/results/ASCAT.sc/", s, "_ASCAT.sc.rds"))),
    #               lapply(samples[3:6], function(s) readRDS(paste0(s, "_ASCAT.sc.rds"))))
    names(res_list) <- NULL
    
    res_comb <- list(allTracks = do.call(c, lapply(1:length(res_list), function(s) {
      res <- res_list[[s]]$allTracks
      names(res) <- paste0(samples[ifelse(mod(s,6)==0, 6, mod(s,6))], "_", ifelse(s<7, "10X", "DLP"), "_", names(res))
      return(res[res_list[[s]]$filters])
    })))
    timetoread_tumours = 1

    all_res=run_sc_sequencing(tumour_bams=names(res_comb[["allTracks"]]),
                              res = res_comb,
                              allchr=paste0("chr", c(1:22,"X")),
                              sex=rep('male',length(res_comb[["allTracks"]])),
                              chrstring_bam='chr',
                              purs = c(0.5,1),
                              ploidies = seq(1.7,5,0.01),
                              maxtumourpsi=5,
                              binsize=500000,
                              segmentation_alpha = (1/gamma),
                              predict_refit = TRUE,
                              print_results = TRUE,
                              build="hg38",
                              MC.CORES=40,
                              barcodes_10x = NULL, 
                              outdir = paste0("./All_mpcf_",gamma,"/"),  #Subfolder to save files
                              probs_filters = 0.1, 
                              path_to_phases = NULL, 
                              list_ac_counts_paths = NULL,
                              sc_filters = TRUE, 
                              projectname = paste0("MPNST_All"),
                              steps = NULL,
                              smooth_sc = FALSE, 
                              multipcf = TRUE)
    saveRDS(all_res,file=paste0("all_ASCAT.sc_mpcf_", gamma,".rds"))
  } else {
    all_res <- readRDS(paste0("all_ASCAT.sc_mpcf_", gamma,".rds"))
  }
  
  if (!file.exists(paste0("MPNST_all_probes_mpcf_", gamma,".rds"))) {
    all_probes <- do.call(rbind, lapply(all_res$allTracks.processed[[1]][["lSegs"]], function(s) return(s[["output"]]))) %>% 
      select(chrom, num.mark, loc.start, loc.end) %>% dplyr::rename(n.probes = num.mark, startpos = loc.start, endpos = loc.end) %>%
      filter(n.probes>1) #Remove segments of just 1/2 probes(these are NAs in ascat.sc outputs)
    
    all_chr_probes <- rbind(data.frame(chrom = 0, total.probes = 0), all_probes[,1:2] %>% group_by(chrom) %>% summarise(total.probes = sum(n.probes))) %>% 
      mutate(cum.probes = cumsum(total.probes))
    saveRDS(all_probes, paste0("MPNST_all_probes_mpcf_", gamma,".rds"))
    saveRDS(all_chr_probes, paste0("MPNST_all_chr_probes_mpcf_", gamma,".rds"))
  } else {
    all_probes <- readRDS(paste0("MPNST_all_probes_mpcf_", gamma,".rds"))
    all_chr_probes <- readRDS(paste0("MPNST_all_chr_probes_mpcf_", gamma,".rds"))
  }
  
  #Note refit results were mostly identical to without refit so used without as a refit put a few cells into WGD state
  if (!file.exists(paste0("MPNST_all_CN_mtx_mpcf_", gamma,".rds"))) {
    scDNA_CN_mtx <- do.call(rbind, lapply(1:length(all_res$allProfiles), function(b){
      rep(na.omit(as.numeric(all_res$allProfiles[[b]][,"total_copy_number"])), all_probes$n.probes)
    }))
    # names(all_barcodes) <- unlist(lapply(barcodes, function(s) return(as.character(s$sample_id)))) %>% unname()
    # rownames(scDNA_CN_mtx) <- all_barcodes[gsub(".bam", "", names(all_res$allProfiles))]
    rownames(scDNA_CN_mtx) <- names(all_res$allProfiles)
    saveRDS(scDNA_CN_mtx, paste0("MPNST_all_CN_mtx_mpcf_", gamma,".rds"))
  } else {
    scDNA_CN_mtx <- readRDS(paste0("MPNST_all_CN_mtx_mpcf_", gamma,".rds"))
  }
  
  #Generate raw heatmap
  MPNST_CN_hclust_man_ward <- hclust_save_load(scDNA_CN_mtx, sample = paste0("MPNST_all_CN_all_mpcf_", gamma,"_raw"), distance = "manhattan", method = "ward.D2" )
  ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx)),
                                         Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx))),
                         col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
                                    Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")), 
                         annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1),
                                                        Tech = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
  png(filename = paste0("MPNST_all_scDNA_CN_mpcf_", gamma,"_raw.png"), width = 4000, height = 4000, res = 200)
  print(sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx, hclust = MPNST_CN_hclust_man_ward, row_ann = ha_row, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X"))))
  dev.off()
  
  K = 70
  named_clusters <- cutree(MPNST_CN_hclust_man_ward, k = K) %>% sort()
  ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[names(named_clusters),])),
                                         Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx[names(named_clusters),]))),
                         col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
                                    Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")), show_annotation_name = F)
  png(filename = paste0("MPNST_all_scDNA_CN_mpcf_", gamma,"_raw_hclust_K", K, ".png"), width = 4000, height = 4000, res = 200)
  print(sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), row_ann = ha_row, title = "MPNST All \nTotal CN "))
  dev.off()

  #Generate pass heatmap
  MPNST_CN_hclust_man_ward <- hclust_save_load(scDNA_CN_mtx[all_res$filters,], sample = paste0("MPNST_all_CN_all_mpcf_", gamma,"_pass"), distance = "manhattan", method = "ward.D2" )
  ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[all_res$filters,])),
                                         Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx[all_res$filters,]))),
                         col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
                                    Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")), show_annotation_name = F)
  png(filename = paste0("MPNST_all_scDNA_CN_mpcf_", gamma,"_pass.png"), width = 4000, height = 4000, res = 200)
  print(sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[all_res$filters,], hclust = MPNST_CN_hclust_man_ward, row_ann = ha_row, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), title = "MPNST All \nTotal CN "))
  dev.off()
  
  #Filter actually doesn't work well given it is two distributions (10X and DLP)
  MPNST_CN_hclust_man_ward <- hclust_save_load(scDNA_CN_mtx, sample = paste0("MPNST_all_CN_all_mpcf_", gamma,"_raw"), distance = "manhattan", method = "ward.D2" )
  nonwgd_cells <- names(cutree(MPNST_CN_hclust_man_ward, k = 10)[cutree(MPNST_CN_hclust_man_ward, k = 10) == 1]) #Cluster 2 and 3 is WGD cells
  
  #Check quality of cells fitted with extra WGD
  wgd_cells <- names(cutree(MPNST_CN_hclust_man_ward, k = 10)[cutree(MPNST_CN_hclust_man_ward, k = 10) %in% c(2:3)]) #Cluster 2 and 3 is WGD cells
  
  #Plot on nreads/noise plot
  if (F) {
    allT <- all_res$allTracks.processed
    allS <- all_res$allSolutions
    getloess <- function(qu, nr) {
      nms <- paste0("n", 1:length(qu))
      names(qu) <- nms
      quo <- qu[order(nr, decreasing = F)]
      fitted <- stats::runmed(quo, k = 31, endrule = "keep")
      names(fitted) <- names(quo)
      list(fitted = fitted[nms], residuals = quo[nms] - 
             fitted[nms])
    }
    getQuality.SD <- function(allT) {
      sapply(allT, function(x) {
        median(abs(diff(unlist(lapply(x$lCTS, function(y) y$smoothed)))))
      })
    }
    nrecords <- sapply(allT, function(x) sum(unlist(lapply(x$lCTS, 
                                                           function(y) y$records))))
    # thresholdNrec <- quantile(nrecords, probs = probs)
    # ambiguous <- sapply(allS, function(x) x$ambiguous)
    # doublet <- sapply(allS, function(x) if (!is.null(x$bestfit)) 
    #   !x$bestfit$ambiguous
    #   else F)
    qualities <- getQuality.SD(allT)
    # thresholdQual <- quantile(qualities, probs = 1 - probs)
    # keep <- qualities <= thresholdQual & nrecords >= thresholdNrec & 
    #   !ambiguous & !doublet
    # keep2 <- !(qualities < thresholdQual & nrecords < thresholdNrec)
    # keep2 <- keep2 & !ambiguous
    # ll <- getloess(qualities[keep2], log2(nrecords)[keep2])
    # ord <- order(log2(nrecords)[keep2], decreasing = F)
    # filters <- (1:length(nrecords)) %in% (which(keep2)[abs(ll$residuals) <= 0.02]) & keep
    
    #Plot
    pdf(paste0("MPNST_WGD_log2nrecord_vs_qualities.pdf"))
    plot(qualities, log2(nrecords), xlab = "Noise logr", ylab = "Total number of reads", pch = 19, 
         cex = ifelse(names(allT) %in% wgd_cells, 0.3, 0.1),
         col = ifelse(names(allT) %in% wgd_cells, "red", "black"))
    dev.off()
  }

  #Get kmeans clusters
  K = 22
  
  png(filename = paste0("MPNST_all_scDNA_CN_mpcf_", gamma,"_K", K, ".png"), width = 4000, height = 4000, res = 200)
  if (!file.exists(paste0("MPNST_all_k_means_K",K,"_clusters.rds"))) {
    set.seed(123)
    hm <- ComplexHeatmap::Heatmap(matrix = scDNA_CN_mtx[nonwgd_cells,], cluster_rows = T, cluster_row_slices = FALSE, row_km = K, row_km_repeats = 1000)
    clusters <- row_order(hm)
    k_means_clusters <- lapply(1:length(clusters), function(j) {
      return(rownames(scDNA_CN_mtx[nonwgd_cells,])[clusters[[j]]])
    })
    k_means_clusters_hclust <- lapply(k_means_clusters, function(k) {
      if(length(k) > 1) {
        hclust <- fastcluster::hclust(dist(scDNA_CN_mtx[k,], method = "manhattan"), method = "ward.D2")
        return(hclust[["labels"]][hclust[["order"]]])
      } else {
        return(k)
      }
    }) #perform hclust on each kmeans cluster
    saveRDS(k_means_clusters, paste0("MPNST_all_k_means_K",K,"_clusters.rds"))
    saveRDS(k_means_clusters_hclust, paste0("MPNST_all_k_means_K",K,"_clusters_hclust.rds"))
  } else {
    k_means_clusters <- readRDS(paste0("MPNST_all_k_means_K",K,"_clusters_hclust.rds")) #using hclust on each kmeans
  }
  named_clusters <- unlist(lapply(1:length(k_means_clusters), function(i) return(rep(i, length(k_means_clusters[[i]]))))) %>% set_names(unlist(k_means_clusters))
  ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[names(named_clusters),])),
                                         Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx[names(named_clusters),]))),
                         col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
                                    Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")), 
                         annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1),
                                                        Tech = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
  sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), row_ann = ha_row)
  dev.off()
  
  #Plot profiles of wgd cells and mixed cluster
  if (F) {
    cluster_name <- "WGD"
    idx_to_print <- which(names(all_res$allProfiles) %in% wgd_cells)
    
    cluster_name <- "kmeans_2"
    idx_to_print <- which(names(all_res$allProfiles) %in% k_means_clusters[[2]])
    
    pdf(paste0(cluster_name, "_profiles.pdf"), width = 15, height = 4)
    for (i in idx_to_print) {
      try({
        plotSolution(all_res$allTracks.processed[[i]], purity = all_res$allSolutions[[i]]$purity, 
                     ploidy = all_res$allSolutions[[i]]$ploidy, gamma = 1, 
                     sol = all_res$allSolutions[[i]])
        title(names(all_res$allTracks)[i])
      })
    }
    dev.off()
  }
  
  #Named kmeans of scDNA_CN_cluster_ids (cleaned kmeans)
  png(filename = paste0("MPNST_all_scDNA_CN_mpcf_", gamma,"_K", K, "_named.png"), width = 4000, height = 4000, res = 200)
  named_clusters <- unlist(lapply(1:length(scDNA_CN_cluster_ids), function(i) return(rep(names(scDNA_CN_cluster_ids)[i], length(scDNA_CN_cluster_ids[[i]]))))) %>% set_names(unlist(scDNA_CN_cluster_ids))
  ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[names(named_clusters),])),
                                         Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx[names(named_clusters),]))),
                         col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
                                    Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")),
                         annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1),
                                                        Tech = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
  sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), row_ann = ha_row)
  dev.off()
  
  #Plot only 10X or DLP
  subset = "DLP"
  png(filename = paste0("MPNST_", subset, "_scDNA_CN_mpcf_", gamma,"_K", K, ".png"), width = 4000, height = 4000, res = 200)
  k_means_clusters <- readRDS(paste0("MPNST_all_k_means_K",K,"_clusters.rds"))
  named_clusters <- unlist(lapply(1:length(scDNA_CN_cluster_ids), function(i) return(rep(names(scDNA_CN_cluster_ids)[i], length(scDNA_CN_cluster_ids[[i]]))))) %>% set_names(unlist(scDNA_CN_cluster_ids))
  named_clusters <- named_clusters[grep(subset, names(named_clusters))]
  ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[names(named_clusters),])),
                                         Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx[names(named_clusters),]))),
                         col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
                                    Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")),
                         annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1),
                                                        Tech = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
  sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), row_ann = ha_row)
  dev.off()
}

####################################################################################################################################
### Part 6: Get allele specifc calls for all 10X and DLP cells
####################################################################################################################################

if (T) {
  gamma = 5
  system(paste0("mkdir -p ", output.dir, "ascn_mpcf_", gamma))
  
  if (!file.exists(paste0("all_ASCAT.sc_ascn_mpcf_", gamma,".rds"))) {
    all_res <- readRDS(file=paste0("all_ASCAT.sc_mpcf_", gamma,".rds"))
    
    path_to_phases <- list(paste0("/camp/project/proj-vanloo/analyses/hyan/mpnst/bulk/results/BB_spiked/All_with_Ext/",
                                  "normal_multisample_haplotypes_chr",c(1:22,"23"),".vcf"))
    DLP_ac_path <- "/camp/project/proj-vanloo/analyses/hyan/mpnst/DLP_plus/results/ASCAT.sc/asCN/Allele_Counts/"
    tenX_ac_path <- "/camp/project/proj-vanloo/analyses/hyan/mpnst/10X_DLP/results/ASCAT.sc/asCN/Allele_Counts/"
    
    #Turn allele freq into per cell for 10X (DLP already done)
    if (F) {
      system(paste0("mkdir -p ", tenX_ac_path))
      for (s in samples[1:6]) {
        for (c in 1:23) {
          print(c)
          allele_counts <- read.delim(paste0("../../../10X_DNA/results/Genotype_SNPs/Allele_Freq/",s,"/", s,"_alleleFrequencies_chr",c,".txt"), stringsAsFactors = F)
          colnames(allele_counts)[1] <- "#CHR"
          pbmclapply(grep(paste0(s, "_10X"), names(all_res$allTracks), value = T), function(b) {
            write.table(allele_counts %>% filter(Barcode == gsub(".*_", "", b)), file = paste0(tenX_ac_path, b, "_alleleFrequencies_chr", c, ".txt"), quote = F, row.names = F, sep = "\t")
          }, mc.cores = 6) #don't do too many cores, runs out of memory
        }
      }
    }
    
    list_ac_counts_paths <- c(lapply(grep("10X", names(all_res$allTracks), value = T), function(b){
      paste0(tenX_ac_path, b, "_alleleFrequencies_chr", c(1:22,"23"), ".txt")
    }) %>% set_names(grep("10X", names(all_res$allTracks), value = T)),
    lapply(grep("DLP", names(all_res$allTracks), value = T), function(b){
      paste0(DLP_ac_path, gsub(".*_DLP_(.*).bam", "\\1", b) , "_alleleFrequencies_chr", c(1:22,"23"), ".txt")
    }) %>% set_names(grep("DLP", names(all_res$allTracks), value = T)))
    
    if (F) {
      # as_res <- getAS_CNA(all_res, path_to_phases = path_to_phases, 
      #                     list_ac_counts_paths = list_ac_counts_paths[], purs = c(0.5,1), 
      #                     ploidies = seq(1.7,5,0.01), outdir = paste0("./ascn_mpcf_",gamma,"/"), projectname = paste0("MPNST_DLP_All"), 
      #                     steps = NULL, mc.cores = 10)
      # names(as_res$allProfiles_AS) <- names(as_res$allProfiles)
      #Ran split as kept crashing
      
      suppressPackageStartupMessages(require(GenomicRanges))
      suppressPackageStartupMessages(require(data.table))
      .searchGrid <- function(baf, logr, sizes, purs = seq(0.9, 
                                                           0.99, 0.01), ploidies = seq(2, 6, 0.01)) {
        getNANB <- function(baf, logr, purity, ploidy) {
          K <- (2^logr * ploidy - 2 * (1 - purity))/purity
          Na <- (1 - purity - baf * (2 - 2 * purity) - (purity * 
                                                          baf - purity) * K)/purity
          Nb <- K - Na
          list(Na = Na, Nb = Nb)
        }
        errors <- function(baf, logr, sizes, purity, ploidy) {
          nanb <- getNANB(baf, logr, purity, ploidy)
          ssize <- sum(sizes)
          sum(((nanb$Na - round(nanb$Na))^2 + (nanb$Nb - round(nanb$Nb))^2) * 
                (sizes/ssize))
        }
        ggprofile <- function(baf, logr, purity, ploidy) {
          nanb <- getNANB(baf, logr, purity, ploidy)
        }
        errs <- t(sapply(purs, function(rho) sapply(ploidies, 
                                                    function(psi) errors(baf = baf, logr = logr, sizes, 
                                                                         purity = rho, ploidy = psi))))
        rownames(errs) <- purs
        colnames(errs) <- ploidies
        purs <- as.numeric(rownames(errs))
        ploidies <- as.numeric(colnames(errs))
        mins <- arrayInd(which.min(errs), dim(errs))
        purity <- purs[mins[1]]
        ploidy <- ploidies[mins[2]]
        return(list(purity = purity, ploidy = ploidy, profile = ggprofile(baf, 
                                                                          logr, purity, ploidy), errs = errs))
      }
      readPhases <- function(phasing_paths) {
        phasing <- lapply(phasing_paths, function(x) as.data.frame(data.table::fread(x)))
        phases <- lapply(phasing, function(x) {
          x <- x[grep("0\\|1|1\\|0", x[, 10]), ]
          phase <- gsub("(.*)\\|(.*)", "\\1", x[, 10])
          phases1 <- x[, "REF"]
          phases1[phase == "1"] <- x[phase == "1", "ALT"]
          phases2 <- x[, "ALT"]
          phases2[phases2 == phases1] <- x[phases2 == phases1, 
                                           "REF"]
          list(chr = x[, 1], pos = x[, 2], phases1 = phases1, 
               phases2 = phases2)
        })
        phasesall <- do.call("rbind", lapply(phases, function(x) data.frame(chr = x[[1]], 
                                                                            pos = x[[2]], phase1 = x[[3]], phase2 = x[[4]])))
      }
      getPhasedInfo <- function(ac, phasing) {
        ac[, 1] <- gsub("chr", "", ac[, 1])
        phasing[, 1] <- gsub("chr", "", phasing[, 1])
        ac.. <- merge(ac, phasing, by.x = c("#CHR", "POS"), by.y = c("chr", 
                                                                     "pos"))
        nmsPH1 <- ac..[, "phase1"]
        nmsPH2 <- ac..[, "phase2"]
        letters <- c("A", "C", "G", "T")
        letters_index <- sapply(paste0("Count_", letters), function(x) which(colnames(ac..) %in% 
                                                                               x))
        names(letters_index) <- letters
        ind1 <- letters_index[nmsPH1]
        ind2 <- letters_index[nmsPH2]
        ind1 <- cbind(seq_along(ind1), ind1)
        ind2 <- cbind(seq_along(ind2), ind2)
        df <- data.frame(chr = ac..[, 1], pos = ac..[, 2], counts1 = as.numeric(as.character(ac..[ind1])), 
                         counts2 = as.numeric(as.character(ac..[ind2])))
        df <- df[rowSums(df[, 3:4]) > 0, ]
        df
      }
      getHet <- function(snp1, snp2) {
        gr1 <- GRanges(gsub("chr", "", snp1[, 1]), IRanges(snp1[, 
                                                                2], snp1[, 2]))
        gr2 <- GRanges(gsub("chr", "", snp2[, 1]), IRanges(snp2[, 
                                                                2], snp2[, 2]))
        ovs <- findOverlaps(gr1, gr2)
        snp1[queryHits(ovs), ]
      }
      getAC <- function(ac_counts_paths, phases) {
        acs <- do.call("rbind", lapply(ac_counts_paths, function(dd) {
          gc()
          ok <- getHet(as.data.frame(data.table::fread(dd)), 
                       phases)
        }))
      }
      fitBinom.1dist <- function(counts, depths, steps = NULL) {
        if (any(depths < counts)) {
          print("Anomaly: total counts smaller than genotype counts - verify allele count files")
          depths[depths > counts] <- counts[depths > counts]
        }
        if (is.null(steps)) 
          steps <- if (length(counts)%/%3 > 10) 
            5
        else 3
        nonas <- !is.na(counts) & !is.na(depths)
        counts <- counts[nonas]
        depths <- depths[nonas]
        haploblocks <- cut(1:length(counts), pmax(length(counts)%/%steps, 
                                                  2))
        if (length(counts) < 10) 
          haploblocks <- rep(1, length(counts))
        lcounts <- split(counts, haploblocks)
        ldepths <- split(depths, haploblocks)
        lcounts <- lapply(1:length(lcounts), function(x) if (rnorm(1) < 
                                                             0) 
          ldepths[[x]] - lcounts[[x]]
          else lcounts[[x]])
        counts <- sapply(lcounts, sum)
        depths <- sapply(ldepths, sum)
        values <- seq(0.5, 1, 0.001)
        llh <- sapply(values, function(x) {
          sum(log(dbinom(counts, size = depths, prob = x, log = F) + 
                    dbinom(counts, size = depths, prob = 1 - x, log = F)))
        })
        llh <- llh - max(llh)
        normalised <- exp(llh)/sum(exp(llh))
        baf <- values[which.max(llh)]
        cs <- cumsum(normalised)
        q95 <- values[c(which(cs > 0.05)[1], which(cs > 0.95)[1])]
        c(q95[1], baf, q95[2])
      }
      fitBinom.1dist.noswitch <- function(counts, depths) {
        nonas <- !is.na(counts) & !is.na(depths)
        counts <- counts[nonas]
        depths <- depths[nonas]
        values <- seq(0, 1, 0.001)
        llh <- sapply(values, function(x) {
          sum(dbinom(counts, size = depths, prob = x, log = T))
        })
        llh <- llh - max(llh)
        normalised <- exp(llh)/sum(exp(llh))
        baf <- values[which.max(llh)]
        cs <- cumsum(normalised)
        q95 <- values[c(which(cs > 0.05)[1], which(cs > 0.95)[1])]
        c(q95[1], baf, q95[2])
      }
      getProfile <- function(df, prof, purity, ploidy, purs, ploidies, 
                             steps = NULL) {
        nprof <- data.frame(chr = as.character(prof[, "chromosome"]), 
                            startpos = as.numeric(prof[, "start"]), endpos = as.numeric(prof[, 
                                                                                             "end"]), total_copy_number = as.numeric(prof[, 
                                                                                                                                          "total_copy_number"]), total_copy_number_logr = as.numeric(prof[, 
                                                                                                                                                                                                          "total_copy_number_logr"]), logr = as.numeric(prof[, 
                                                                                                                                                                                                                                                             "logr"]), logr.sd = as.numeric(prof[, "logr.sd"]), 
                            fitted = as.numeric(prof[, "total_copy_number"]), 
                            q05 = as.numeric(NA), BAF = as.numeric(NA), q95 = as.numeric(NA), 
                            q05_noswitch = as.numeric(NA), BAF_noswitch = as.numeric(NA), 
                            q95_noswitch = as.numeric(NA), stringsAsFactors = F)
        nprof[nprof[, 1] == "chr23", 1] <- "chrX"
        nprof[nprof[, 1] == "23", 1] <- "X"
        nprof[nprof[, 1] == "chr24", 1] <- "chrY"
        nprof[nprof[, 1] == "24", 1] <- "Y"
        grseg <- GRanges(gsub("chr", "", nprof[, "chr"]), IRanges(as.numeric(nprof[, 
                                                                                   "startpos"]), as.numeric(nprof[, "endpos"])))
        grsnp <- GRanges(gsub("chr", "", df[, 1]), IRanges(df[, 
                                                              2], df[, 2]))
        ovs <- findOverlaps(grseg, grsnp)
        qH <- queryHits(ovs)
        sH <- subjectHits(ovs)
        df[, 3] <- as.numeric(as.character(df[, 3]))
        df[, 4] <- as.numeric(as.character(df[, 4]))
        for (i in unique(qH)) {
          inds <- 1:nrow(df) %in% sH[qH == i]
          nprof[i, c("q05", "BAF", "q95")] <- fitBinom.1dist(df[inds, 
                                                                3], rowSums(df[inds, c(3, 4)]), steps = steps)
          nprof[i, c("q05_noswitch", "BAF_noswitch", "q95_noswitch")] <- fitBinom.1dist.noswitch(df[inds, 
                                                                                                    3], rowSums(df[inds, c(3, 4)]))
        }
        nona <- !is.na(nprof[, "BAF"])
        sG <- .searchGrid(as.numeric(nprof[nona, "BAF"]), as.numeric(nprof[nona, 
                                                                           "logr"]), as.numeric(nprof[nona, "endpos"]) - as.numeric(nprof[nona, 
                                                                                                                                          "startpos"]), purs = purs, ploidies = ploidies)
        sG.fixed <- .searchGrid(nprof[nona, "BAF"], nprof[nona, 
                                                          "logr"], nprof[nona, "endpos"] - nprof[nona, "startpos"], 
                                purs = purity, ploidies = ploidy)
        exceptNA <- function(vec, nona) {
          newvec <- rep(NA, length(nona))
          newvec[nona] <- vec
          newvec
        }
        nprof[, "ntot_free"] <- transform_bulk2tumour(nprof[, 
                                                            "logr"], sG$purity, sG$ploidy)
        nprof[, "ntot_fixed"] <- transform_bulk2tumour(nprof[, 
                                                             "logr"], sG$purity, sG$ploidy)
        list(nprof.free = cbind(nprof, nA = round(nprof[, "fitted"] * 
                                                    nprof[, "BAF"]), nB = nprof[, "fitted"] - round(nprof[, 
                                                                                                          "fitted"] * nprof[, "BAF"]), nA_sc_raw = exceptNA(sG$profile$Nb, 
                                                                                                                                                            nona), nB_sc_raw = exceptNA(sG$profile$Na, nona), 
                                nA_sc = exceptNA(round(sG$profile$Nb), nona), nB_sc = exceptNA(round(sG$profile$Na), 
                                                                                               nona)), purity.free = sG$purity, ploidy.free = sG$ploidy, 
             nprof.fixed = cbind(nprof, nA = round(nprof[, "fitted"] * 
                                                     nprof[, "BAF"]), nB = nprof[, "fitted"] - round(nprof[, 
                                                                                                           "fitted"] * nprof[, "BAF"]), nA_sc_raw = exceptNA(sG.fixed$profile$Nb, 
                                                                                                                                                             nona), nB_sc_raw = exceptNA(sG.fixed$profile$Na, 
                                                                                                                                                                                         nona), nA_sc = exceptNA(round(sG.fixed$profile$Nb), 
                                                                                                                                                                                                                 nona), nB_sc = exceptNA(round(sG.fixed$profile$Na), 
                                                                                                                                                                                                                                         nona)), purity.fixed = sG.fixed$purity, ploidy.fixed = sG.fixed$ploidy)
      }
      getAS_CNA_sample <- function(track, profile, ac_counts_paths, 
                                   purs, ploidies, purity, ploidy, phases = NULL, path_to_phases = NULL, 
                                   steps = NULL) {
        if (is.null(phases)) 
          phases <- readPhases(path_to_phases)
        ac <- getAC(ac_counts_paths, phases)
        ac.ph <- getPhasedInfo(ac, phases)
        prof <- getProfile(ac.ph, prof = profile, steps = steps, 
                           purity = purity, ploidy = ploidy, purs = purs, ploidies = ploidies)
        prof
      }
      phases <- NULL
      if (length(path_to_phases) == 1) {
        print("## read-in Phases")
        phases <- readPhases(path_to_phases[[1]])
      }
      print("## derive Allele-specific Profiles")
      #Run allele specific calling in batches
      cells <- split(1:length(all_res$allTracks.processed), ceiling(seq_along(1:length(all_res$allTracks.processed))/150))
      for (i in 1:length(cells)) {
        allProfiles_AS <- pbmclapply(cells[[i]], 
                                     function(x) {
                                       getAS_CNA_sample(track = all_res$allTracks.processed[[x]], 
                                                        profile = all_res$allProfiles[[x]], ac_counts_paths = list_ac_counts_paths[[x]], 
                                                        phases = phases, purity = if (any(grepl("refitted", 
                                                                                                names(all_res)))) 
                                                          all_res$allSolutions.refitted.auto[[x]]$purity
                                                        else all_res$allSolutions[[x]]$purity, ploidy = if (any(grepl("refitted", 
                                                                                                                      names(all_res)))) 
                                                          all_res$allSolutions.refitted.auto[[x]]$ploidy
                                                        else all_res$allSolutions[[x]]$ploidy, purs = c(0.5,1), 
                                                        ploidies = seq(1.7,5,0.01), path_to_phases = if (length(path_to_phases) > 
                                                                                                         1) 
                                                          path_to_phases[[x]]
                                                        else NULL, steps = NULL)
                                     }, mc.cores = 15)
        saveRDS(allProfiles_AS, paste0("allProfiles_AS_", min(cells[[i]]), "_", max(cells[[i]]), ".rds"))
      }
      all_res$allProfiles_AS <- do.call(c, lapply(1:length(cells), function(i) readRDS(paste0("allProfiles_AS_", min(cells[[i]]), "_", max(cells[[i]]), ".rds")))) #combine
      
      print("## plot Allele-specific Profiles")
      pdf(paste0("all_as_cna_profiles.pdf"), width = 15, height = 5)
      tnull <- lapply(1:length(all_res$allProfiles_AS), function(x) {
        try({
          plot_AS_profile(all_res$allProfiles_AS[[x]]$nprof.fixed)
          title(paste0(names(all_res$allTracks)[x], " - bam", x), 
                cex = 0.5)
        })
      })
      dev.off()
      res
    }
    names(all_res$allProfiles_AS) <- names(all_res$allProfiles)
    #3106
    
    # as_res <- run_sc_sequencing(tumour_bams=names(all_res[["allTracks"]]),
    #                             res = all_res,
    #                             allchr=paste0("chr", c(1:22,"X")),
    #                             sex=rep('male',length(all_res[["allTracks"]])),
    #                             chrstring_bam='chr',
    #                             purs = c(0.5,1),
    #                             ploidies = seq(1.7,5,0.01),
    #                             maxtumourpsi=5,
    #                             binsize=500000,
    #                             segmentation_alpha = (1/gamma),
    #                             predict_refit = FALSE,
    #                             print_results = FALSE,
    #                             build="hg38",
    #                             MC.CORES=40,
    #                             barcodes_10x = NULL, 
    #                             outdir = paste0("./ascn_mpcf_",gamma,"/"),  #Subfolder to save files
    #                             probs_filters = 0.1, 
    #                             path_to_phases = path_to_phases, 
    #                             list_ac_counts_paths = list_ac_counts_paths[names(all_res[["allTracks"]])],
    #                             sc_filters = TRUE, 
    #                             projectname = paste0("MPNST_DLP_All"),
    #                             steps = NULL,
    #                             smooth_sc = FALSE, 
    #                             multipcf = TRUE)
    as_res <- all_res
    saveRDS(as_res,file=paste0("all_ASCAT.sc_ascn_mpcf_", gamma,".rds"))
    as_res <- getAS_CNA_smoothed(as_res, mc.cores = 30) #didn't work so not run for now
    
  }
  else {
    as_res <- readRDS(paste0("all_ASCAT.sc_ascn_mpcf_", gamma,".rds"))
  }
  
  if (!file.exists(paste0("MPNST_all_asCN_probes_mpcf_", gamma,".rds"))) {
    asCN_probes <- do.call(rbind, lapply(as_res$allTracks.processed[[1]][["lSegs"]], function(s) return(s[["output"]]))) %>% 
      select(chrom, num.mark, loc.start, loc.end) %>% dplyr::rename(n.probes = num.mark, startpos = loc.start, endpos = loc.end) %>%
      filter(n.probes>1) #Remove segments of just 1/2 probes(these are NAs in ascat.sc outputs)
    
    asCN_chr_probes <- rbind(data.frame(chrom = 0, total.probes = 0), asCN_probes[,1:2] %>% group_by(chrom) %>% summarise(total.probes = sum(n.probes))) %>% 
      mutate(cum.probes = cumsum(total.probes))
    saveRDS(asCN_probes, paste0("MPNST_all_asCN_probes_mpcf_", gamma,".rds"))
    saveRDS(asCN_chr_probes, paste0("MPNST_all_asCN_chr_probes_mpcf_", gamma,".rds"))
  } else {
    asCN_probes <- readRDS(paste0("MPNST_all_asCN_probes_mpcf_", gamma,".rds"))
    asCN_chr_probes <- readRDS(paste0("MPNST_all_asCN_chr_probes_mpcf_", gamma,".rds"))
  }
  
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
  
  if (!file.exists("MPNST_all_asCN_hclust_ward.rds")) {
    #Function to calculate allele specific distance
    ascn_dist <- function(cna1, cna2) {
      # ind1 <- is.na(cna1[,"nA"]) | is.na(cna2[,"nA"])
      # ind2 <- !ind1 
      ind1 <- (is.na(cna1[,"nA"]) | is.na(cna2[,"nB"])) #Edited as one segment has NA
      ind2 <- !(is.na(cna1[,"nA"]) | is.na(cna2[,"nB"]))
      sizes <- cna1[,"endpos"]-cna1[,"startpos"]
      sizes <- sizes/500000/sum(sizes/500000)
      dists <- sum(pmin(2,abs(cna1[ind2,"nA"]-cna2[ind2,"nA"])+abs(cna1[ind2,"nB"]-cna2[ind2,"nB"]))*sizes[ind2])
      # dists <- dists+sum(pmin(2,abs(cna1[ind1,"fitted"]-cna2[ind1,"fitted"])*sizes[ind1]))
      dists
    }
    
    alldists <- matrix(NA,length(as_res[["allProfiles_AS"]]),length(as_res[["allProfiles_AS"]]))
    for(i in 1:(length(as_res[["allProfiles_AS"]])-1)) {
      cat(".")
      for(j in (i+1):length(as_res[["allProfiles_AS"]]))
        alldists[i,j] <- ascn_dist(as_res[["allProfiles_AS"]][[i]][["nprof.fixed"]],as_res[["allProfiles_AS"]][[j]][["nprof.fixed"]])
    }
    rownames(alldists) <- names(as_res[["allProfiles_AS"]])
    colnames(alldists) <- names(as_res[["allProfiles_AS"]])
    alldists.dist <- as.dist(t(alldists))
    saveRDS(alldists.dist, paste("MPNST_all_asCN_dist.rds"))
    #Manhattan clustering
    saveRDS(hclust(alldists.dist, method = "ward.D2"), 
            paste0("MPNST_all_asCN_hclust_ward.rds"))#maxime's pairwise distance calc
    MPNST_asCN_hclust_ward <- readRDS(paste0("MPNST_all_asCN_hclust_ward.rds"))
  } else {
    MPNST_asCN_hclust_ward <- readRDS(paste0("MPNST_all_asCN_hclust_ward.rds"))
  }
  
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
  
  # named_clusters <- unlist(lapply(1:length(k_means_clusters), function(i) return(rep(i, length(k_means_clusters[[i]]))))) %>% set_names(unlist(k_means_clusters))
  # ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[names(named_clusters),])),
  #                                        Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx[names(named_clusters),]))),
  #                        col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
  #                                   Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")), 
  #                        annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1),
  #                                                       Tech = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
  # sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), row_ann = ha_row)
  # dev.off()
  # 
  # 
  # #Named kmeans of scDNA_CN_cluster_ids (cleaned kmeans)
  # png(filename = paste0("MPNST_all_scDNA_CN_mpcf_", gamma,"_K", K, "_named.png"), width = 4000, height = 4000, res = 200)
  # named_clusters <- unlist(lapply(1:length(scDNA_CN_cluster_ids), function(i) return(rep(names(scDNA_CN_cluster_ids)[i], length(scDNA_CN_cluster_ids[[i]]))))) %>% set_names(unlist(scDNA_CN_cluster_ids))
  # ha_row = rowAnnotation(df = data.frame(Region = gsub("_.*", "", rownames(scDNA_CN_mtx[names(named_clusters),])),
  #                                        Tech = gsub(".*_(10X|DLP)_.*", "\\1", rownames(scDNA_CN_mtx[names(named_clusters),]))),
  #                        col = list(Region = c("R1" = "#B79F00", "R2" = "#00BA38", "R3" = "#00BFC4", "R4" = "#619CFF", "R5" = "#F564E3", "P" = "#F8766D"),
  #                                   Tech = c("10X" = "mediumpurple1", "DLP" = "olivedrab3")),
  #                        annotation_legend_param = list(Region = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1),
  #                                                       Tech = list(title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 18), grid_height = unit(0.8, "cm"), grid_width = unit(0.8, "cm"), gap = unit(2, "cm"), nrow = 1)), show_annotation_name = F)
  # sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx[names(named_clusters),], hclust = F, row_split = named_clusters, probes = all_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X")), row_ann = ha_row)
  # dev.off()
  
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

}
