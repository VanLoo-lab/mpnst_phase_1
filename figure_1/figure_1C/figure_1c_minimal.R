### title: minimal code for figure 1c

### Data import

### Make the final figure
clusters_to_plot <- c(3, 1, 12, 10, "Subclone_P", 9, 6, "Subclone_R1", 7, "Subclone_R5", 11, 13, 5, "Subclone_R2", 4, 2, 8, "Subclone_R3")
cluster_colours <- c("black", "grey70", "grey30", "#F8766D", "pink", "coral", "#B79F00", "gold", "#F564E3", "mediumpurple1", 
                     "darkslategray4", "darkgreen", "#00BA38", "chartreuse", "#619CFF", "deepskyblue", "#00BFC4", "cyan")
names(cluster_colours) <- clusters_to_plot
consensus_clusters_CCF_k_plot <- do.call(rbind, lapply(1:length(samples), function(r) {
  consensus_clusters_CCF_k[[r]] %>% mutate(region = samples[r], cluster = consensus_clusters_k_mod$cluster) %>% 
    arrange(factor(cluster, levels = clusters_to_plot)) %>% rowid_to_column(var="ID")
})) %>% filter(cluster %in% clusters_to_plot)


png(filename = paste0("MPNST_all_region_all_clusters_CCF.png"), width = 12000, height = 6000, res = 200)
consensus_clusters_CCF_k_plot %>% arrange(region) %>% 
  ggplot(aes(x = ID, y = subclonal.fraction, fill = as.factor(cluster))) + geom_bar(stat = "identity", width = 1) + scale_fill_manual(values = cluster_colours) +
  scale_y_continuous(breaks = c(0,0.5,1)) + coord_cartesian(ylim=c(0, 1)) + xlab("Variant") + ylab("CCF") + facet_wrap(~region, ncol = 1, strip.position = "right") +
  theme(text=element_text(size=60), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.y.left = element_line(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background = element_blank(), panel.spacing = unit(2, "lines"))
dev.off()


clusters_to_plot <- c(3, 10, "Subclone_P", 9, 6, "Subclone_R1", 7, "Subclone_R5", 11, 13, 5, "Subclone_R2", 4, 2, 8, "Subclone_R3")
cluster_colours <- c("black", "#F8766D", "pink", "coral", "#B79F00", "gold", "#F564E3", "mediumpurple1", 
                     "darkslategray4", "darkgreen", "#00BA38", "chartreuse", "#619CFF", "deepskyblue", "#00BFC4", "cyan")
names(cluster_colours) <- clusters_to_plot
consensus_clusters_CCF_k_plot <- do.call(rbind, lapply(1:length(samples), function(r) {
  consensus_clusters_CCF_k[[r]] %>% mutate(region = samples[r], cluster = consensus_clusters_k_mod$cluster) %>% 
    arrange(factor(cluster, levels = clusters_to_plot)) %>% rowid_to_column(var="ID")
})) %>% filter(cluster %in% clusters_to_plot)


png(filename = paste0("MPNST_all_region_good_clusters_CCF.png"), width = 12000, height = 6000, res = 200)
consensus_clusters_CCF_k_plot %>% arrange(region) %>% 
  ggplot(aes(x = ID, y = subclonal.fraction, fill = as.factor(cluster))) + geom_bar(stat = "identity", width = 1) + scale_fill_manual(values = cluster_colours) +
  scale_y_continuous(breaks = c(0,0.5,1)) + coord_cartesian(ylim=c(0, 1)) + xlab("Variant") + ylab("CCF") + facet_wrap(~region, ncol = 1, strip.position = "right") +
  theme(text=element_text(size=60), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.y.left = element_line(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background = element_blank(), panel.spacing = unit(2, "lines"))
dev.off() 

# pdf(file = paste0("MPNST_all_region_good_clusters_CCF.pdf"), width = 14, height = 7)
# consensus_clusters_CCF_k_plot %>% arrange(region) %>% 
#   ggplot(aes(x = ID, y = subclonal.fraction, fill = as.factor(cluster))) + geom_bar(stat = "identity", width = 1) + scale_fill_manual(values = cluster_colours) +
#   scale_y_continuous(breaks = c(0,1)) + coord_cartesian(ylim=c(0, 1)) + xlab("Variant") + ylab("CCF") + facet_wrap(~region, ncol = 1, strip.position = "right") +
#   theme(text=element_text(size=28), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.y.left = element_line(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background = element_blank(), panel.spacing = unit(2, "lines"))
# dev.off() 
