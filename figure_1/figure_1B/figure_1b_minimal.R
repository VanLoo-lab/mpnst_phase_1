### title: minimal code for figure 1b


### Data import




### Make the final figure

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











