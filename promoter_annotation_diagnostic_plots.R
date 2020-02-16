###promoter annotation diagnostic plots
canonical_polr2a <- sapply(proms_list_with_canonical, function(xx) {
  xx$polr2a_peak_count[which(xx$canonical == TRUE)]
})
max_polr2a <- sapply(proms_list_with_canonical, function(xx) {
  xx$polr2a_peak_count[which(xx$is_max_polr2a_signal == TRUE)[1]]
})



quartz(file = "canonical_polr2a_chip_signal.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4, 4, 1, 1)+0.1)
hist(canonical_polr2a, col = alpha("red", 0.6), lty = 0, cex.lab = 1.15,
     breaks = 50, freq = FALSE, xlab = "# of ENCODE experiments\nwith POLR2A peak", 
     yaxt = 'n', xaxt = 'n', main = "", cex.main = 0.8)
axis(1, at = c(0, 30, 60), cex.axis = 1.1)
axis(2, at = c(0, 0.15), cex.axis = 1.)
dev.off()


quartz(file = "tissue_specificity_histogram.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
hist(expr_properties_df_with_canonical$tissue_specificity, col = alpha("red", 0.6), lty = 0, 
     breaks = 50, freq = FALSE, xlab = expression(paste("tissue specificity (",  tau, ")")), cex.lab = 1.35,
     yaxt = 'n', xaxt = 'n', main = "", cex.main = 0.8, xlim =  c(0, 1))
axis(1, at = c(0, 0.5, 1), cex.axis = 1.2)
axis(2, at = c(0, 1.5, 3), cex.axis = 1.2)
dev.off()

quartz(file = "polr2a_canonical_vs_max.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4, 5, 1, 2.2)+0.1)
plot(canonical_polr2a[which(expr_properties_df_with_canonical$tissue_specificity < 0.6 & 
                              canonical_polr2a < 10)], cex.lab = 1.1,
     max_polr2a[which(expr_properties_df_with_canonical$tissue_specificity < 0.6 & 
                        canonical_polr2a < 10)],
     cex = 0.5, col = alpha("red", 0.15), pch = 19, bty = 'l', yaxt = 'n', xaxt = 'n', 
     xlab = "# of ENCODE experiments\nwith POLR2A peak (canonical)", 
     ylab = "# of ENCODE experiments\nwith POLR2A peak (max)", 
     main = "")
axis(1, at = c(0, 9), cex.axis = 1)
axis(2, at = c(0, 30, 60), cex.axis = 1)
dev.off()

quartz(file = "polr2a_canonical_vs_tau.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4, 5.7, 2, 1)+0.1)
plot(expr_properties_df_with_canonical$tissue_specificity, canonical_polr2a, cex.lab = 1.35,
     cex = 0.5, col = alpha("red", 0.01), pch = 19, bty = 'l', yaxt = 'n', xaxt = 'n', xlim = c(0,1),
     xlab = expression(paste("tissue specificity (",  tau, ")")), ylab = "# of ENCODE experiments\nwith POLR2A peak", main = "")
axis(1, at = c(0, 0.5, 1), cex.axis = 1.2)
axis(2, at = c(0, 30, 60), cex.axis = 1.2)

dev.off()