quartz(file = "transcripts_lof_point_esimate_vs_upper_ci_bound_all.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(constraint$oe_lof[which(constraint$exp_lof >= 20 & constraint$oe_lof <= 2)], 
     constraint$oe_lof_upper[which(constraint$exp_lof >= 20 & constraint$oe_lof <= 2)], 
     pch = 19, col = alpha("red", 0.01), xlab = "Obs/Exp LoF point estimate", ylab = "LOEUF", cex.lab = 1.35,
     main = "", xlim = c(0, 2), bty = 'l', cex = 0.5, xaxt = 'n', yaxt = 'n')
points(constraint$oe_lof[which(constraint$exp_lof < 20 & constraint$exp_lof > 10 & constraint$oe_lof <= 2)], 
       constraint$oe_lof_upper[which(constraint$exp_lof < 20 & constraint$exp_lof > 10 & constraint$oe_lof <= 2)], 
       pch = 19, col = alpha("yellow", 0.01), cex = 0.5)
points(constraint$oe_lof[which(constraint$exp_lof <= 10 & constraint$oe_lof <= 2)], 
       constraint$oe_lof_upper[which(constraint$exp_lof <= 10 & constraint$oe_lof <= 2)], 
       pch = 19, col = alpha("royalblue1", 0.01), cex = 0.5)
legend <- legend("bottomright", legend = c("genes with >= 20\nexp LoF variants",
                                           "genes in between",
                                           "genes with <= 10\nexp LoF variants"), 
                 pch = 19, bty = 'n', col = c(alpha("red", 0.8), "yellow", "royalblue1"), cex = 0.75)
abline(h = 0.35, lty = "longdash", col = rgb(0,0,0,0.7))
axis(1, at = c(0, 1, 2), cex.axis = 1.2)
axis(2, at = c(0, 1, 2), cex.axis = 1.2)
dev.off()



quartz(file = "transcripts_lof_point_esimate_vs_upper_ci_bound_canonical.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(constraint$oe_lof[which(constraint$exp_lof >= 20 & constraint$canonical == TRUE & constraint$oe_lof <= 2)], 
     constraint$oe_lof_upper[which(constraint$exp_lof >= 20 & constraint$canonical == TRUE & constraint$oe_lof <= 2)], 
     pch = 19, col = alpha("red", 0.01), xlab = "Obs/Exp LoF point estimate", ylab = "LOEUF", cex.lab = 1.35,
     main = "", xlim = c(0, 2), bty = 'l', cex = 0.75, xaxt = 'n', yaxt = 'n')
points(constraint$oe_lof[which(constraint$exp_lof < 20 & constraint$exp_lof > 10 & constraint$canonical == TRUE & constraint$oe_lof <= 2)], 
       constraint$oe_lof_upper[which(constraint$exp_lof < 20 & constraint$exp_lof > 10 & constraint$canonical == TRUE & constraint$oe_lof <= 2)], 
       pch = 19, col = alpha("yellow", 0.01), cex = 0.75)
points(constraint$oe_lof[which(constraint$exp_lof <= 10 & constraint$canonical == TRUE & constraint$oe_lof <= 2)], 
       constraint$oe_lof_upper[which(constraint$exp_lof <= 10 & constraint$canonical == TRUE & constraint$oe_lof <= 2)], 
       pch = 19, col = alpha("royalblue1", 0.01), cex = 0.75)
legend <- legend("bottomright", legend = c("genes with >= 20\nexp lof variants",
                                           "genes in between",
                                           "genes with <= 10\nexp lof variants"), 
                 pch = 19, bty = 'n', col = c(alpha("red", 0.8), "royalblue1"), cex = 0.75)
abline(h = 0.35, lty = "longdash", col = rgb(0,0,0,0.7))
axis(1, at = c(0, 1, 2), cex.axis = 1.2)
axis(2, at = c(0, 1, 2), cex.axis = 1.2)
dev.off()

