###plots on EZH2 promoter binding and the relationship with constraint
quartz(file = "ezh2_binding.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4,5.7,4,1.5)+0.1)
plot(0.5, col = rgb(1,1,1), xlim = c(0.7, 8.4), 
     ylab = "Median # ENCODE experiments\nwith EZH2 peak", xaxt = 'n',
     yaxt = 'n', main = "", bty = 'l', ylim = c(0, 6), xlab = "CpG density quartile", 
     cex.lab = 1.35, cex = 1.2)
points(seq(1, 8, by = 2), 
       sapply(5:2, 
              function(xx) median(proms_to_use$ezh2_peak_count[which(
                proms_to_use$oe_cpg_ratio >= quantile(proms_to_use$oe_cpg_ratio)[xx-1] & 
                  proms_to_use$oe_cpg_ratio < quantile(proms_to_use$oe_cpg_ratio)[xx] & 
                  proms_to_use$tau < 0.6)])), pch = 19, col = rgb(0,0,0,0.7))
points(seq(2, 8, by = 2), 
       sapply(5:2, 
              function(xx) median(proms_to_use$ezh2_peak_count[which(
                proms_to_use$oe_cpg_ratio >= quantile(proms_to_use$oe_cpg_ratio)[xx-1] & 
                  proms_to_use$oe_cpg_ratio < quantile(proms_to_use$oe_cpg_ratio)[xx] & 
                  proms_to_use$tau > 0.6)])), pch = 19, col = "hotpink")
abline(v = c(2.5, 4.5, 6.5), lty = "longdash", col = rgb(0,0,0,0.7))
axis(2, at = c(0, 2, 4, 6), cex.axis = 1.2)
axis(1, at = c(1.5, 3.5, 5.5, 7.8), labels = c("1", "2", "3", "4"), cex.axis = 1.2)
dev.off()

###
quartz(file = "loeuf_ezh2.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
plot(density(proms_good_power$oe_lof_upper[
  which(proms_good_power$ezh2_peak_count >= 2 & 
          proms_good_power$tau > 0.6 & proms_good_power$oe_cpg_ratio > 
          quantile(proms_to_use$oe_cpg_ratio, 0.75))], from = 0, to = 2), 
  lwd = 2.5, bty = 'l', col = "hotpink", ylab = "Density", xlab = "LOEUF",
  xaxt = 'n', yaxt = 'n', main = "tissue-specific genes w/\nhigh-CpG-density promoters", 
  font.main = 1, xlim = c(0, 2), cex.main = 1.37, cex.lab = 1.35)
lines(density(proms_good_power$oe_lof_upper[
  which(proms_good_power$ezh2_peak_count < 2 & 
          proms_good_power$tau > 0.6 & proms_good_power$oe_cpg_ratio > 
          quantile(proms_to_use$oe_cpg_ratio, 0.75))], from = 0, to = 2), 
  lwd = 2.5, col = rgb(0,0,0,0.7))
axis(1, at = c(0, 0.75, 1.5), cex.axis = 1.2)
axis(2, at = c(0, 1.5), cex.axis = 1.2)
legend("topright", legend = c("with ezh2-bound\npromoter", 
                              "without ezh2-bound\npromoter"), 
       bty = 'n', lwd = 2.1, cex = 0.75, col = c("hotpink", rgb(0,0,0,0.7)))
dev.off()

summary(lm(proms_good_power$oe_lof_upper[which(proms_good_power$tau > 0.6 & proms_good_power$oe_cpg_ratio > 
             quantile(proms_to_use$oe_cpg_ratio, 0.75))] ~ proms_good_power$oe_cpg_ratio[which(proms_good_power$tau > 0.6 & 
                                                           proms_good_power$oe_cpg_ratio > quantile(proms_to_use$oe_cpg_ratio, 0.75))]*
             as.factor(proms_good_power$ezh2_peak_count[which(proms_good_power$tau > 0.6 & 
                             proms_good_power$oe_cpg_ratio > quantile(proms_to_use$oe_cpg_ratio, 0.75))]>=2)))

