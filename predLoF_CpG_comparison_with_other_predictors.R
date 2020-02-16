test_set_proms <- proms_good_power[as.numeric(rownames(test_set))]
score_list <- list(test_set_proms$haplo_score, test_set_proms$haplo_score_2, test_set_proms$haplo_score_3)
citations <- c("Huang et al.", "Steinberg et al.", "Han et al.")
colors <- c(alpha("cornflowerblue", 0.7), rgb(0,0,0,0.7), alpha("orange", 0.7))



quartz(file = "prediction_comparison.pdf", width = 3.6, height = 6.2, pointsize = 8, type = "pdf")
par(mfrow = c(3,2))
plot(0.5, col = rgb(1,1,1), xlim = c(0, 700), ylab = "precision (positive predictive value)", xaxt = 'n',
     yaxt = 'n', main = "", bty = 'l', ylim = c(0, 1), 
     xlab = "# genes correctly classified\nas highly LoF-intolerant", cex.lab = 1.35)
points(sapply(seq(0, 1, by = 0.05), function(xx) {
  table(test_set_proms$is_intolerant[which(
    score_list[[1]] >= xx)])[2]
}), sapply(seq(0, 1, by = 0.05), function(xx) {
  prop.table(table(test_set_proms$is_intolerant[which(
    score_list[[1]] >= xx)]))[2]
}), col = colors[1], pch = 19, cex = 1)

points(sapply(seq(0, 1, by = 0.05), function(xx) {
  table(test_set_proms$is_intolerant[which(
    predictions >= xx)])[2]
}), sapply(seq(0, 1, by = 0.05), function(xx) {
  prop.table(table(test_set_proms$is_intolerant[which(
    predictions >= xx)]))[2]
}), col = alpha("red", 0.7), pch = 19, cex = 1)
axis(2, at = c(0, 0.5, 1), cex.axis = 1.2)

legend <- legend("bottomleft", legend = c("CpG-density-based", citations[1]), 
                 col = c(alpha("red", 0.7), colors[1]), pch = 19, bty = 'n', cex = 1)
#
plot(0.5, col = rgb(1,1,1), xlim = c(0, 1290), ylab = "negative predictive value", xaxt = 'n',
     yaxt = 'n', main = "", bty = 'l', ylim = c(0, 1), cex.lab = 1.35,
     xlab = "# genes correctly classified\nas non-highly LoF-intolerant")
points(sapply(seq(0, 1, by = 0.05), function(xx) {
  table(test_set_proms$is_intolerant[which(
    score_list[[1]] <= xx)])[1]
}), sapply(seq(0, 1, by = 0.05), function(xx) {
  prop.table(table(test_set_proms$is_intolerant[which(
    score_list[[1]] <= xx)]))[1]
}), col = colors[1], pch = 19, cex = 1)

points(sapply(seq(0, 1, by = 0.05), function(xx) {
  table(test_set_proms$is_intolerant[which(
    predictions <= xx)])[1]
}), sapply(seq(0, 1, by = 0.05), function(xx) {
  prop.table(table(test_set_proms$is_intolerant[which(
    predictions <= xx)]))[1]
}), col = alpha("red", 0.7), pch = 19, cex = 1)
axis(2, at = c(0, 0.5, 1), cex.axis = 1.2)

#
sapply(2:length(score_list), function(xx) {
  score_to_compare <- score_list[[xx]]
  legend_name <- citations[[xx]]
  color <- colors[xx]
  
  plot(0.5, col = rgb(1,1,1), xlim = c(0, 700), ylab = "", xaxt = 'n',
       yaxt = 'n', main = "", bty = 'l', ylim = c(0, 1), 
       xlab = "", cex.lab = 1.35)
  points(sapply(seq(0, 1, by = 0.05), function(xx) {
    table(test_set_proms$is_intolerant[which(
      score_to_compare >= xx)])[2]
  }), sapply(seq(0, 1, by = 0.05), function(xx) {
    prop.table(table(test_set_proms$is_intolerant[which(
      score_to_compare >= xx)]))[2]
  }), col = color, pch = 19, cex = 1)
  axis(2, at = c(0, 0.5, 1), cex.axis = 1.2)
  axis(1, at = c(0, 300, 600), cex.axis = 1.2)
  
  points(sapply(seq(0, 1, by = 0.05), function(xx) {
    table(test_set_proms$is_intolerant[which(
      predictions >= xx)])[2]
  }), sapply(seq(0, 1, by = 0.05), function(xx) {
    prop.table(table(test_set_proms$is_intolerant[which(
      predictions >= xx)]))[2]
  }), col = alpha("red", 0.7), pch = 19, cex = 1)
  
  legend <- legend("bottomleft", legend = c("CpG-density-based", legend_name), 
                   col = c(alpha("red", 0.7), color), pch = 19, bty = 'n', cex = 1)
  
  #
  plot(0.5, col = rgb(1,1,1), xlim = c(0, 1290), ylab = "", xaxt = 'n',
       yaxt = 'n', main = "", bty = 'l', ylim = c(0, 1), 
       xlab = "")
  points(sapply(seq(0, 1, by = 0.05), function(xx) {
    table(test_set_proms$is_intolerant[which(
      score_to_compare <= xx)])[1]
  }), sapply(seq(0, 1, by = 0.05), function(xx) {
    prop.table(table(test_set_proms$is_intolerant[which(
      score_to_compare <= xx)]))[1]
  }), col = color, pch = 19, cex = 1)
  
  points(sapply(seq(0, 1, by = 0.05), function(xx) {
    table(test_set_proms$is_intolerant[which(
      predictions <= xx)])[1]
  }), sapply(seq(0, 1, by = 0.05), function(xx) {
    prop.table(table(test_set_proms$is_intolerant[which(
      predictions <= xx)]))[1]
  }), col = alpha("red", 0.7), pch = 19, cex = 1)
  axis(2, at = c(0, 0.5, 1), cex.axis = 1.2)
  axis(1, at = c(0, 600, 1200), cex.axis = 1.2)
}) 

dev.off()


lm_cpg_based_predictor <- summary(lm(test_set_proms$oe_lof_upper ~ predictions))
lm_haplo_score <- summary(lm(test_set_proms$oe_lof_upper ~ test_set_proms$haplo_score))
lm_haplo_score_2 <- summary(lm(test_set_proms$oe_lof_upper ~ test_set_proms$haplo_score_2))
lm_haplo_score_3 <- summary(lm(test_set_proms$oe_lof_upper ~ test_set_proms$haplo_score_3))

quartz(file = "prediction_variance_explained.pdf", height = 2.5, width = 1.8, pointsize = 8, type = "pdf")
par(mar = c(10, 7, 1, 1))
plot(1, 100*lm_haplo_score_2$adj.r.squared, pch = 19, cex = 1.2, col = alpha("red", 0.62), bty = 'l', 
     xlab = "", ylab = "% explained\nout-of-sample\nLOEUF variance", cex.lab = 1.25, ylim = c(0, 35), xlim = c(0.8, 4.2), xaxt = 'n', yaxt = 'n')
points(2, 100*lm_haplo_score_3$adj.r.squared, pch = 19, cex = 1.2, col = alpha("red", 0.62))
points(3, 100*lm_haplo_score$adj.r.squared, pch = 19, cex = 1.2, col = alpha("red", 0.62))
points(4, 100*lm_cpg_based_predictor$adj.r.squared, pch = 19, cex = 1.2, col = alpha("red", 0.62))
axis(2, at = c(0, 15, 30), cex.axis = 1.2)
axis(1, at = c(1,2,3,4), labels = c("Steinberg et al.", 
                                    "Han et al.",
                                    "Huang et al.",
                                    "CpG-density-based"), las = 2, cex.axis = 1.25)
abline(v = c(1.5, 2.5, 3.5, 4.5), lty = "longdash", col = rgb(0,0,0,0.7))
dev.off()

####

