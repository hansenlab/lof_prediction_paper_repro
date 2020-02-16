###
low_power_set <- data.frame(lof_int = proms_low_power$is_intolerant, 
                            oe_cpg = proms_low_power$oe_cpg_ratio, 
                            score = proms_low_power$score,
                            exon_score = proms_low_power$exon_score)


predictions <- predict(logistic_regression_fit, newdata = low_power_set, type = "response")

quartz(file = "oe_lof_small_protein_coding.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
plot(density(na.omit(proms_low_power$oe_lof[
  which(predictions > 0.75)]), from = 0), col = "dark orange", lwd = 2.5, 
  main = "", cex.main = 1.5, cex.lab = 1.35,
  bty = "l", font.main = 1,
  xlab = "Obs/Exp LoF point estimate", xaxt = 'n', yaxt = 'n')
lines(density(na.omit(proms_low_power$oe_lof[
  which(predictions < 0.25)]), from = 0), col = "cornflowerblue", lwd = 2.5)
axis(1, at = c(0, 0.75, 1.5), cex.axis = 1.2)
axis(2, at = c(0, 1.5), cex.axis = 1.2)

dev.off()

quartz(file = "sv_deletion_size_small_protein_coding_1.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
hist(getDistributionOfOverlappingSVSizes(proms_low_power[
  which(predictions > 0.75)]), breaks = 20, freq = FALSE, ylab = "",
  lty = 0, col = "dark orange", main = "", xlab = "", xaxt= 'n', 
  yaxt = 'n', ylim = c(0, 0.003), cex.main = 1.5, cex.lab = 1.35)
axis(2, at = c(0, 0.002), labels = c("0", ""), cex.axis = 1.2)
#legend <- legend("top", legend = c(""), col = c("white"), bty = 'n')
dev.off()

quartz(file = "sv_deletion_size_small_protein_coding_2.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
hist(getDistributionOfOverlappingSVSizes(proms_low_power[
  which(predictions < 0.25)]), breaks = 20, freq = FALSE, ylab = "",
  lty = 0, col = alpha("cornflowerblue"), main = "", xlab = "deletion size", xaxt= 'n', 
  yaxt = 'n', ylim = c(0, 0.003), cex.main = 1.5, cex.lab = 1.35)
#legend <- legend("top", legend = c("predicted\nhaplosufficient"), col = c("white"), bty = 'n')
axis(2, at = c(0, 0.002), labels = c("0", ""), cex.axis = 1.2)
axis(1, at = c(0, 2000, 4000), cex.axis = 1.2)

dev.off()


quartz(file = "sv_overlap_small_protein_coding.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
par(mar = c(1,5.7,2.5,1.5)+0.1)
barplot(c(100*getProportionOfRangesOverlappingSVs(proms_low_power[
  which(predictions > 0.75)]),
  100*getProportionOfRangesOverlappingSVs(proms_low_power[
    which(predictions < 0.25)])), 
  ylab = "% of promoters with deletions\nin healthy individuals", cex.lab = 1.35,
  names.arg = c("", ""),  
  border = "white", las = 2, ylim = c(0, 40),
  col = c("dark orange", "cornflowerblue"), 
  main = "", cex.main = 1.5, font.main = 1,
  yaxt = 'n', xlim = c(0, 3))
axis(2, at = c(0, 20, 40), cex.axis = 1.2)
dev.off()

permutation_dist <- replicate(10000, { # assess the significance of the difference using a permutation test
  getProportionOfRangesOverlappingSVs(proms_low_power[which(predictions > 0.75 | predictions < 0.25)][
    sample(1:length(proms_low_power[which(predictions > 0.75 | predictions < 0.25)]), length(which(predictions > 0.75)))]) - 
    getProportionOfRangesOverlappingSVs(proms_low_power[which(predictions > 0.75 | predictions < 0.25)][
      sample(1:length(proms_low_power[which(predictions > 0.75 | predictions < 0.25)]), length(which(predictions < 0.25)))])
})



####check for enrichment of TF genes in the genes predicted as highly LoF-intolerant
positive_prediction <- proms_low_power[which(predictions >= 0.75)]$gene_id
negative_prediction <- proms_low_power[which(predictions <= 0.75)]$gene_id

tf_genes <- read_csv('tf_gene_list.csv') #table from Barrera et al., 2016


n11 <- length(which(positive_prediction %in% tf_genes$`Enseml Gene ID`))
n12 <- length(positive_prediction) - n11
n21 <- length(which(negative_prediction %in% tf_genes$`Enseml Gene ID`))
n22 <- length(negative_prediction) - n21
fisher.test(matrix(c(n11, n12, n21, n22), nrow = 2))



