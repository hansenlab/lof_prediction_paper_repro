###train the predictive model for LoF-intolerance model using the good power transcripts
full_df <- data.frame(lof_int = proms_good_power$is_intolerant, oe_cpg = proms_good_power$oe_cpg_ratio,
                      score = proms_good_power$score, exon_score = proms_good_power$exon_score)

sampled_indices <- sample(1:dim(full_df)[1], 3000)
train_set <- full_df[sampled_indices, ]
test_set <- full_df[-sampled_indices, ]
#save the indices and the train and test sets for easy reproducibility of the figures
save(sampled_indices, train_set, test_set, file = "prediction_sets.rda")

load(file = "prediction_sets.rda")
logistic_regression_fit <- glm(lof_int ~ oe_cpg + score + exon_score, data = train_set, 
                               family = "binomial")

predictions <- predict(logistic_regression_fit, newdata = test_set, type = "response")

###now make the validation plots

quartz(file = "loeuf_prediction_test_set.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
plot(density(na.omit(proms_good_power[-sampled_indices]$oe_lof_upper[which(predictions > 0.75)]), 
             from = 0, to = 2), col = "dark orange", lwd = 2.5, xlim = c(0, 2), cex.lab = 1.35,
     main = "test set (genes\nwith known constraint)", font.main = 1,
     cex.main = 1.5, bty = "l", xlab = "LOEUF", xaxt = 'n', yaxt = 'n')
lines(density(na.omit(proms_good_power[-sampled_indices]$oe_lof_upper[which(predictions < 0.25)]), 
              from = 0), col = alpha("cornflowerblue", 0.7), lwd = 2.5)
abline(v = 0.35, lty = "longdash", col = rgb(0,0,0,0.5))
legend <- legend("topright", legend = c("predicted label:", "highly LoF-intolerant",
                                        "non-highly LoF-intolerant"), lwd = 2.1, 
                 col = c("white", "dark orange", "cornflowerblue"), lty = "solid", bty = 'n', cex = 0.65)
axis(2, at = c(0, 1.5, 3), cex.axis = 1.2)
axis(1, at = c(0, 0.75, 1.5), cex.axis = 1.2)
dev.off()

##now validate using structural variants. First, load the structural variants
structural_variants_df <- read.delim('structural_variants_15k_genomes.bed', header = TRUE, stringsAsFactors = FALSE)
structural_variants_df <- structural_variants_df[which(structural_variants_df$SVTYPE == "DEL" & 
                                                         structural_variants_df$FILTER == "PASS"), ]

structural_variants_granges <- makeGRangesFromDataFrame(structural_variants_df, starts.in.df.are.0based = TRUE)
genome(seqinfo(structural_variants_granges)) <- "hg19"
seqlevelsStyle(structural_variants_granges) <- "ucsc"
structural_variants_granges$AC <- structural_variants_df[which(structural_variants_df$SVTYPE == "DEL" &
                                                                 structural_variants_df$FILTER == "PASS"), 30]
structural_variants_granges$AC <- as.numeric(structural_variants_granges$AC)


#only keep deletions overlapping exactly one promoter
all_overlaps <- findOverlaps(structural_variants_granges, proms_to_use)
structural_variants_granges <- structural_variants_granges[-unique(
  queryHits(all_overlaps)[which(duplicated(queryHits(all_overlaps)))])]



low_loeuf <- replicate(1000, getBootstrapProportionOfRangesOverlappingSVs(
  proms_good_power[which(proms_good_power$oe_lof_upper <= quantile(
    proms_good_power$oe_lof_upper, 0.25))]))
middle_loeuf_1 <- replicate(1000, getBootstrapProportionOfRangesOverlappingSVs(
  proms_good_power[which(proms_good_power$oe_lof_upper > 
                           quantile(proms_good_power$oe_lof_upper, 0.25) & 
                           proms_good_power$oe_lof_upper <= quantile(
                             proms_good_power$oe_lof_upper, 0.5))]))
middle_loeuf_2 <- replicate(1000, getBootstrapProportionOfRangesOverlappingSVs(
  proms_good_power[which(proms_good_power$oe_lof_upper > 
                           quantile(proms_good_power$oe_lof_upper, 0.5) & 
                           proms_good_power$oe_lof_upper <= quantile(
                             proms_good_power$oe_lof_upper, 0.75))]))
high_loeuf <- replicate(1000, getBootstrapProportionOfRangesOverlappingSVs(
  proms_good_power[which(proms_good_power$oe_lof_upper > quantile(
    proms_good_power$oe_lof_upper, 0.75))]))

quartz(file = "sv_validation_deletion_overlap.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4.5,5.5,2.9,1)+0.1)
boxplot(low_loeuf, frame = FALSE, lty = "solid", col = rgb(1,0,0,0.59), yaxt = 'n',
        xlim = c(0.8, 4.2), medlty = 0, ylim = c(0, 0.4), at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, xlab = "LOEUF quartile",
        ylab = "% of promoters with deletions\nin healthy individuals", cex.lab = 1.35)
points(1, getProportionOfRangesOverlappingSVs(proms_good_power[which(proms_good_power$oe_lof_upper <= 
                                                                       quantile(proms_good_power$oe_lof_upper, 0.25))]), pch = 19, cex= 0.45)
boxplot(middle_loeuf_1, frame = FALSE, lty = "solid", col = rgb(1,0,0,0.59), 
        medlty = 0, at = 2, yaxt = 'n',
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE)
points(2, getProportionOfRangesOverlappingSVs(proms_good_power[which(proms_good_power$oe_lof_upper > 
                                                                       quantile(proms_good_power$oe_lof_upper, 0.25) & 
                                                                       proms_good_power$oe_lof_upper <= 
                                                                       quantile(proms_good_power$oe_lof_upper, 0.5))]), pch = 19, cex= 0.45)
boxplot(middle_loeuf_2, frame = FALSE, lty = "solid", col = rgb(1,0,0,0.59), 
        medlty = 0, at = 3, yaxt = 'n',
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE)
points(3, getProportionOfRangesOverlappingSVs(proms_good_power[which(proms_good_power$oe_lof_upper > 
                                                                       quantile(proms_good_power$oe_lof_upper, 0.55) & 
                                                                       proms_good_power$oe_lof_upper <= 
                                                                       quantile(proms_good_power$oe_lof_upper, 0.75))]), pch = 19, cex= 0.45)
boxplot(high_loeuf, frame = FALSE, lty = "solid", col = rgb(1,0,0,0.59), 
        medlty = 0, at = 4, yaxt = 'n',
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE)
points(4, getProportionOfRangesOverlappingSVs(proms_good_power[which(proms_good_power$oe_lof_upper > 
                                                                       quantile(proms_good_power$oe_lof_upper, 0.75))]), pch = 19, cex= 0.45)
axis(1, at = c(1, 2, 3, 4), labels = c("1", "2", "3", "4"), 
     las = 2, cex.axis = 1.2, las = 1)
axis(2, at = c(0, 0.2, 0.40), labels = c("0", "20", "40"), cex.axis = 1.2)
dev.off()








####promoter deletion size vs LOEUF plots
quartz(file = "sv_validation_deletion_size.pdf", width = 4.4, height = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(2,2))
hist(getDistributionOfOverlappingSVSizes(proms_good_power[which(proms_good_power$oe_lof_upper < 
                                                                  quantile(proms_good_power$oe_lof_upper, 0.25))]), breaks = 20, freq = FALSE, 
     lty = 0, col = alpha("red", 0.6), xlab = "", xaxt = 'n', yaxt = 'n', main = "", cex.lab = 1.35)
legend <- legend("top", legend = c("1st LOEUF quartile"), col = c("white"), bty = 'n', cex = 1.5)
axis(2, at = c(0, 0.002), cex.axis = 1.2)

hist(getDistributionOfOverlappingSVSizes(proms_good_power[which(proms_good_power$oe_lof_upper > 
                                                                  quantile(proms_good_power$oe_lof_upper, 0.25) & 
                                                                  proms_good_power$oe_lof_upper < quantile(proms_good_power$oe_lof_upper, 0.5))]), breaks = 20, freq = FALSE, 
     lty = 0, col = alpha("red", 0.6), ylim = c(0, 0.002), xlab = "", xaxt = 'n', yaxt = 'n', 
     main = "", cex.lab = 1.35)
legend <- legend("top", legend = c("2nd LOEUF quartile"), col = c("white"), bty = 'n', cex = 1.5)
axis(2, at = c(0, 0.002), cex.axis = 1.2)

hist(getDistributionOfOverlappingSVSizes(proms_good_power[which(proms_good_power$oe_lof_upper > 
                                                                  quantile(proms_good_power$oe_lof_upper, 0.5) & 
                                                                  proms_good_power$oe_lof_upper < quantile(proms_good_power$oe_lof_upper, 0.75))]), breaks = 20, freq = FALSE, 
     lty = 0, col = alpha("red", 0.6), xlab = "deletion size", ylim = c(0, 0.002), xaxt = 'n', yaxt = 'n', 
     main = "", cex.lab = 1.35)
legend <- legend("top", legend = c("3rd LOEUF quartile"), col = c("white"), bty = 'n', cex = 1.5)
axis(2, at = c(0, 0.002), cex.axis = 1.2)

hist(getDistributionOfOverlappingSVSizes(proms_good_power[which(proms_good_power$oe_lof_upper > 
                                                                  quantile(proms_good_power$oe_lof_upper, 0.75))]), breaks = 20, freq = FALSE, 
     lty = 0, col = alpha("red", 0.6), xlab = "deletion size", ylim = c(0, 0.002), 
     xaxt = 'n', yaxt = 'n', 
     main = "", cex.lab = 1.35)
legend <- legend("top", legend = c("4th LOEUF quartile"), col = c("white"), bty = 'n', cex = 1.5)
axis(1, at = c(0, 2000, 4000), cex.axis = 1.2)
axis(2, at = c(0, 0.002), cex.axis = 1.2)

dev.off()











