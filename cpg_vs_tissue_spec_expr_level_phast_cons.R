###this script contains all plots that assess the relationship between promoter CpG density and tissue specificity/expression level/across species conservation


#####plot for tissue specificity
proms_granges_for_plot <- proms_good_power
proms_granges_for_plot$group <- "least tissue specific quartile"
proms_granges_for_plot$group[which(proms_granges_for_plot$tau < quantile(proms_granges_for_plot$tau, 0.5) 
                                   & proms_granges_for_plot$tau >= quantile(proms_granges_for_plot$tau, 0.25))] <- "2nd quartile"
proms_granges_for_plot$group[which(proms_granges_for_plot$tau < quantile(proms_granges_for_plot$tau, 0.75) 
                                   & proms_granges_for_plot$tau >= quantile(proms_granges_for_plot$tau, 0.5))] <- "3rd quartile"
proms_granges_for_plot$group[which(proms_granges_for_plot$tau >= quantile(proms_granges_for_plot$tau, 0.75))] <- "most tissue specific quartile"

proms_granges_for_plot$oe_cpg_category <- "top 25% promoter CpG density"
proms_granges_for_plot$oe_cpg_category[which(
  proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.75))] <- "in between"
proms_granges_for_plot$oe_cpg_category[which(
  proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.25))] <- "bottom 25% promoter CpG density"
proms_granges_for_plot$oe_cpg_category <- as.factor(proms_granges_for_plot$oe_cpg_category)

proms_granges_for_plot$group <- factor(proms_granges_for_plot$group, 
                                       levels = c("least tissue specific quartile", "2nd quartile", "3rd quartile", "most tissue specific quartile"), 
                                       labels = c("4", "3", "2", "1"))

proms_df_for_plot <- data.frame(constraint = proms_granges_for_plot$oe_lof_upper, 
                                group1 = proms_granges_for_plot$group, 
                                group2 = proms_granges_for_plot$oe_cpg_category)

quartz(file = "promoter_cpg_tissue_specificity.pdf", height = 2.5, width = 1.8, pointsize = 8, type = "pdf")
ggplot(proms_df_for_plot, aes(x = constraint, y = group1, color = group2, 
                              point_color = group2, fill = group2)) + 
  geom_density_ridges2(scale = 0.80) +
  scale_fill_manual(values = c(alpha("forest green", 0.57), alpha("yellow", 0.57), alpha("brown", 0.57)), 
                    labels = c("top 25% CpG density", "in between", "bottom 25% CpG density"), name = "") +  
  scale_color_manual(values = c(alpha("forest green", 0.57), alpha("yellow", 0.57), alpha("brown", 0.57)), guide = "none") +
  scale_discrete_manual("point_color", values = c(alpha("brown", 0.57), alpha("forest green", 0.57), alpha("royalblue1", 0.57)), guide = "none") + 
  xlim(0,2) + labs(x = "LOEUF", 
                   y = expression(paste("tissue spec (", tau, ")", " quartile"))) + 
  guides(fill = guide_legend(
    override.aes = list(
      fill = c(alpha("brown", 0.7), alpha("yellow", 0.7), alpha("forest green", 0.7)),
      color = NA, point_color = NA))) +
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size=6.2)) + theme(legend.position = c(0.35, 0.98)) + theme(axis.title=element_text(size=9)) +
  theme(legend.text = element_text(size = 3.7))
dev.off()
######



#####plot for expression level
proms_granges_for_plot <- proms_good_power
proms_granges_for_plot$group <- "least expression"
proms_granges_for_plot$group[which(proms_granges_for_plot$median_expr_in_maximum_expression_tissue < quantile(proms_granges_for_plot$median_expr_in_maximum_expression_tissue, 0.5) 
                                   & proms_granges_for_plot$median_expr_in_maximum_expression_tissue >= quantile(proms_granges_for_plot$median_expr_in_maximum_expression_tissue, 0.25))] <- "2nd quartile"
proms_granges_for_plot$group[which(proms_granges_for_plot$median_expr_in_maximum_expression_tissue < quantile(proms_granges_for_plot$median_expr_in_maximum_expression_tissue, 0.75) 
                                   & proms_granges_for_plot$median_expr_in_maximum_expression_tissue >= quantile(proms_granges_for_plot$median_expr_in_maximum_expression_tissue, 0.5))] <- "3rd quartile"
proms_granges_for_plot$group[which(proms_granges_for_plot$median_expr_in_maximum_expression_tissue >= quantile(proms_granges_for_plot$median_expr_in_maximum_expression_tissue, 0.75))] <- "greatest expression"

proms_granges_for_plot$oe_cpg_category <- "top 25% promoter CpG density"
proms_granges_for_plot$oe_cpg_category[which(
  proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.75))] <- "in between"
proms_granges_for_plot$oe_cpg_category[which(
  proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.25))] <- "bottom 25% promoter CpG density"
proms_granges_for_plot$oe_cpg_category <- as.factor(proms_granges_for_plot$oe_cpg_category)


proms_granges_for_plot$group <- factor(proms_granges_for_plot$group, 
                                       levels = c("least expression", "2nd quartile", "3rd quartile", "greatest expression"), 
                                       labels = c("4", "3", "2", "1"))


proms_df_for_plot <- data.frame(constraint = proms_granges_for_plot$oe_lof_upper, 
                                group1 = proms_granges_for_plot$group, 
                                group2 = proms_granges_for_plot$oe_cpg_category)


quartz(file = "promoter_cpg_expression_level.pdf", height = 2.5, width = 1.8, pointsize = 8, type = "pdf")
ggplot(proms_df_for_plot, aes(x = constraint, y = group1, color = group2, 
                              point_color = group2, fill = group2)) + 
  geom_density_ridges2(scale = 0.57) +
  scale_fill_manual(values = c(alpha("forest green", 0.57), alpha("yellow", 0.57), alpha("brown", 0.57)), 
                    labels = c("", "", ""), name = "") +  
  scale_color_manual(values = c(alpha("forest green", 0.57), alpha("yellow", 0.57), alpha("brown", 0.57)), guide = "none") +
  scale_discrete_manual("point_color", values = c(alpha("brown", 0.57), alpha("forest green", 0.57), alpha("royalblue1", 0.57)), guide = "none") + 
  xlim(0,2) + labs(x = "LOEUF", y = "expr level quartile") + 
  guides(fill = FALSE) +
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size=6.2)) + theme(legend.position = c(0.35, 0.98)) + theme(axis.title=element_text(size=9)) +
  theme(legend.text = element_text(size = 3.7))
dev.off()
######

####variance explained
lm_cpg_only <- summary(lm(proms_good_power$oe_lof_upper ~ proms_good_power$oe_cpg_ratio))
lm_tau_only <- summary(lm(proms_good_power$oe_lof_upper ~ proms_good_power$tau))
lm_expr_level_only <- summary(lm(proms_good_power$oe_lof_upper ~ proms_good_power$median_expr_in_maximum_expression_tissue))
lm_expr_level_and_tau <- summary(lm(proms_good_power$oe_lof_upper ~ proms_good_power$median_expr_in_maximum_expression_tissue + 
                                      proms_good_power$tau))

quartz(file = "promoter_cpg_variance_explained.pdf", height = 2.5, width = 1.8, pointsize = 8, type = "pdf")
par(mar = c(10, 5.5, 1, 1))
plot(1, 100*lm_expr_level_only$adj.r.squared, pch = 19, cex = 1.2, col = alpha("red", 0.62), bty = 'l', 
     xlab = "", ylab = "% explained\nLOEUF variance", cex.lab = 1.25, ylim = c(0, 21), xlim = c(0.8, 4.2), xaxt = 'n', yaxt = 'n')
points(2, 100*lm_tau_only$adj.r.squared, pch = 19, cex = 1.2, col = alpha("red", 0.62))
points(3, 100*lm_expr_level_and_tau$adj.r.squared, pch = 19, cex = 1.2, col = alpha("red", 0.62))
points(4, 100*lm_cpg_only$adj.r.squared, pch = 19, cex = 1.2, col = alpha("red", 0.62))
axis(2, at = c(0, 10, 20), cex.axis = 1.2)
axis(1, at = c(1,2,3,4), labels = c("expr level", expression(tau),
                                    expression(paste(tau, " x expr level")), 
                                    "prom CpG density"), las = 2, cex.axis = 1.25)
#axis(1, at = c(1,2,3,4), labels = c("expr level", "tissue spec",
#                                      "tissue spec x expr level", 
#                                    "prom CpG density"), las = 2, cex.axis = 1.35)
abline(v = c(1.5, 2.5, 3.5), lty = "longdash", col = rgb(0,0,0,0.7))
dev.off()




######supplemental plots for the relationship between tissue specificity/expression level with promoter CpG density
proms_granges_for_plot <- proms_good_power
proms_granges_for_plot$group <- 10
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.9) 
                                   & proms_granges_for_plot$oe_cpg_ratio > quantile(proms_granges_for_plot$oe_cpg_ratio, 0.8))] <- 9
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.8) 
                                   & proms_granges_for_plot$oe_cpg_ratio > quantile(proms_granges_for_plot$oe_cpg_ratio, 0.7))] <- 8
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.7) 
                                   & proms_granges_for_plot$oe_cpg_ratio > quantile(proms_granges_for_plot$oe_cpg_ratio, 0.6))] <- 7
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.6) 
                                   & proms_granges_for_plot$oe_cpg_ratio > quantile(proms_granges_for_plot$oe_cpg_ratio, 0.5))] <- 6
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.5) 
                                   & proms_granges_for_plot$oe_cpg_ratio > quantile(proms_granges_for_plot$oe_cpg_ratio, 0.4))] <- 5
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.4) 
                                   & proms_granges_for_plot$oe_cpg_ratio > quantile(proms_granges_for_plot$oe_cpg_ratio, 0.3))] <- 4
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.3) 
                                   & proms_granges_for_plot$oe_cpg_ratio > quantile(proms_granges_for_plot$oe_cpg_ratio, 0.2))] <- 3
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.2) 
                                   & proms_granges_for_plot$oe_cpg_ratio > quantile(proms_granges_for_plot$oe_cpg_ratio, 0.1))] <- 2
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.1) 
)] <- 1

proms_granges_for_plot$group <- factor(proms_granges_for_plot$group, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 
                                       labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))



proms_df_for_plot <- data.frame(group1 = proms_granges_for_plot$group, 
                                tissue_spec = proms_granges_for_plot$n_tissues_with_detectable_expression)

quartz(file = "promoter_cpg_tissue_specificity_deciles_1.pdf", height = 2.5, width = 1.8, pointsize = 8, type = "pdf")
ggplot(proms_df_for_plot, aes(x = tissue_spec, y = group1)) + 
  geom_density_ridges2(fill = alpha("red", 0.57), color = rgb(0,0,0,0.7), from = 0, to = 53) + 
  labs(x = "# tissues with\ndetectable expression", y = "") + 
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size=6.2)) + theme(axis.title=element_text(size=9))
dev.off()


proms_df_for_plot <- data.frame(group1 = proms_granges_for_plot$group, 
                                tissue_spec = proms_granges_for_plot$tau)

quartz(file = "promoter_cpg_tissue_specificity_deciles_2.pdf", height = 2.5, width = 1.8, pointsize = 8, type = "pdf")
ggplot(proms_df_for_plot, aes(x = tissue_spec, y = group1)) + 
  geom_density_ridges2(fill = alpha("red", 0.57), color = rgb(0,0,0,0.7), from = 0, to = 1) + 
  labs(x = expression(paste("tissue spec (", tau, ")")), y = "") + 
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size=6.2)) + theme(axis.title=element_text(size=9))
dev.off()


proms_df_for_plot <- data.frame(group1 = proms_granges_for_plot$group, 
                                tissue_spec = proms_granges_for_plot$median_expr_in_maximum_expression_tissue)

quartz(file = "promoter_expression_level_deciles.pdf", height = 2.5, width = 1.8, pointsize = 8, type = "pdf")
ggplot(proms_df_for_plot, aes(x = tissue_spec, y = group1)) + 
  geom_density_ridges2(fill = alpha("red", 0.57), color = rgb(0,0,0,0.7), from = 0, to = 15) + 
  labs(x = "expr level (log2(TPM+1))", y = "CpG density decile") + 
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size=6.2)) + theme(axis.title=element_text(size=9))
dev.off()





###########relationship between promoter CpG density and across-species conservation
#####phast cons
proms_granges_for_plot <- proms_good_power[-which(is.na(proms_good_power$exon_score))]
proms_granges_for_plot$group <- "least constrained quartile"
proms_granges_for_plot$group[which(proms_granges_for_plot$exon_score < quantile(proms_granges_for_plot$exon_score, 0.5) 
                                   & proms_granges_for_plot$exon_score >= quantile(proms_granges_for_plot$exon_score, 0.25))] <- "2nd quartile"
proms_granges_for_plot$group[which(proms_granges_for_plot$exon_score < quantile(proms_granges_for_plot$exon_score, 0.75) 
                                   & proms_granges_for_plot$exon_score >= quantile(proms_granges_for_plot$exon_score, 0.5))] <- "3rd quartile"
proms_granges_for_plot$group[which(proms_granges_for_plot$exon_score >= quantile(proms_granges_for_plot$exon_score, 0.75))] <- "most constrained quartile"

proms_granges_for_plot$oe_cpg_category <- "top 25% promoter CpG density"
proms_granges_for_plot$oe_cpg_category[which(
  proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.75))] <- "in between"
proms_granges_for_plot$oe_cpg_category[which(
  proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.25))] <- "bottom 25% promoter CpG density"
proms_granges_for_plot$oe_cpg_category <- as.factor(proms_granges_for_plot$oe_cpg_category)

proms_granges_for_plot$group <- factor(proms_granges_for_plot$group, 
                                       levels = c("least constrained quartile", "2nd quartile", "3rd quartile", "most constrained quartile"), 
                                       labels = c("1", "2", "3", "4"))

proms_df_for_plot <- data.frame(constraint = proms_granges_for_plot$oe_lof_upper, 
                                group1 = proms_granges_for_plot$group, 
                                group2 = proms_granges_for_plot$oe_cpg_category)

quartz(file = "promoter_cpg_vs_exon_phastcons.pdf", height = 2.5, width = 1.8, pointsize = 8, type = "pdf")
ggplot(proms_df_for_plot, aes(x = constraint, y = group1, color = group2, 
                              point_color = group2, fill = group2)) + 
  geom_density_ridges2(scale = 0.95) +
  scale_fill_manual(values = c(alpha("forest green", 0.57), alpha("yellow", 0.57), alpha("brown", 0.57)), 
                    labels = c("", "", ""), name = "") +  
  scale_color_manual(values = c(alpha("forest green", 0.57), alpha("yellow", 0.57), alpha("brown", 0.57)), guide = "none") +
  scale_discrete_manual("point_color", values = c(alpha("brown", 0.57), alpha("forest green", 0.57), alpha("royalblue1", 0.57)), guide = "none") + 
  xlim(0,2) + labs(x = "LOEUF", y = "Exon PhastCons quartile") + 
  guides(fill = FALSE) +
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size=6.2)) + theme(legend.position = c(0.35, 0.98)) + theme(axis.title=element_text(size=9)) +
  theme(legend.text = element_text(size = 3.7))
dev.off()
######

proms_granges_for_plot <- proms_good_power[-which(is.na(proms_good_power$score))]
proms_granges_for_plot$group <- "least constrained quartile"
proms_granges_for_plot$group[which(proms_granges_for_plot$score < quantile(proms_granges_for_plot$score, 0.5) 
                                   & proms_granges_for_plot$score >= quantile(proms_granges_for_plot$score, 0.25))] <- "2nd quartile"
proms_granges_for_plot$group[which(proms_granges_for_plot$score < quantile(proms_granges_for_plot$score, 0.75) 
                                   & proms_granges_for_plot$score >= quantile(proms_granges_for_plot$score, 0.5))] <- "3rd quartile"
proms_granges_for_plot$group[which(proms_granges_for_plot$score >= quantile(proms_granges_for_plot$score, 0.75))] <- "most constrained quartile"

proms_granges_for_plot$oe_cpg_category <- "top 25% promoter CpG density"
proms_granges_for_plot$oe_cpg_category[which(
  proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.75))] <- "in between"
proms_granges_for_plot$oe_cpg_category[which(
  proms_granges_for_plot$oe_cpg_ratio < quantile(proms_granges_for_plot$oe_cpg_ratio, 0.25))] <- "bottom 25% promoter CpG density"
proms_granges_for_plot$oe_cpg_category <- as.factor(proms_granges_for_plot$oe_cpg_category)

proms_granges_for_plot$group <- factor(proms_granges_for_plot$group, 
                                       levels = c("least constrained quartile", "2nd quartile", "3rd quartile", "most constrained quartile"), 
                                       labels = c("1", "2", "3", "4"))

proms_df_for_plot <- data.frame(constraint = proms_granges_for_plot$oe_lof_upper, 
                                group1 = proms_granges_for_plot$group, 
                                group2 = proms_granges_for_plot$oe_cpg_category)

quartz(file = "promoter_cpg_vs_phastcons.pdf", height = 2.5, width = 1.8, pointsize = 8, type = "pdf")
ggplot(proms_df_for_plot, aes(x = constraint, y = group1, color = group2, 
                              point_color = group2, fill = group2)) + 
  geom_density_ridges2(scale = 0.95) +
  scale_fill_manual(values = c(alpha("forest green", 0.57), alpha("yellow", 0.57), alpha("brown", 0.57)), 
                    labels = c("top 25% CpG density", "in between", "bottom 25% CpG density"), name = "") +  
  scale_color_manual(values = c(alpha("forest green", 0.57), alpha("yellow", 0.57), alpha("brown", 0.57)), guide = "none") +
  scale_discrete_manual("point_color", values = c(alpha("brown", 0.57), alpha("forest green", 0.57), alpha("royalblue1", 0.57)), guide = "none") + 
  xlim(0,2) + labs(x = "LOEUF", y = "Promoter PhastCons quartile") + 
  guides(fill = guide_legend(
    override.aes = list(
      fill = c(alpha("brown", 0.7), alpha("yellow", 0.7), alpha("forest green", 0.7)),
      color = NA, point_color = NA))) +
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size=6.2)) + theme(legend.position = c(0.35, 0.98)) + theme(axis.title=element_text(size=9)) +
  theme(legend.text = element_text(size = 3.7))
dev.off()
######





####variance explained
lm_cpg_only <- summary(lm(proms_good_power$oe_lof_upper ~ proms_good_power$oe_cpg_ratio))
lm_promoter_only <- summary(lm(proms_good_power$oe_lof_upper ~ proms_good_power$score))
lm_exon_only <- summary(lm(proms_good_power$oe_lof_upper ~ proms_good_power$exon_score))
lm_promoter_and_exon <- summary(lm(proms_good_power$oe_lof_upper ~ proms_good_power$exon_score + proms_good_power$score))
lm_cpg_and_promoter_and_exon <- summary(lm(proms_good_power$oe_lof_upper ~ proms_good_power$oe_cpg_ratio + 
                                             proms_good_power$score + 
                                             proms_good_power$exon_score))

#try dotplot instead
quartz(file = "promoter_cpg_variance_explained_2.pdf", height = 2.5, width = 1.8, pointsize = 8, type = "pdf")
par(mar = c(12, 5.5, 1, 1))
plot(1, 100*lm_exon_only$adj.r.squared, pch = 19, cex = 1.2, col = alpha("red", 0.62), bty = 'l', 
     xlab = "", ylab = "% explained\nLOEUF variance", cex.lab = 1.25, ylim = c(0, 35), xlim = c(0.8, 5.2), xaxt = 'n', yaxt = 'n')
points(2, 100*lm_promoter_only$adj.r.squared, pch = 19, cex = 1.2, col = alpha("red", 0.62))
points(3, 100*lm_cpg_only$adj.r.squared, pch = 19, cex = 1.2, col = alpha("red", 0.62))
points(4, 100*lm_promoter_and_exon$adj.r.squared, pch = 19, cex = 1.2, col = alpha("red", 0.62))
points(5, 100*lm_cpg_and_promoter_and_exon$adj.r.squared, pch = 19, cex = 1.2, col = alpha("red", 0.62))
axis(2, at = c(0, 15, 30), cex.axis = 1.2)
axis(1, at = c(1,2,3,4,5), labels = c("exon PhastCons", "prom PhastCons",
                                      "prom CpG density", "prom+exon PhastCons", 
                                      "all three combined"), las = 2, cex.axis = 1.25)
abline(v = c(1.5, 2.5, 3.5, 4.5), lty = "longdash", col = rgb(0,0,0,0.7))
dev.off()

