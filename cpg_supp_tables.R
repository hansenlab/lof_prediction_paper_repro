###make the two supplementary tables

###first, the table with predictions for unascertained genes

##predicted highly LoF-intolerant genes
proms_predicted_high <- proms_low_power[which(predictions > 0.75)]

edb <- EnsDb.Hsapiens.v75
proms_all <- promoters(edb, upstream = 2000, downstream = 2000, 
                       columns = c("gene_name", "tx_id", "tx_cds_seq_start", "tx_cds_seq_end",
                                   "tx_biotype", "gene_id"))
genome(seqinfo(proms_all)) <- "hg19"
seqlevelsStyle(proms_all) <- "ucsc"

proms_all_prot_coding <- proms_all[which(proms_all$tx_biotype == "protein_coding")]
#proms_all_prot_coding <- proms_all_prot_coding[-which(proms_all_prot_coding$tx_id %in% proms_predicted_high$tx_id)]

overlaps <- findOverlaps(proms_predicted_high, proms_all_prot_coding)
split_overlaps_by_gene <- split(overlaps, queryHits(overlaps))
overlapping_promoters_gene_ids <- lapply(split_overlaps_by_gene, function(xx) unique(proms_all_prot_coding$gene_id[subjectHits(xx)]))
names(overlapping_promoters_gene_ids) <- proms_predicted_high$gene_id #I use gene names here because there's no duplicates
number_of_overlapping_genes <- lengths(overlapping_promoters_gene_ids)
table(number_of_overlapping_genes) #check how many have more than 1. All genes have at least 1 overlapping gene, because they overlap with themselves

which(number_of_overlapping_genes > 1)
overlapping_promoters_gene_ids[which(number_of_overlapping_genes > 1)] <- sapply(overlapping_promoters_gene_ids[which(number_of_overlapping_genes > 1)], 
                                                                                 function(xx) paste0(xx[1], ", ", xx[2]))

#create supp table data frame
supp_table_data_frame_high_pred <- data.frame(gene_name = proms_predicted_high$gene_name, gene_id = proms_predicted_high$gene_id, 
                                              transcript_id = proms_predicted_high$tx_id,
                                              prediction_probability_of_LoF_intolerance_by_predLoF_CpG = predictions[which(predictions > 0.75)], 
                                              chr = seqnames(proms_predicted_high),
                                              promoter_start_coord = start(proms_predicted_high), promoter_end_coord = end(proms_predicted_high),
                                              other_genes_with_overlapping_promoter = "none")
supp_table_data_frame_high_pred$other_genes_with_overlapping_promoter <- as.character(supp_table_data_frame_high_pred$other_genes_with_overlapping_promoter)

supp_table_data_frame_high_pred$other_genes_with_overlapping_promoter[
  which(number_of_overlapping_genes > 1)] <- sapply(overlapping_promoters_gene_ids[which(number_of_overlapping_genes > 1)], 
                                                    function(xx) as.character(xx))


##
proms_predicted_low <- proms_low_power[which(predictions < 0.25)]

overlaps <- findOverlaps(proms_predicted_low, proms_all_prot_coding)
split_overlaps_by_gene <- split(overlaps, queryHits(overlaps))
overlapping_promoters_gene_ids <- lapply(split_overlaps_by_gene, function(xx) unique(proms_all_prot_coding$gene_id[subjectHits(xx)]))
names(overlapping_promoters_gene_ids) <- proms_predicted_low$gene_id #I use gene names here because there's no duplicates
number_of_overlapping_genes <- lengths(overlapping_promoters_gene_ids)
table(number_of_overlapping_genes) #check how many have more than 1. All genes have at least 1 overlapping gene, because they overlap with themselves

which(number_of_overlapping_genes > 1)
overlapping_promoters_gene_ids[which(number_of_overlapping_genes > 1)] <- sapply(overlapping_promoters_gene_ids[which(number_of_overlapping_genes > 1)], 
                                                                                 function(xx) paste0(xx[1], ", ", xx[2]))

#create supp table data frame
supp_table_data_frame_low_pred <- data.frame(gene_name = proms_predicted_low$gene_name, gene_id = proms_predicted_low$gene_id, 
                                             transcript_id = proms_predicted_low$tx_id,
                                             prediction_probability_of_LoF_intolerance_by_predLoF_CpG = predictions[which(predictions < 0.25)], 
                                             chr = seqnames(proms_predicted_low),
                                             promoter_start_coord = start(proms_predicted_low), promoter_end_coord = end(proms_predicted_low),
                                             other_genes_with_overlapping_promoter = "none")
supp_table_data_frame_low_pred$other_genes_with_overlapping_promoter <- as.character(supp_table_data_frame_low_pred$other_genes_with_overlapping_promoter)

supp_table_data_frame_low_pred$other_genes_with_overlapping_promoter[
  which(number_of_overlapping_genes > 1)] <- sapply(overlapping_promoters_gene_ids[which(number_of_overlapping_genes > 1)], 
                                                    function(xx) as.character(xx))

supp_table_all_predictions <- rbind(supp_table_data_frame_high_pred, supp_table_data_frame_low_pred)
colnames(supp_table_all_predictions)[4] <- "prediction_probability_of_high_LoF_intolerance_by_predLoF-CpG"
write_csv(supp_table_all_predictions, "pred_LoF_CpG_predictions.csv")


#####table with promoter annotations before and after

load(file = "proms_list_with_canonical_with_mouse_info_no_criteria.rda")
proms_list_with_canonical <- endoapply(proms_list_with_canonical, function(xx) {
  xx[which(xx$tx_biotype == "protein_coding")]
})

proms_unlisted <- unlist(proms_list_with_canonical)
proms_only_canonical <- proms_unlisted[which(proms_unlisted$canonical == TRUE)]

proms_good_power_ids_with_non_canonical_prom <- proms_good_power[which(proms_good_power$canonical == FALSE)]$gene_id
proms_low_power_ids_with_non_canonical_prom <- proms_low_power[which(proms_low_power$canonical == FALSE)]$gene_id

proms_with_non_canonical_prom <- c(proms_good_power_ids_with_non_canonical_prom, proms_low_power_ids_with_non_canonical_prom)

proms_canonical_and_alt <- data.frame(gene_id = proms_with_non_canonical_prom, 
      gene_name = as.character(sapply(proms_with_non_canonical_prom, 
                     function(xx) proms_to_use$gene_name[which(proms_to_use$gene_id == xx)])),                               
      tx_id_canonical = as.character(sapply(proms_with_non_canonical_prom, 
                               function(xx) proms_only_canonical$tx_id[which(proms_only_canonical$gene_id == xx)])), 
      chr = as.character(sapply(proms_with_non_canonical_prom, 
                   function(xx) as.character(seqnames(proms_only_canonical[which(proms_only_canonical$gene_id == xx)])))),
      promoter_start_coord_canonical = as.numeric(sapply(proms_with_non_canonical_prom, 
                                     function(xx) start(proms_only_canonical[which(proms_only_canonical$gene_id == xx)]))), 
      promoter_end_coord_canonical = as.numeric(sapply(proms_with_non_canonical_prom, 
                                     function(xx) end(proms_only_canonical[which(proms_only_canonical$gene_id == xx)]))),
      tx_id_used = as.character(sapply(proms_with_non_canonical_prom, 
                               function(xx) proms_to_use$tx_id[which(proms_to_use$gene_id == xx)])), 
      promoter_start_coord_used = as.numeric(sapply(proms_with_non_canonical_prom, 
                                     function(xx) start(proms_to_use[which(proms_to_use$gene_id == xx)]))), 
      promoter_end_coord_used = as.numeric(sapply(proms_with_non_canonical_prom, 
                                   function(xx) end(proms_to_use[which(proms_to_use$gene_id == xx)]))))

write_csv(proms_canonical_and_alt, "promoter_coordinates_canonical_and_used.csv")





