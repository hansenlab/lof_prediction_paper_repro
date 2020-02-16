####
load(file = "proms_list_with_canonical_with_mouse_info_no_criteria.rda")

#also load the data frama containing the expression properties for each ENSEMBL gene id
#we will use tissue specificity for the annotation that will follow
load(file = "expression_properties_df_gene_level.rda")
expr_properties_df_with_canonical <- expr_properties_df[which(rownames(expr_properties_df) %in% 
                                                                names(proms_list_with_canonical)), ]
proms_list_with_canonical <- proms_list_with_canonical[which(names(proms_list_with_canonical) %in% 
                                                               rownames(expr_properties_df_with_canonical))]
#the following is necessary, and ensures the order of gene ids in the expression df is the same as the order in the list of promoters
expr_properties_df_with_canonical <- expr_properties_df_with_canonical[names(proms_list_with_canonical), ] 


#define some agreement criteria that will be used
pol_canonical_agree <- sapply(proms_list_with_canonical, function(xx) {
  xx$is_max_polr2a_signal[which(xx$canonical == TRUE)]
})
cpg_island_canonical_agree <- sapply(proms_list_with_canonical, function(xx) {
  xx$overlaps_cpg_island[which(xx$canonical == TRUE)]
})
canonical_polr2a <- sapply(proms_list_with_canonical, function(xx) {
  xx$polr2a_peak_count[which(xx$canonical == TRUE)]
})
max_polr2a <- sapply(proms_list_with_canonical, function(xx) {
  xx$polr2a_peak_count[which(xx$is_max_polr2a_signal == TRUE)[1]]
})
has_cpg_island <- sapply(proms_list_with_canonical, function(xx) {
  TRUE %in% xx$overlaps_cpg_island
})


####now annotate
#case 1
case_1_indices <- which(pol_canonical_agree == TRUE & 
                   expr_properties_df_with_canonical$tissue_specificity < 0.6 & 
                   canonical_polr2a > 35)
proms_list_with_canonical[case_1_indices] <- endoapply(proms_list_with_canonical[case_1_indices], function(xx) {
       xx$use[which(xx$canonical == TRUE)] <- TRUE
       xx
})


#case 2
case_2_indices <- which(pol_canonical_agree == FALSE & 
                   expr_properties_df_with_canonical$tissue_specificity < 0.6 & 
                   canonical_polr2a < 10 & max_polr2a > 35)
proms_list_with_canonical[case_2_indices] <- endoapply(proms_list_with_canonical[case_2_indices], function(xx) {
       xx$use[which(xx$is_max_polr2a_signal == TRUE)] <- TRUE
       xx
})


#case 3
case_3_indices <- which(expr_properties_df_with_canonical$tissue_specificity > 0.6 & 
                   cpg_island_canonical_agree == TRUE & canonical_polr2a < 10)
proms_list_with_canonical[case_3_indices] <- endoapply(proms_list_with_canonical[case_3_indices], function(xx) {
       xx$use[which(xx$canonical == TRUE)] <- TRUE
       xx                                    
})


#case 4
case_4_indices <- which(expr_properties_df_with_canonical$tissue_specificity > 0.6 & 
                          has_cpg_island == FALSE & canonical_polr2a < 10)
proms_list_with_canonical[case_4_indices] <- endoapply(proms_list_with_canonical[case_4_indices], function(xx) {
       xx$use[which(xx$canonical == TRUE & xx$distance_from_closest_mouse_TSS <= 500)] <- TRUE
       xx                                    
})


#case 5
unique_prom_ranges_count <- sapply(proms_list_with_canonical, function(xx) length(unique(xx)))
case_5_indices <- which(expr_properties_df_with_canonical$tissue_specificity > 0.6 & 
                          has_cpg_island == FALSE & canonical_polr2a < 10 & 
                          unique_prom_ranges_count == 1)
proms_list_with_canonical[case_5_indices] <- endoapply(proms_list_with_canonical[case_5_indices], function(xx) {
       if (is.na(xx$distance_from_closest_mouse_TSS)) {
        xx$use <- TRUE
       }                                          
       xx                                    
})


proms_with_canonical_to_use <- unlist(proms_list_with_canonical)
proms_with_canonical_to_use <- proms_with_canonical_to_use[which(proms_with_canonical_to_use$use == TRUE)]
proms_with_canonical_to_use_list <- split(proms_with_canonical_to_use, proms_with_canonical_to_use$gene_id)

#now in some instances corresponding to case 2 above, there are multiple alternative promoters
#satisfying the criterion, so we pick the one with the trascript that has the highest exp_lof in gnomAD
proms_with_canonical_to_use_list[which(lengths(proms_with_canonical_to_use_list) > 1)] <- endoapply(
  proms_with_canonical_to_use_list[which(lengths(proms_with_canonical_to_use_list) > 1)], function(xx) {
    if (TRUE %in% xx$overlaps_cpg_island){
      xx[which(xx$overlaps_cpg_island == TRUE)][which.max(xx[which(xx$overlaps_cpg_island == TRUE)]$exp_lof)]
    } else {
      xx[which.max(xx$exp_lof)] 
    }
})

save(proms_with_canonical_to_use_list, file = "proms_with_canonical_to_use_list.rda")




