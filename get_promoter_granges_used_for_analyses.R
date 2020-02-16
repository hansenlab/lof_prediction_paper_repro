###
library(readr)
library(rtracklayer)
library(EnsDb.Hsapiens.v75)
library(phastCons100way.UCSC.hg19)
###
####load the promoters with with high-confidence annotation
load(file = "proms_with_canonical_to_use_list.rda")

proms_to_use <- unlist(proms_with_canonical_to_use_list)
names(proms_to_use) <- proms_to_use$tx_id
proms_to_use$high_confidence_tolerant <- FALSE
proms_to_use$high_confidence_tolerant[which(proms_to_use$oe_lof > 0.8 & proms_to_use$exp_lof >= 20)] <- TRUE
proms_to_use$low_power_likely_intolerant <- FALSE
proms_to_use$low_power_likely_intolerant[which(proms_to_use$oe_lof < 0.5 & proms_to_use$exp_lof < 20)] <- TRUE

gsco <- phastCons100way.UCSC.hg19

##add promoter phast cons info
proms_to_use$score <- gscores(gsco, proms_to_use)$default

##add exonic phast cons info
edb <- EnsDb.Hsapiens.v75
all_coding_exons <- cdsBy(edb, filter = SeqNameFilter(c(1:22)), columns = "tx_id")
all_coding_exons <- unlist(all_coding_exons)
exons_to_use <- all_coding_exons[which(all_coding_exons$tx_id %in% proms_to_use$tx_id)]

exons_to_use$score <- gscores(gsco, exons_to_use)$default
exons_by_tx <- split(exons_to_use, exons_to_use$tx_id)
exons_by_tx <- exons_by_tx[which(names(exons_by_tx) %in% proms_to_use$tx_id)]
exons_by_tx <- exons_by_tx[proms_to_use$tx_id]

proms_to_use$exon_score <- as.numeric(sapply(exons_by_tx, function(xx) {
  gr <- xx
  widths <- width(gr)
  weighted_average_numerator <- sum(gr$score*widths)
  weighted_average_denominator <- sum(widths)
  weighted_average_numerator/weighted_average_denominator
}))

###add other haploinsufficiency scores. this will be used later for the comparison
hapl_predictions <- import.bed("HI_Predictions_Version3.bed") #haploinsufficiency predictions from Huang et al. 2010
hapl_predictions$gene_name <- gsub("[|].*", "", hapl_predictions$name)
#hapl_predictions <- hapl_predictions[which(hapl_predictions$gene_name %in% proms_to_use$gene_name)]
hapl_predictions <- hapl_predictions[-which(hapl_predictions$gene_name == "CKS1B")] #that's a name that appears twice. Since the LOEUF assignment will be done using the gene names, we cannot trust this

proms_to_use$haplo_score <- NA
proms_to_use$haplo_score[which(proms_to_use$gene_name %in% 
                                 hapl_predictions$gene_name)] <- unlist(sapply(proms_to_use$gene_name[
                                   which(proms_to_use$gene_name %in% hapl_predictions$gene_name)], 
                                   function(xx) hapl_predictions$score[which(hapl_predictions$gene_name == xx)]))

###########
hapl_predictions <- read_csv('haplo_prediction_nar.csv') #haploinsufficiency predictions from Steinberg et al., 2015
colnames(hapl_predictions)[1] <- "gene_id"

proms_to_use$haplo_score_2 <- NA
proms_to_use$haplo_score_2[which(proms_to_use$gene_id %in% 
                                   hapl_predictions$gene_id)] <- unlist(sapply(proms_to_use$gene_id[
                                     which(proms_to_use$gene_id %in% hapl_predictions$gene_id)], 
                                     function(xx) hapl_predictions$GHIS[which(hapl_predictions$gene_id == xx)]))


###########
hapl_predictions <- read_csv('haplo_prediction_natcom.csv') #haploinsufficiency predictions from Han et al., 2018
colnames(hapl_predictions)[1] <- "gene_id" #that's a name that appears twice. Since the LOEUF assignment will be done using the gene names, we cannot trust this

proms_to_use$haplo_score_3 <- NA
proms_to_use$haplo_score_3[which(proms_to_use$gene_id %in% 
                                   hapl_predictions$gene_id)] <- unlist(sapply(proms_to_use$gene_id[
                                     which(proms_to_use$gene_id %in% hapl_predictions$gene_id)], 
                                     function(xx) hapl_predictions$Episcore[which(hapl_predictions$gene_id == xx)]))


###add metadata columns about expression properties
load(file = "expression_properties_df_gene_level.rda")
expr_properties <- expr_properties_df[which(rownames(expr_properties_df) %in% proms_to_use$gene_id), ]

proms_to_use$tau <- sapply(proms_to_use$gene_id, function(xx) 
  expr_properties$tissue_specificity[which(rownames(expr_properties) == xx)])

proms_to_use$median_expr_in_maximum_expression_tissue <- sapply(proms_to_use$gene_id, function(xx) 
  expr_properties$median_expr_in_maximum_expression_tissue[which(rownames(expr_properties) == xx)])

proms_to_use$median_expr_across_all_tissues <- sapply(proms_to_use$gene_id, function(xx) 
  expr_properties$median_expr_across_all_tissues[which(rownames(expr_properties) == xx)])

proms_to_use$n_tissues_with_detectable_expression <- sapply(proms_to_use$gene_id, function(xx) 
  expr_properties$n_tissues_with_detectable_expression[which(rownames(expr_properties) == xx)])

proms_to_use$n_tissues_with_detectable_expression <- sapply(proms_to_use$gene_id, function(xx) 
  expr_properties$n_tissues_with_detectable_expression[which(rownames(expr_properties) == xx)])


###add info on EZH2 binding
session <- browserSession("UCSC")
genome(session) <- "hg19"


query <- ucscTableQuery(session, "Txn Factor ChIP")
tableName(query) <- "wgEncodeRegTfbsClusteredV3"
tf_table <- getTable(query)
ezh2_table <- tf_table[which(tf_table$name == "EZH2"), ]

colnames(ezh2_table)[c(2,3,4)] <- c("chr", "start", "end")
ezh2_table <- ezh2_table[, c("chr", "start", "end", "name", "expCount")]
ezh2 <- makeGRangesFromDataFrame(ezh2_table, keep.extra.columns = TRUE)
genome(ezh2) <- "hg19"

overlaps <- findOverlaps(proms_to_use, ezh2)
overlaps_counts <- ezh2$expCount[subjectHits(overlaps)]
peak_counts <- sapply(split(overlaps_counts, queryHits(overlaps)), max)
proms_to_use$ezh2_peak_count <- 0
proms_to_use$ezh2_peak_count[as.integer(names(peak_counts))] <- peak_counts
###


###get the granges with the transcripts with sufficient power to estimate LOEUF
proms_good_power <- proms_to_use[which(proms_to_use$exp_lof >= 20)]
proms_good_power <- keepNonOverlappingRanges(proms_good_power)

overlapping_low_power_intolerant <- proms_good_power$tx_id[
  unique(queryHits(findOverlaps(proms_good_power, proms_to_use[which(proms_to_use$low_power_likely_intolerant == TRUE)])))]

proms_good_power <- proms_good_power[-which(proms_good_power$oe_lof_upper > 0.5 & 
                                              proms_good_power$tx_id %in% overlapping_low_power_intolerant)]

###get the granges with the transcripts with low power to estimate LOEUF
proms_low_power <- proms_to_use[which(proms_to_use$exp_lof <= 10)]
proms_low_power <- proms_low_power[-unique(queryHits(findOverlaps(proms_low_power, 
                                                                  proms_to_use[which(proms_to_use$exp_lof > 10 & proms_to_use$high_confidence_tolerant == FALSE)], 
                                                                  ignore.strand = TRUE)))]
#now exclude the following genes because they are likely to represent false annotation of tx end (because there is a different isoform with more exp_lof, therefore bigger - see end of this script)
max_exp_lof_for_each_unascertained_gene <- sapply(proms_low_power$gene_name, 
                                                  function(xx) max(constraint$exp_lof[which(constraint$gene == xx)])) #can do this based on gene_name because I have checked it is unique
diff_from_exp_lof_used <- max_exp_lof_for_each_unascertained_gene - proms_low_power$exp_lof
which(diff_from_exp_lof_used > 20)
proms_low_power <- proms_low_power[-which(proms_low_power$gene_name %in% names(which(diff_from_exp_lof_used > 20)))]

proms_low_power <- keepNonOverlappingRanges2(proms_low_power)




