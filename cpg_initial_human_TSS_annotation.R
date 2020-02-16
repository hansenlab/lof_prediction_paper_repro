#########################
#####the purpose of this set of scripts is to systematically look at the TSS coordinates provided by ensembl, 
#####and classify a subset of them as "high-confidence" TSSs. 

##Dependencies
library(BSgenome.Hsapiens.UCSC.hg19)
library(EnsDb.Hsapiens.v75)
library(rtracklayer)
library(AnnotationHub)
library(readr)
##


##first load the gnomAD data frame 
constraint <- read_tsv('constraint.txt')
constraint <- as(constraint, "DataFrame")

##get all human promoters in Ensembl
edb <- EnsDb.Hsapiens.v75
proms_all <- promoters(edb, upstream = 2000, downstream = 2000, 
                       columns = c("gene_name", "tx_id", "tx_cds_seq_start", "tx_cds_seq_end",
                                   "tx_biotype", "gene_id"))
genome(seqinfo(proms_all)) <- "hg19"
seqlevelsStyle(proms_all) <- "ucsc"

chrs <- names(Hsapiens)[1:22]
proms_all <- proms_all[which(seqnames(proms_all) %in% chrs[1:22])] #exclude promoters on the sex chromosomes because it is hard to talk about LoF-intolerance for those
seqlevels(proms_all) <- seqlevels(proms_all)[1:22]
##

##now for each of these promoters, add a new column corresponding to the number of ENCODE ChIP-seq experiments
##where that promoter has a POLR2A peak and another column corresponding to the number of 
##ENCODE experiments where that promoter has a DHS
session <- browserSession("UCSC")
genome(session) <- "hg19"

#POLR2A
query <- ucscTableQuery(session, "Txn Factor ChIP")
tableName(query) <- "wgEncodeRegTfbsClusteredV3"
tf_table <- getTable(query)
polr2a_table <- tf_table[which(tf_table$name == "POLR2A"), ]
colnames(polr2a_table)[c(2,3,4)] <- c("chr", "start", "end")
polr2a_table <- polr2a_table[, c("chr", "start", "end", "name", "expCount")]
polr2a <- makeGRangesFromDataFrame(polr2a_table, keep.extra.columns = TRUE)
genome(polr2a) <- "hg19"

overlaps <- findOverlaps(proms_all, polr2a)
overlaps_counts <- polr2a$expCount[subjectHits(overlaps)]
peak_counts <- sapply(split(overlaps_counts, queryHits(overlaps)), max)
proms_all$polr2a_peak_count <- 0
proms_all$polr2a_peak_count[as.integer(names(peak_counts))] <- peak_counts

#DHS
query <- ucscTableQuery(session, "DNase Clusters")
tableName(query) <- "wgEncodeRegDnaseClusteredV3"
dhs_table <- getTable(query)
colnames(dhs_table)[c(2,3,4)] <- c("chr", "start", "end")
dhs_table <- dhs_table[, c("chr", "start", "end", "name")]
dhs <- makeGRangesFromDataFrame(dhs_table, keep.extra.columns = TRUE)
genome(dhs) <- "hg19"

overlaps <- findOverlaps(proms_all, dhs)
overlaps_counts <- dhs$name[subjectHits(overlaps)]
peak_counts <- sapply(split(overlaps_counts, queryHits(overlaps)), max)
proms_all$dhs_peak_count <- 0
proms_all$dhs_peak_count[as.integer(names(peak_counts))] <- peak_counts

#now add the constraint info as metadata to the proms_all granges object
anyDuplicated(constraint$transcript)
anyDuplicated(proms_all$tx_id)

rownames(constraint) <- constraint$transcript
names(proms_all) <- proms_all$tx_id

merged_df  <- merge(as.data.frame(values(proms_all)),
                    as.data.frame(constraint[, c("canonical", "oe_lof", "oe_lof_upper", "exp_lof")]),
                    by = 0, sort = FALSE, all.x = TRUE)
rownames(merged_df) <- merged_df$Row.names

merged_df$Row.names <- NULL
values(proms_all) <- merged_df[names(proms_all), ]
proms_all$oe_cpg_ratio <- getOECpGRatio(proms_all, Hsapiens) #function defined in the "cpg_functions" script

#also add metadata column for CpG island overlap
hub <- AnnotationHub()
query(hub, c("cpg","hg19"))
cpg_islands <- hub[["AH5086"]]
genome(seqinfo(cpg_islands)) <- "hg19"
seqlevelsStyle(cpg_islands) <- "ucsc"
cpg_islands <- cpg_islands[which(seqnames(cpg_islands) %in% unique(seqnames(proms_all)))]
seqlevels(cpg_islands) <- seqlevels(cpg_islands)[1:22]

proms_all$overlaps_cpg_island <- FALSE
proms_all$overlaps_cpg_island[unique(queryHits(findOverlaps(proms_all, cpg_islands)))] <- TRUE
#

#add some extra metadata columns that will be important downstream
proms_all$TSS <- start(proms_all) + 1999
proms_all$use <- FALSE #this will be updated when the criteria for high-confidence annotation are applied
proms_all$is_intolerant <- FALSE
proms_all$is_intolerant[which(proms_all$oe_lof_upper < 0.35)] <- TRUE
proms_all$distance_from_closest_mouse_TSS <- NA #this will be updated when the mouse info is added on top
#
save(proms_all, file = "proms_all.rda")

##the following loop excludes promoters that are found in subtelomeric regions 
##(that is, regions near the two ends of chromosomes - here defined as 2 Mb from the start/end of a chromosome). 
##That's because the CpG islands in subtelomeric regions are different than the CpG islands in the rest of the genome (Tanay's paper states that as well), 
##and therefore their presence will not be used as a marker of promoters 
##(and we also don't expect it to be related to LoF-intolerance) 
proms <- proms_all[-which(proms_all$gene_id %in% proms_all$gene_id[which(start(proms_all) < 2000000)])]
for (i in 1:22){
  indices <- which(end(proms) > (length(Hsapiens[[i]])-2000000) & seqnames(proms) == paste0("chr", i))
  if (length(indices) > 0){
    proms <- proms[-which(proms$gene_id %in% proms$gene_id[indices])]
  }
  proms
}

proms_list <- split(proms, proms$gene_id)

has_canonical <- sapply(proms_list, function(xx) {
  TRUE %in% xx$canonical
})
proms_list_with_canonical <- proms_list[which(has_canonical == TRUE)]

proms_list_with_canonical <- endoapply(proms_list_with_canonical, function(xx) {
  xx$is_max_polr2a_signal <- FALSE
  xx$is_max_polr2a_signal[which(xx$polr2a_peak_count == max(xx$polr2a_peak_count))] <- TRUE
  
  xx$is_max_oe_cpg <- FALSE
  xx$is_max_oe_cpg[which(xx$oe_cpg_ratio == max(xx$oe_cpg_ratio))] <- TRUE
  
  xx$is_max_dhs_signal <- FALSE
  xx$is_max_dhs_signal[which(xx$dhs_peak_count == max(xx$dhs_peak_count))] <- TRUE
  
  xx$diff_from_max_polr2a_signal <- 0
  xx$diff_from_max_polr2a_signal[which(xx$is_max_polr2a_signal == FALSE)] <- xx$polr2a_peak_count[
    which.max(xx$polr2a_peak_count)] - xx$polr2a_peak_count[which(xx$is_max_polr2a_signal == FALSE)]
  
  xx$diff_from_max_dhs_signal <- 0
  xx$diff_from_max_dhs_signal[which(xx$is_max_dhs_signal == FALSE)] <- xx$dhs_peak_count[
    which.max(xx$dhs_peak_count)] - xx$dhs_peak_count[which(xx$is_max_dhs_signal == FALSE)]
  xx
})


save(proms_list_with_canonical, file = "proms_list_with_canonical_no_mouse_info_no_criteria.rda")