#####the goal is to get the TSS coordinates of all mouse refesq transcripts corresponding to mouse genes 
#####with human homologs that I use in the list of promoters.
#####I need to do sequential biomaRt queries because 
#####a) I first need to get mouse homologs of human transcripts, and the use the mouse database to get refseq tx ids for these
#####b) biomaRt doesn't allow simultaneous queries for 
#####both transcript-level features (refseq tx ids) and gene-level features (human homologs)

#load dependencies
library(rtracklayer)
library(biomaRt)
library(dplyr)
library(purrr)

#
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

load(file = "proms_list_with_canonical_no_mouse_info_no_criteria.rda")
human_mouse_homologs <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                             "mmusculus_homolog_ensembl_gene", 
                                             "mmusculus_homolog_associated_gene_name", 
                                             "mmusculus_homolog_orthology_confidence", 
                                             "mmusculus_homolog_perc_id_r1"),
                              filters = "ensembl_gene_id",
                              values = names(proms_list_with_canonical),
                              mart = human)


human_mouse_homologs <- human_mouse_homologs[
  which(human_mouse_homologs$mmusculus_homolog_orthology_confidence == 1), ]

mouse_refseq_df <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                        "refseq_mrna"),
                         filters = "ensembl_gene_id",
                         values = unique(human_mouse_homologs$mmusculus_homolog_ensembl_gene),
                         mart = mouse)
mouse_refseq_df <- mouse_refseq_df[-which(mouse_refseq_df$refseq_mrna == ""), ]
mouse_refseq_list <- split(mouse_refseq_df, mouse_refseq_df$ensembl_gene_id)

mouse_homology_df <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                          "hsapiens_homolog_ensembl_gene", 
                                          "hsapiens_homolog_orthology_confidence", 
                                          "hsapiens_homolog_perc_id_r1"),
                           filters = "ensembl_gene_id",
                           values = unique(mouse_refseq_df$ensembl_gene_id),
                           mart = mouse)


mouse_homology_df <- mouse_homology_df[which(mouse_homology_df$hsapiens_homolog_orthology_confidence == 1), ]
mouse_homology_df <- mouse_homology_df[which(mouse_homology_df$hsapiens_homolog_ensembl_gene %in% names(proms_list_with_canonical)), ]
mouse_homology_list <- split(mouse_homology_df, mouse_homology_df$ensembl_gene_id)
mouse_homology_list <- mouse_homology_list[names(mouse_refseq_list)]

#now merge the two lists
human_mouse_homology_df <- map2_df(mouse_refseq_list, mouse_homology_list, inner_join, by = "ensembl_gene_id")
human_mouse_homology_df <- human_mouse_homology_df[-which(human_mouse_homology_df$refseq_mrna %in% 
                                      human_mouse_homology_df$refseq_mrna[duplicated(human_mouse_homology_df$refseq_mrna)]), ]
colnames(human_mouse_homology_df)[which(colnames(human_mouse_homology_df) == "refseq_mrna")] <- "refseq_id"

#now get Tx coords for these refseq ids
session <- browserSession()
genome(session) <- "hg19"

query <- ucscTableQuery(session, "Other RefSeq")
tableName(query) <- "xenoRefGene"
coords_table <- getTable(query)


coords_table_mouse <- coords_table[, c("name", "name2", "chrom", "strand", "txStart", "txEnd")]
colnames(coords_table_mouse) <- c("refseq_id", "gene_name", "chr", "strand", "start", "end")
coords_table_mouse$refseq_id <- as.character(coords_table_mouse$refseq_id)
coords_table_mouse <- coords_table_mouse[which(coords_table_mouse$refseq_id %in% human_mouse_homology_df$refseq_id), ]
coords_table_mouse <- coords_table_mouse[-which(coords_table_mouse$refseq_id %in% 
                         coords_table_mouse$refseq_id[duplicated(coords_table_mouse$refseq_id)]), ] #remove refseq ids which don't have unique positions

save(human_mouse_homology_df, coords_table_mouse, file = "mouse_refseq_hg19_coords_objects.rda") #save so we don't have to redo these biomart/ucsc queries

mouse_refseq_coords_list <- split(coords_table_mouse, coords_table_mouse$refseq_id)
human_mouse_homology_df_list <- split(human_mouse_homology_df, human_mouse_homology_df$refseq_id)
human_mouse_homology_df_list <- human_mouse_homology_df_list[which(names(human_mouse_homology_df_list) %in% 
                                                                     names(mouse_refseq_coords_list))]
mouse_refseq_coords_list <- mouse_refseq_coords_list[names(human_mouse_homology_df_list)]


mouse_refseq_tx_human_coords <- map2_df(human_mouse_homology_df_list, mouse_refseq_coords_list, inner_join, by = "refseq_id")
mouse_refseq_tx_human_coords$TSS <- NA
mouse_refseq_tx_human_coords$TSS[which(mouse_refseq_tx_human_coords$strand == "+")] <- mouse_refseq_tx_human_coords$start[
  which(mouse_refseq_tx_human_coords$strand == "+")]
mouse_refseq_tx_human_coords$TSS[which(mouse_refseq_tx_human_coords$strand == "-")] <- mouse_refseq_tx_human_coords$end[
  which(mouse_refseq_tx_human_coords$strand == "-")]


mouse_refseq_tx_human_coords <- makeGRangesFromDataFrame(mouse_refseq_tx_human_coords, keep.extra.columns = TRUE, 
                                                         starts.in.df.are.0based = TRUE)

human_genes_mouse_refseq_TSS <- split(mouse_refseq_tx_human_coords, 
                                      mouse_refseq_tx_human_coords$hsapiens_homolog_ensembl_gene)

proms_list_with_canonical[which(names(proms_list_with_canonical) 
                                %in% names(human_genes_mouse_refseq_TSS))] <- endoapply(proms_list_with_canonical[
                                  which(names(proms_list_with_canonical)  %in% names(human_genes_mouse_refseq_TSS))], 
                         function(xx) {
        xx$distance_from_closest_mouse_TSS <- sapply(1:length(xx), 
                  function(x) min(abs(xx$TSS[x] - human_genes_mouse_refseq_TSS[which(names(human_genes_mouse_refseq_TSS) == xx$gene_id[1])][[1]]$TSS)))
        xx
})

save(proms_list_with_canonical, file = "proms_list_with_canonical_with_mouse_info_no_criteria.rda")

