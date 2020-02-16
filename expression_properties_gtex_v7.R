library(CePa)
library(matrixStats)
#expr_mat <- read.gct('GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct') #version 8 tx-level
#expr_mat <- read_tsv('GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt.gz') #version 7 tx-level
#expr_mat <- read.gct('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct') #version 6 tx-level

rownames(expr_mat) <- gsub("[.].*", "", rownames(expr_mat))
#save(expr_mat, file = "tx_level_expr_mat.rda")
save(expr_mat, file = "gene_level_expr_mat.rda")

#gtex_samples <- read.delim("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", stringsAsFactors = FALSE) v8
gtex_samples <- read.delim("GTEx_v7_Annotations_SampleAttributesDS.txt", stringsAsFactors = FALSE)
gtex_samples$SAMPID <- gsub("-", ".", gtex_samples$SAMPID)
#for some reason, the annotation file contains sample IDs that are not in the gtex matrix column names. remove those sample ids
notinmat <- which(!(gtex_samples$SAMPID %in% colnames(expr_mat)))
sampids <- gtex_samples$SAMPID[-notinmat]             
tissueids <- gtex_samples$SMTSD[-notinmat]
tissue_name <- sapply(colnames(expr_mat), function(xx) tissueids[which(sampids == xx)])

gtex_tissues <- unique(tissue_name)
sapply(gtex_tissues, function(xx) length(which(tissue_name == xx))) #check how many samples in each tissue


###get a list of expression matrices for each tissue
expr_mat_list <- lapply(gtex_tissues, function(xx) expr_mat[, which(tissue_name == xx)])
names(expr_mat_list) <- gtex_tissues

#get matrices for the median expression, and the variance in expression for every gene within each tissue
expr_median <- lapply(expr_mat_list, function(xx) rowMedians(xx))
expr_median_matrix <- do.call(cbind, expr_median)
rownames(expr_median_matrix) <- rownames(expr_mat)
expr_median_df <- as.data.frame(expr_median_matrix)
expr_median_df$`Brain - Cerebellum median` <- rowMedians(expr_median_matrix[, 
                                                                        grep("Brain - Cerebell", colnames(expr_median_matrix))]) #take the median of the two cerebellum measurements
expr_median_df$`Brain - Median expression` <- rowMedians(expr_median_matrix[, grep("Brain", colnames(expr_median_matrix))[-c(4:5)]]) #I exclude cerebellum because usually it's different
expr_median_matrix <- as.matrix(expr_median_df)

expr_variance <- lapply(expr_mat_list, function(xx) {
  log_mat <- log2(xx+1)
  rowVars(log_mat)
})
expr_variance_matrix <- do.call(cbind, expr_variance)
rownames(expr_variance_matrix) <- rownames(expr_mat)
expr_variance_df <- as.data.frame(expr_variance_matrix)
expr_variance_df$`Brain - Median variance` <- rowMedians(expr_variance_matrix[, grep("Brain", colnames(expr_variance_matrix))[-c(4:5)]])
expr_variance_matrix <- as.matrix(expr_variance_df)


#get tissue specificity tau
getTau <- function(medians_matrix){ #tau is thought to be the most robust measure of tissue specificity
  x_i_hat_vector_list <- lapply(1:(dim(medians_matrix)[2]), function(xx) {
    tissue_vector <- medians_matrix[, xx]
    x_i_hat_vector <- tissue_vector/rowMaxs(medians_matrix)
    x_i_hat_vector #this is a vector containing the x_i_hat values for all genes
  })
  x_i_hat_matrix <- do.call(cbind, x_i_hat_vector_list)
  tau_numerators <- rowSums(1 - x_i_hat_matrix)
  tau_vector <- tau_numerators/(dim(medians_matrix)[2] - 1) #this gives a vector of taus, for each gene
  tau_vector
}

taus <- getTau(log2(expr_median_matrix[, -head(grep("Brain", colnames(expr_median_matrix)), -2)]+1)) 
#calling the head function in the above with -2 as the offset means that I compute the specificity using the median expression across
#all brain tissues (minus the two cerebellum tissues) and the median across the two cerebellum tissues, 
#and not the expression in each of these tissues separately

tissue_with_max_expression <- apply(expr_median_matrix, 1, function(xx) colnames(expr_median_matrix)[which.max(xx)])
tissue_with_second_max_expression <- apply(expr_median_matrix, 1, 
                                           function(xx) colnames(expr_median_matrix)[order(xx, decreasing = TRUE)[2]])

median_expr_across_all_tissues <- rowMedians(expr_median_matrix[, 1:53])
median_expr_in_maximum_expression_tissue <- sapply(1:dim(expr_median_matrix)[1], function(xx) {
  tissue_index <- which(colnames(expr_median_matrix) == tissue_with_max_expression[xx])
  expr_median_matrix[xx, tissue_index]                      
})
median_expr_in_brain <- expr_median_matrix[, 54]

median_var_across_all_tissues <- rowMedians(expr_variance_matrix[, 1:53])
var_in_maximum_expression_tissue <- sapply(1:dim(expr_median_matrix)[1], function(xx) {
  tissue_index <- which(colnames(expr_variance_matrix) == tissue_with_max_expression[xx])
  expr_variance_matrix[xx, tissue_index]                      
})
var_in_brain <- expr_variance_matrix[, 54]

n_tissues_with_detectable_expression <- apply(expr_median_matrix[, 1:53], 1, function(xx) length(which(xx > 0.3))) #the 0.3 cutoff is what's used in the gnomAD paper

expr_properties_df <- data.frame(median_expr_across_all_tissues = log2(median_expr_across_all_tissues+1), 
                                 median_expr_in_maximum_expression_tissue = log2(median_expr_in_maximum_expression_tissue+1), 
                                 median_expr_in_brain = log2(median_expr_in_brain+1),
                                 median_var_across_all_tissues = median_var_across_all_tissues, 
                                 var_in_maximum_expression_tissue = var_in_maximum_expression_tissue, 
                                 var_in_brain = var_in_brain, 
                                 tissue_specificity = taus, 
                                 n_tissues_with_detectable_expression = n_tissues_with_detectable_expression,
                                 tissue_with_max_expression = tissue_with_max_expression, 
                                 tissue_with_second_max_expression = tissue_with_second_max_expression)
rownames(expr_properties_df) <- rownames(expr_mat)
save(expr_properties_df, file = "expression_properties_df_gene_level.rda")
save(expr_median_matrix, expr_variance_matrix, file = "median_expression_and_expression_variance_matrices_gene_level.rda")








