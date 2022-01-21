#' @include spatial_cor.R
#' @include iSpatial.R
NULL

#' 10 fold cross validation
#'
#'
#' @param spRNA seurat object of spatial transcriptome
#' @param scRNA seurat object of scRNA-seq
#' @param fold number of X fold crossing validation
#' @param correct.spRNA Whether to stabilize expression in spRNA.
#' @param correct.scRNA Whether to stabilize expression in scRNA.
#' 
#' 
#' @return correlation between inferred spatial expression and ground true.
#' 
#' @export

cross_validation = function(spRNA, scRNA, fold = 10, correct.spRNA = TRUE, correct.scRNA = FALSE){
  target_genes = intersect(rownames(spRNA), rownames(scRNA))
  gene_number = length(target_genes)
  
  set.seed(3)
  gene_group = suppressWarnings(base::split(sample(target_genes), 1:as.numeric(fold))) # group genes into 10 
  
  spatial_cor_10x_cross = lapply(gene_group, function(x){
    validation_gene = x
    train_gene = setdiff(target_genes, validation_gene)
    
    # train
    spRNA_group = iSpatial(spRNA[train_gene, ], scRNA, correct.spRNA = correct.spRNA, correct.scRNA = correct.scRNA)
    
    # validation
    suppressWarnings(spatial_cor(spRNA_group, "enhanced", spRNA, "RNA", feature = validation_gene))
  })
  
  spatial_cor_10x_cross
}
