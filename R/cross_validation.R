#' @include spatial_cor.R
#' @include iSpatial.R
NULL

#' 10 fold cross validation
#'
#'
#' @param spRNA seurat object of spatial transcriptome
#' @param scRNA seurat object of scRNA-seq
#' 
#' @return correlation between inferred spatial expression and ground true.
#' 
#' @export

cross_validation = function(spRNA, scRNA){
  
  target_genes = rownames(spRNA)
  gene_number = length(target_genes)
  
  set.seed(3)
  gene_group = base::split(sample(target_genes), 1:10) # group genes into 10 
  
  spatial_cor_10x_cross = lapply(gene_group, function(x){
    validation_gene = x
    train_gene = setdiff(target_genes, validation_gene)
    
    # train
    spRNA_group = iSpatial(spRNA[train_gene, ], scRNA)
    
    # validation
    spatial_cor(spRNA_group, "enhanced", spRNA, "RNA", feature = validation_gene)
  })
  
  spatial_cor_10x_cross
}