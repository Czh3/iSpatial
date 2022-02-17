#' @import Seurat
#' @import SeuratObject
#' @import ggplot2
#' @import harmony
#' @import Matrix
#'
NULL

############# infer whole genome spatial transcriptome from single cell RNAseq

#' run harmony integration
#' 
#' This function calls harmony package(https://github.com/immunogenomics/harmony) 
#' to integrate the scRNA-seq and spatial transcriptome data.
#' 
#' @param seurat_obj merged seurat object
#' @param genes_select genes used to integration, usually targeted genes.
#' @param npc number of PC
#' @return seurat object
#' 
run_harmony = function(seurat_obj, genes_select, npc = 20){
  SeuratObject::VariableFeatures(seurat_obj) = genes_select
  seurat_obj = Seurat::ScaleData(seurat_obj, features = genes_select, vars.to.regress = "tech", verbose = FALSE)
  seurat_obj = Seurat::RunPCA(seurat_obj, npcs = npc, features = genes_select, verbose = FALSE)
  
  # harmony
  seurat_obj <- harmony::RunHarmony(
    object = seurat_obj,
    group.by.vars = 'tech',
    plot_convergence = F,
    theta = 2,
    lambda = 0.5,
    verbose = FALSE,
    max.iter.harmony = 20
  )
  return(seurat_obj)
}

#' stabilize expression
#' 
#' This function use nearest neighbor to correct the expression. 
#' This could remove some  False Positive / False negative expressions.
#' This function corrects "RNA" assay "data" slot. 
#' 
#' @param obj seurat object
#' @param neighbor number of nearest neighbors used to correct
#' @param npcs number of principal component (PC) used to RunPCA
#' @param n.core number of CPU cores used to parallel
#' @return seurat object with corrected expression
#' 
#' @export
#' 
#' 
stabilize_expr = function(obj, neighbor = 3, npcs = 10, n.core = 10){
  
  neighbor = neighbor + 1 # include cell itself

  obj = Seurat::FindVariableFeatures(obj, selection.method = "vst", verbose = FALSE) 
  obj = Seurat::ScaleData(obj, verbose = FALSE)
  obj = Seurat::RunPCA(obj, npcs = npcs, verbose = FALSE)
  obj = suppressMessages(Seurat::FindNeighbors(obj, return.neighbor = T, k.param = neighbor, verbose = FALSE))
  
  # imputation expr value by KNN
  enhancer_expr_tmp = expm1(obj@assays$RNA@data)
  enhancer_expr = parallel::mclapply(colnames(obj), function(cell){
    cell_neighbors = Seurat::TopNeighbors(obj@neighbors$RNA.nn, cell, n=neighbor)
    0.5 * matrixStats::rowMedians(as.matrix(enhancer_expr_tmp[,cell_neighbors[-1]])) + 0.5 * enhancer_expr_tmp[,cell_neighbors[1]]
    #Matrix::rowMeans(obj@assays$RNA@data[,cell_neighbors])
  }, mc.cores = n.core)
  enhancer_expr = do.call(cbind, enhancer_expr)
  colnames(enhancer_expr) = colnames(obj)
  rownames(enhancer_expr) = rownames(obj)
  obj@assays$RNA@data = as(log1p(enhancer_expr), "dgCMatrix")
  rm(enhancer_expr, enhancer_expr_tmp)
  gc()
  return(obj)
}

#' calculate correlation
#' 
#'
#' @param x sparse matrix
#' @return correlation matrix
#' 
#' @export
#' 
sparse.cor <- function(x){
  n <- nrow(x)
  m <- ncol(x)
  ii <- unique(x@i)+1 # rows with a non-zero element
  
  Ex <- Matrix::colMeans(x)
  nozero <- as.vector(x[ii,]) - rep(Ex,each=length(ii)) 
  
  covmat <- ( Matrix::crossprod(matrix(nozero,ncol=m)) +
                Matrix::crossprod(t(Ex))*(n-length(ii))
  )/(n-1)
  sdvec <- sqrt(diag(covmat))
  covmat/Matrix::crossprod(t(sdvec))
}


#' Main function of iSpatial
#' 
#' This function integrates scRNA-seq and spatial transcriptome. With limited
#' genes in merFISH/SeqFISH, iSpatial infer/enhance all genes from scRNA-seq.
#' Finally, this function generates spatially expression of all genes.
#' 
#' @param spRNA seurat object of spatial transcriptome data
#' @param scRNA seurat object of single cell RNA-seq data
#' @param dims which dimensions to use when find the nearest neighbors
#' @param k.neighbor number of neighbors to use when infer the expression
#' @param infered.assay names of output assay in seurat object
#' @param weighted.KNN we use a dynamics weight to assign scRNA values to spRNA.
#' The dynamics weight based on the correlation between scRNA cell and spRNA 
#' cell. 
#' @param RNA.weight the weight of scRNA-seq expression when genes detected both
#' in spRNA and scRNA. RNA.weight should be 0 to 1.
#' The inferred expression = (1 - RNA.weight) * spRNA + RNA.weight * scRNA.
#' @param n.core number of CPU cores used to parallel.
#' @param correct.scRNA Whether to stabilize expression in scRNA.
#' @param correct.spRNA Whether to stabilize expression in spRNA.
#' @param correct.neighbor number of nearest neighbors used to correct expr.
#' 
#' @return returns a seurat object with inferred expression in infered.assay.
#' 
#' @export
#' 

infer_v0.1 = function(
  spRNA,
  scRNA, 
  dims = 1:30,
  k.neighbor = 30,
  infered.assay = "enhanced",
  weighted.KNN = TRUE,
  RNA.weight = 0.5,
  n.core = 10,
  correct.spRNA = TRUE,
  correct.scRNA = FALSE,
  correct.neighbor = 3
){
  spRNA$tech = "spatial"
  scRNA$tech = "scRNA"
 
  if(is.null(spRNA@assays$RNA)){
    stop(paste(spRNA, " do not have 'RNA' assay."))
  }
  
  if(is.null(scRNA@assays$RNA)){
    stop(paste(scRNA, " do not have 'RNA' assay."))
  }
  
  if(length(spRNA@assays$RNA@data) == 0){
    stop(paste(spRNA, " is not normlized. Run Seurat::NormalizeData."))
  }
  
  if(length(scRNA@assays$RNA@data) == 0){
    stop(paste(scRNA, " is not normlized. Run Seurat::NormalizeData."))
  }
  
  if(correct.spRNA){
    message("Stablize spatial transcriptome.")
    spRNA = stabilize_expr(spRNA, neighbor = correct.neighbor, n.core = n.core, npcs = length(dims))
  }
  
  if(correct.scRNA){
    message("Stablize single cell RNAseq.")
    scRNA = stabilize_expr(scRNA, neighbor = correct.neighbor, n.core = n.core, npcs = length(dims))
  }
  
  genes_select = intersect(rownames(scRNA), rownames(spRNA))
  
  if(length(genes_select) < 10){
    stop("Too few intesected genes between scRNA and spRNA.")
  }
  
  # merge two objects
  message("integrating.")
  integrated = merge(scRNA, spRNA)
  
  # count level normalization
  message("normalization.")
  norm_data = expm1(integrated@assays$RNA@data)
  
  
  # trim high expression gene, which affect the mean of genes
  trim_quantil1 = Matrix::rowMeans(scRNA@assays$RNA@data[genes_select, ])
  trim_quantil1 = trim_quantil1[trim_quantil1 < quantile(trim_quantil1, prob = 0.98) ] 
  
  trim_quantil2 = Matrix::rowMeans(spRNA@assays$RNA@data[genes_select, ])
  trim_quantil2 = trim_quantil2[trim_quantil2 < quantile(trim_quantil2, prob = 0.98) ] 
  
  trim_genes = unique(c(names(trim_quantil1), names(trim_quantil2)))
  
  # normalize factor based on trimmed genes
  norm_factor = Matrix::colMeans(norm_data[trim_genes, ]) 
  norm_factor = as.numeric(norm_factor)
  norm_factor = norm_factor/mean(norm_factor)
  
  norm_data@x <- norm_data@x / rep.int(norm_factor, diff(norm_data@p))
  integrated@assays$RNA@data = log1p(norm_data)
  
  # remove cell with number of expressed gene < 98%
  integrated = integrated[, norm_factor != 0]
  
  rm(norm_data, norm_factor)
  gc()
  
  # check normalization
  #avg_expr = Seurat::AverageExpression(integrated, group.by="tech", slot="data")
  #boxplot(log1p(avg_expr$RNA[genes_select, ]))
  
  # global level integration
  
  integrated = suppressWarnings(run_harmony(integrated, genes_select, npc = length(dims)))
  
  # find neighbors
  integrated = Seurat::FindNeighbors(integrated, k.param = k.neighbor, reduction="harmony", dims=dims, return.neighbor = T)
  neigbors = sapply(colnames(subset(integrated, subset = tech == "scRNA", invert = TRUE)), function(i){
    Seurat::TopNeighbors(integrated@neighbors$RNA.nn, i, n = k.neighbor)
  })
  
  neigbors = as.data.frame(neigbors)
  
  integrated_scRNA = subset(integrated, subset = tech == "scRNA")
  integrated_merFISH = subset(integrated, subset = tech == "scRNA", invert = TRUE)
  cells_name = colnames(integrated_scRNA)
  
  neigbors = as.list(neigbors)
  
  # only select neighbor in scRNA data
  neigbors = lapply(neigbors, function(x){
    x = x[x %in% cells_name]
  })
  
  
  enhancer_expr = integrated_merFISH@assays$RNA@data
  
  # infer expression via scRNA
  message("infer expression.")
  if(weighted.KNN){
    # dynamics proportion to assign scRNA values to merFISH
    enhancer_expr = parallel::mclapply(colnames(enhancer_expr), function(cell){
      cell_neighbors = neigbors[[cell]]
      if (length(cell_neighbors) == 0){
        integrated@assays$RNA@data[, cell]
      }else{
        cor_dist = sparse.cor(integrated@assays$RNA@data[genes_select, c(cell, cell_neighbors)])[,1]
        cor_dist[is.na(cor_dist)] <- 0 
        cor_dist[cor_dist < 0] <- 0
        cor_dist = cor_dist ** 2
        # normalized correlation distance matrix
        cor_dist = cor_dist / sum(cor_dist)
        cor_dist = c((1-RNA.weight) * cor_dist[1], RNA.weight * cor_dist[-1]) # cor_dist[1] is the cell in spRNA, here = 1
        
        # inner produce
        infer_expr = integrated@assays$RNA@data[, c(cell, cell_neighbors)] %*% cor_dist 
        infer_expr[,1] 
      }
    }, mc.cores = n.core)
  }else{
    # constant proportion to assign scRNA values to merFISH
    enhancer_expr = parallel::mclapply(colnames(enhancer_expr), function(cell){
      cell_neighbors = neigbors[[cell]]
      if (length(cell_neighbors) == 0){
        enhancer_expr[, cell]
      }else if (length(cell_neighbors) == 1){
        (1-RNA.weight) * enhancer_expr[, cell] + RNA.weight * integrated_scRNA@assays$RNA@data[,cell_neighbors]
      }else{
        (1-RNA.weight) * enhancer_expr[, cell] + RNA.weight * Matrix::rowMeans(integrated_scRNA@assays$RNA@data[,cell_neighbors])
      }
    }, mc.cores = n.core)
  }
  
  enhancer_expr = do.call(cbind, enhancer_expr)
  colnames(enhancer_expr) = colnames(integrated_merFISH)
  
  # assign inferred expression value to new assay
  integrated_merFISH[[infered.assay]] = SeuratObject::CreateAssayObject(data = as(enhancer_expr, "dgCMatrix"))
  
  DefaultAssay(integrated_merFISH) <- infered.assay
  
  # add spatial information
  integrated_merFISH@images = spRNA@images

  for(img in names(integrated_merFISH@images)){
    integrated_merFISH@images[[img]]@coordinates <- integrated_merFISH@images[[img]]@coordinates[
        rownames(integrated_merFISH@images[[img]]@coordinates) %in% colnames(integrated_merFISH), ]
  }

  return(integrated_merFISH)
}


infer = function(
  spRNA,
  scRNA, 
  dims = 1:30,
  k.neighbor = 30,
  infered.assay = "enhanced",
  weighted.KNN = TRUE,
  RNA.weight = 0.5,
  n.core = 10,
  correct.spRNA = TRUE,
  correct.scRNA = FALSE,
  correct.neighbor = 3
){
  spRNA$tech = "spatial"
  scRNA$tech = "scRNA"
  
  if(is.null(spRNA@assays$RNA)){
    stop(paste(spRNA, " do not have 'RNA' assay."))
  }
  
  if(is.null(scRNA@assays$RNA)){
    stop(paste(scRNA, " do not have 'RNA' assay."))
  }
  
  if(length(spRNA@assays$RNA@data) == 0){
    stop(paste(spRNA, " is not normlized. Run Seurat::NormalizeData."))
  }
  
  if(length(scRNA@assays$RNA@data) == 0){
    stop(paste(scRNA, " is not normlized. Run Seurat::NormalizeData."))
  }
  
  if(correct.spRNA){
    message("Stablize spatial transcriptome.")
    spRNA = stabilize_expr(spRNA, neighbor = correct.neighbor, n.core = n.core, npcs = length(dims))
  }
  
  if(correct.scRNA){
    message("Stablize single cell RNAseq.")
    scRNA = stabilize_expr(scRNA, neighbor = correct.neighbor, n.core = n.core, npcs = length(dims))
  }
  
  genes_select = intersect(rownames(scRNA), rownames(spRNA))
  
  if(length(genes_select) < 10){
    stop("Too few intesected genes between scRNA and spRNA.")
  }
  
  # run pca
  SeuratObject::VariableFeatures(scRNA) = genes_select
  scRNA = Seurat::ScaleData(scRNA, verbose = FALSE)
  scRNA = Seurat::RunPCA(scRNA, npcs = length(dims), verbose = FALSE)
  
  SeuratObject::VariableFeatures(spRNA) = genes_select
  spRNA = Seurat::ScaleData(spRNA, verbose = FALSE)
  spRNA = Seurat::RunPCA(spRNA, npcs = length(dims), verbose = FALSE)
  
  
  # merge two objects
  message("1st level integration")
  anchors <- Seurat::FindIntegrationAnchors(object.list = list(spRNA, scRNA), normalization.method = "LogNormalize", 
                                            reduction = "rpca", anchor.features = genes_select,
                                            dims = dims,  k.anchor = 20, verbose = FALSE)
  integrated <- Seurat::IntegrateData(anchorset = anchors, dims = dims, 
                                      normalization.method = "LogNormalize", verbose = TRUE)
  rm(anchors)
  gc()

  
  # count level normalization
  message("normalization")
  norm_data = expm1(integrated@assays$RNA@data)
  
  
  # trim high expression gene, which affect the mean of genes
  trim_quantil1 = Matrix::rowMeans(scRNA@assays$RNA@data[genes_select, ])
  trim_quantil1 = trim_quantil1[trim_quantil1 < quantile(trim_quantil1, prob = 0.98) ] 
  
  trim_quantil2 = Matrix::rowMeans(spRNA@assays$RNA@data[genes_select, ])
  trim_quantil2 = trim_quantil2[trim_quantil2 < quantile(trim_quantil2, prob = 0.98) ] 
  
  trim_genes = unique(c(names(trim_quantil1), names(trim_quantil2)))
  
  # normalize factor based on trimmed genes
  norm_factor = Matrix::colMeans(norm_data[trim_genes, ]) 
  norm_factor = as.numeric(norm_factor)
  norm_factor = norm_factor/mean(norm_factor)
  
  norm_data@x <- norm_data@x / rep.int(norm_factor, diff(norm_data@p))
  integrated@assays$RNA@data = log1p(norm_data)
  
  # remove cell with number of expressed gene < 98%
  integrated = integrated[, norm_factor != 0]
  
  spRNA_images = spRNA@images
  
  rm(norm_data, norm_factor, scRNA, spRNA)
  gc()
  

  # check normalization
  #avg_expr = Seurat::AverageExpression(integrated, group.by="tech", slot="data")
  #boxplot(log1p(avg_expr$RNA[genes_select, ]))
  
  # 2nd level integration
  message("2nd level integration")
  SeuratObject::VariableFeatures(integrated, assay = "integrated") = genes_select
  integrated = Seurat::ScaleData(integrated, features = genes_select, vars.to.regress = "tech",
                                 assay = "integrated", verbose = FALSE)
  integrated = Seurat::RunPCA(integrated, npcs = length(dims), 
                              assay = "integrated", verbose = FALSE)

  integrated <- suppressWarnings( harmony::RunHarmony(
    object = integrated,
    group.by.vars = 'tech',
    plot_convergence = F,
    theta = 3,
    lambda = 0.5,
    assay.use = "integrated",
    verbose = FALSE,
    max.iter.harmony = 20
  ))
  
  # find neighbors
  integrated = Seurat::FindNeighbors(integrated, k.param = k.neighbor, reduction="harmony",
                                     dims=dims, return.neighbor = T, assay = "integrated")
  neigbors = sapply(colnames(subset(integrated, subset = tech == "scRNA", invert = TRUE)), function(i){
    Seurat::TopNeighbors(integrated@neighbors$integrated.nn, i, n = k.neighbor)
  })
  
  neigbors = as.data.frame(neigbors)
  
  DefaultAssay(integrated) <- "RNA"
  integrated = Seurat::DietSeurat(integrated, counts = F, assays = "RNA")
  integrated_scRNA = suppressWarnings(subset(integrated, subset = tech == "scRNA"))
  integrated_merFISH = suppressWarnings(subset(integrated, subset = tech == "scRNA", invert = TRUE))
  cells_name = colnames(integrated_scRNA)
  
  neigbors = as.list(neigbors)
  
  # only select neighbor in scRNA data
  neigbors = lapply(neigbors, function(x){
    x = x[x %in% cells_name]
  })
  
  
  enhancer_expr = integrated_merFISH@assays$RNA@data
  
  # infer expression via scRNA
  message("infer expression.")
  if(weighted.KNN){
    # dynamics proportion to assign scRNA values to merFISH
    enhancer_expr = parallel::mclapply(colnames(enhancer_expr), function(cell){
      cell_neighbors = neigbors[[cell]]
      if (length(cell_neighbors) == 0){
        integrated@assays$RNA@data[, cell]
      }else{
        cor_dist = sparse.cor(integrated@assays$RNA@data[genes_select, c(cell, cell_neighbors)])[,1]
        cor_dist[is.na(cor_dist)] <- 0 
        cor_dist[cor_dist < 0] <- 0
        cor_dist = cor_dist ** 2
        # normalized correlation distance matrix
        cor_dist = cor_dist / sum(cor_dist)
        cor_dist = c((1-RNA.weight) * cor_dist[1], RNA.weight * cor_dist[-1]) # cor_dist[1] is the cell in spRNA, here = 1
        
        # inner produce
        infer_expr = integrated@assays$RNA@data[, c(cell, cell_neighbors)] %*% cor_dist 
        infer_expr[,1] 
      }
    }, mc.cores = n.core)
  }else{
    # constant proportion to assign scRNA values to merFISH
    enhancer_expr = parallel::mclapply(colnames(enhancer_expr), function(cell){
      cell_neighbors = neigbors[[cell]]
      if (length(cell_neighbors) == 0){
        enhancer_expr[, cell]
      }else if (length(cell_neighbors) == 1){
        (1-RNA.weight) * enhancer_expr[, cell] + RNA.weight * integrated_scRNA@assays$RNA@data[,cell_neighbors]
      }else{
        (1-RNA.weight) * enhancer_expr[, cell] + RNA.weight * Matrix::rowMeans(integrated_scRNA@assays$RNA@data[,cell_neighbors])
      }
    }, mc.cores = n.core)
  }
  
  enhancer_expr = do.call(cbind, enhancer_expr)
  colnames(enhancer_expr) = colnames(integrated_merFISH)
  
  # assign inferred expression value to new assay
  integrated_merFISH[[infered.assay]] = SeuratObject::CreateAssayObject(data = as(enhancer_expr, "dgCMatrix"))
  
  DefaultAssay(integrated_merFISH) <- infered.assay
  
  integrated_merFISH = Seurat::DietSeurat(integrated_merFISH, assays = infered.assay)
  
  # add spatial information
  integrated_merFISH@images = spRNA_images
  
  for(img in names(integrated_merFISH@images)){
    integrated_merFISH@images[[img]]@coordinates <- integrated_merFISH@images[[img]]@coordinates[
      rownames(integrated_merFISH@images[[img]]@coordinates) %in% colnames(integrated_merFISH), ]
  }
  
  return(integrated_merFISH)
}

#' Hierarchy version of iSpatial

infer_Hierarchy = function(
  spRNA,
  scRNA, 
  dims = 1:30,
  k.neighbor = 30,
  cluster.resolution = 0.1,
  infered.assay = "enhanced",
  weighted.KNN = TRUE,
  RNA.weight = 0.3,
  n.core = 10,
  correct.spRNA = TRUE,
  correct.scRNA = FALSE,
  correct.neighbor = 3
){
  spRNA$tech = "spatial"
  scRNA$tech = "scRNA"
  
  if(is.null(spRNA@assays$RNA)){
    stop(paste(spRNA, " do not have 'RNA' assay."))
  }
  
  if(is.null(scRNA@assays$RNA)){
    stop(paste(scRNA, " do not have 'RNA' assay."))
  }
  
  if(length(spRNA@assays$RNA@data) == 0){
    stop(paste(spRNA, " is not normlized. Run Seurat::NormalizeData."))
  }
  
  if(length(scRNA@assays$RNA@data) == 0){
    stop(paste(scRNA, " is not normlized. Run Seurat::NormalizeData."))
  }
  
  if(correct.spRNA){
    message("Stablize spatial transcriptome.")
    spRNA = stabilize_expr(spRNA, neighbor = correct.neighbor, n.core = n.core, npcs = length(dims))
  }
  
  if(correct.scRNA){
    message("Stablize single cell RNAseq.")
    scRNA = stabilize_expr(scRNA, neighbor = correct.neighbor, n.core = n.core, npcs = length(dims))
  }
  
  genes_select = intersect(rownames(scRNA), rownames(spRNA))
  
  if(length(genes_select) < 10){
    stop("Too few intesected genes between scRNA and spRNA.")
  }
  
  # merge two objects
  message("integrating.")
  integrated = merge(scRNA, spRNA)
  
  # count level normalization
  message("normalization.")
  norm_data = expm1(integrated@assays$RNA@data)
  trim_quantil = Matrix::rowMeans(norm_data[genes_select, ])
  
  # trim high expression gene, which affect the mean of genes
  trim_quantil = trim_quantil[trim_quantil < quantile(trim_quantil, prob = 0.95) ] 
  
  # normalize factor based on trimmed genes
  norm_factor = Matrix::colMeans(norm_data[names(trim_quantil), ]) 
  norm_factor = as.numeric(norm_factor)
  norm_factor = norm_factor/mean(norm_factor)
  
  norm_data@x <- norm_data@x / rep.int(norm_factor, diff(norm_data@p))
  integrated@assays$RNA@data = log1p(norm_data)
  
  # remove cell with number of expressed gene < 95%
  integrated = integrated[, norm_factor != 0]
  
  rm(norm_data, norm_factor)
  gc()
  
  # check normalization
  #avg_expr = Seurat::AverageExpression(integrated, group.by="tech", slot="data")
  #boxplot(log1p(avg_expr$RNA[genes_select, ]))
  
  # global level integration
  
  integrated = suppressWarnings(run_harmony(integrated, genes_select, npc = length(dims)))
  

  # cluster
  integrated <- Seurat::FindNeighbors(integrated, reduction="harmony")
  integrated <- Seurat::FindClusters(integrated, resolution = cluster.resolution)
  

  integrated_cluster = Seurat::SplitObject(integrated, split.by = "seurat_clusters")
  

    #doParallel::registerDoParallel(n.core)
    #neigbors = foreach::foreach(obj = integrated_cluster, .combine = cbind) foreach::`%dopar%` {
    #  # cluster level integration
    #  obj = suppressWarnings(run_harmony(obj, genes_select))
    #  # find neighbor
    #  obj = Seurat::FindNeighbors(obj, k.param = k.neighbor, reduction="harmony", dims=dims, return.neighbor = T)
    #  res = sapply(colnames(subset(obj, subset = tech == "scRNA", invert = TRUE)), function(i){
    #    Seurat::TopNeighbors(obj@neighbors$RNA.nn, i, n = k.neighbor)
    #  })
    #  res
    #}
    #neigbors = as.data.frame(neigbors)
    #doParallel::stopImplicitCluster()  }
  
  
  neigbors = parallel::mclapply(integrated_cluster, function(obj){
    # cluster level integration
    obj = suppressWarnings(run_harmony(obj, genes_select, npc = length(dims)))
    # find neighbor
    obj = Seurat::FindNeighbors(obj, k.param = k.neighbor, reduction="harmony", dims=dims, return.neighbor = T)
    res = sapply(colnames(subset(obj, subset = tech == "scRNA", invert = TRUE)), function(i){
      Seurat::TopNeighbors(obj@neighbors$RNA.nn, i, n = k.neighbor)
    })
    res
  }, mc.cores = n.core)
  neigbors = do.call(cbind, neigbors)
  neigbors = as.data.frame(neigbors)
  
  
  integrated_scRNA = subset(integrated, subset = tech == "scRNA")
  integrated_merFISH = subset(integrated, subset = tech == "scRNA", invert = TRUE)
  cells_name = colnames(integrated_scRNA)
  
  neigbors = as.list(neigbors)
  neigbors = lapply(neigbors, function(x){
    x = x[x %in% cells_name]
  })
  
  enhancer_expr = integrated_merFISH@assays$RNA@data
  
  # infer expression via scRNA
  message("infer expression.")
  if(weighted.KNN){
    # dynamics proportion to assign scRNA values to merFISH
    enhancer_expr = parallel::mclapply(colnames(enhancer_expr), function(cell){
      cell_neighbors = neigbors[[cell]]
      if (length(cell_neighbors) == 0){
        integrated@assays$RNA@data[, cell]
      }else{
        cor_dist = sparse.cor(integrated@assays$RNA@data[genes_select, c(cell, cell_neighbors)])[,1]
        cor_dist[is.na(cor_dist)] <- 0 
        cor_dist[cor_dist < 0] <- 0
        cor_dist = cor_dist ** 2
        # normalized correlation distance matrix
        cor_dist = cor_dist / sum(cor_dist)
        cor_dist = c((1-RNA.weight) * cor_dist[1], RNA.weight * cor_dist[-1]) # cor_dist[1] is the cell in spRNA, here = 1
        
        # inner produce
        infer_expr = integrated@assays$RNA@data[, c(cell, cell_neighbors)] %*% cor_dist 
        infer_expr[,1] 
      }
    }, mc.cores = n.core)
  }else{
    # constant proportion to assign scRNA values to merFISH
    enhancer_expr = parallel::mclapply(colnames(enhancer_expr), function(cell){
      cell_neighbors = neigbors[[cell]]
      if (length(cell_neighbors) == 0){
        enhancer_expr[, cell]
      }else if (length(cell_neighbors) == 1){
        (1-RNA.weight) * enhancer_expr[, cell] + RNA.weight * integrated_scRNA@assays$RNA@data[,cell_neighbors]
      }else{
        (1-RNA.weight) * enhancer_expr[, cell] + RNA.weight * Matrix::rowMeans(integrated_scRNA@assays$RNA@data[,cell_neighbors])
      }
    }, mc.cores = n.core)
  }
  
  enhancer_expr = do.call(cbind, enhancer_expr)
  colnames(enhancer_expr) = colnames(integrated_merFISH)
  
  # assign inferred expression value to new assay
  integrated_merFISH[[infered.assay]] = SeuratObject::CreateAssayObject(data = as(enhancer_expr, "dgCMatrix"))
  
  integrated_merFISH = Seurat::DietSeurat(integrated_merFISH, assays = infered.assay)
  
  DefaultAssay(integrated_merFISH) <- infered.assay
  
  # add spatial information
  integrated_merFISH@images = spRNA@images
  
  for(img in names(integrated_merFISH@images)){
    integrated_merFISH@images[[img]]@coordinates <- integrated_merFISH@images[[img]]@coordinates[
      rownames(integrated_merFISH@images[[img]]@coordinates) %in% colnames(integrated_merFISH), ]
  }
  
  return(integrated_merFISH)
  
}

