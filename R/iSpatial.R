

############# infer whole genome spatial transcriptome from single cell RNAseq

#' run harmony integration
#' 
#' This function calls harmony package(https://github.com/immunogenomics/harmony) 
#' to integrate the scRNA-seq and spatial transcriptome data.
#' 
#' @param seurat_obj merged seurat object
#' @param genes_select genes used to integration, usually targeted genes.
#' @return seurat object
#' 
run_harmony = function(seurat_obj, genes_select){
  VariableFeatures(seurat_obj) = genes_select
  seurat_obj = ScaleData(seurat_obj, features = genes_select, verbose = FALSE)
  seurat_obj = RunPCA(seurat_obj, npcs = 30, features = genes_select, verbose = FALSE)
  
  # harmony
  seurat_obj <- harmony::RunHarmony(
    object = seurat_obj,
    group.by.vars = 'tech',
    plot_convergence = F,
    theta = 3,
    lambda = 0.1,
    verbose = FALSE,
    max.iter.harmony = 15
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
stabilize_expr = function(obj, neighbor = 5, npcs = 10, n.core = 10){
  obj = FindVariableFeatures(obj, selection.method = "vst", verbose = FALSE) 
  obj = ScaleData(obj, verbose = FALSE)
  obj = RunPCA(obj, npcs = npcs, verbose = FALSE)
  obj <- FindNeighbors(obj, return.neighbor = T, k.param = neighbor, verbose = FALSE)
  
  # imputation expr value by KNN
  enhancer_expr = obj@assays$RNA@data
  enhancer_expr = mclapply(colnames(enhancer_expr), function(cell){
    cell_neighbors = TopNeighbors(obj@neighbors$RNA.nn, cell, n=neighbor)
    0.7*rowMedians(obj@assays$RNA@data[,cell_neighbors[-1]]) + 0.3 * obj@assays$RNA@data[,cell_neighbors[1]]
    #Matrix::rowMeans(obj@assays$RNA@data[,cell_neighbors])
  }, mc.cores = n.core)
  enhancer_expr = do.call(cbind, enhancer_expr)
  colnames(enhancer_expr) = colnames(obj)
  rownames(enhancer_expr) = rownames(obj)
  obj@assays$RNA@data = as(enhancer_expr, "dgCMatrix")
  rm(enhancer_expr)
  gc()
  return(obj)
}

#' calculate correlation
#' 
#' @param x sparse matrix
#' @return correlation matrix
#' 
sparse.cor <- function(x){
  n <- nrow(x)
  m <- ncol(x)
  ii <- unique(x@i)+1 # rows with a non-zero element
  
  Ex <- colMeans(x)
  nozero <- as.vector(x[ii,]) - rep(Ex,each=length(ii)) 
  
  covmat <- ( crossprod(matrix(nozero,ncol=m)) +
                crossprod(t(Ex))*(n-length(ii))
  )/(n-1)
  sdvec <- sqrt(diag(covmat))
  covmat/crossprod(t(sdvec))
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
#' @param RNA.weight the weight of scRNA-seq expression when genes detected both
#' in spRNA and scRNA. When choose defaulted "NA", we use a dynamics weight
#' to assign scRNA values to spRNA. The dynamics weight based on the correlation
#' between scRNA cell and spRNA cell. If a RNA.weight (from 0 to 1) were given,
#' the inferred expression = (1 - RNA.weight) \* spRNA + RNA.weight \* scRNA.
#' @param n.core number of CPU cores used to parallel.
#' @param correct.expr Whether to stabilize expression in scRNA and spRNA.
#' @param correct.neighbor number of nearest neighbors used to correct expr.
#' 
#' @return returns a seurat object with inferred expression in infered.assay.
#' 
#' @export
#' 

iSpatial = function(
  spRNA,
  scRNA, 
  dims = 1:30,
  k.neighbor = 30,
  infered.assay = "enhanced",
  RNA.weight = NA,
  n.core = 10,
  correct.expr = TRUE,
  correct.neighbor = 5
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
  
  
  if(correct.expr){
    message("Stablize spatial transcriptome.")
    spRNA = stabilize_expr(spRNA, neighbor = correct.neighbor, n.core = n.core, npcs = length(dims))
    message("Stablize single cell RNAseq.")
    scRNA = stabilize_expr(scRNA, neighbor = correct.neighbor, n.core = n.core, npcs = length(dims))
  }
  
  genes_select = intersect(rownames(scRNA), rownames(spRNA))
  
  if(length(genes_select) < 10){
    stop("Too few intesected genes between scRNA and spRNA.")
  }
  
  
  scRNA = NA_scRNA
  spRNA = NA_merFISH
  
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
  #avg_expr = AverageExpression(integrated, group.by="tech", slot="data")
  #boxplot(log1p(avg_expr$RNA[genes_select, ]))
  
  # global level integration
  integrated = run_harmony(integrated, genes_select)
  
  # find neighbors
  integrated = FindNeighbors(integrated, k.param = k.neighbor, reduction="harmony", dims=dims, return.neighbor = T)
  neigbors = sapply(colnames(subset(integrated, subset = tech == "scRNA", invert = TRUE)), function(i){
    TopNeighbors(integrated@neighbors$RNA.nn, i, n = k.neighbor)
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
  if(is.na(RNA.weight)){
    # dynamics proportion to assign scRNA values to merFISH
    enhancer_expr = parallel::mclapply(colnames(enhancer_expr), function(cell){
      cell_neighbors = neigbors[[cell]]
      if (length(cell_neighbors) == 0){
        integrated@assays$RNA@data[, cell]
      }else{
        cor_dist = sparse.cor(integrated@assays$RNA@data[genes_select, c(cell, cell_neighbors)])[,1]
        cor_dist[is.na(cor_dist)] <- 0 
        cor_dist[cor_dist < 0] <- 0
        cor_dist = cor_dist**3
        # normalized correlation distance matrix
        cor_dist = cor_dist/sum(cor_dist) 
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
  integrated_merFISH[[infered.assay]] = CreateAssayObject(data = as(enhancer_expr, "dgCMatrix"))
  
  DefaultAssay(integrated_merFISH) <- infered.assay
  
  # add spatial information
  coord.df = spRNA@images$image@coordinates
  integrated_merFISH@images$image =  new(
    Class = 'SlideSeq',
    assay = infered.assay,
    key = "image_",
    coordinates = coord.df
  )
  
  # remove NA gene
  #integrated_merFISH = integrated_merFISH[!is.na(rowSums(integrated_merFISH[[infered.assay]]@data)),]
  
  return(integrated_merFISH)
}

#' Hierarchy version of iSpatial

iSpatial_Hierarchy = function(
  spRNA,
  scRNA, 
  dims = 1:30,
  n.core = 10,
  RNA.weight = 0.3,
  cluster.resolution = 0.1,
  infered.assay = "enhanced"
){
  spRNA$tech = "spatial"
  scRNA$tech = "scRNA"
  
  genes_select = intersect(rownames(scRNA), rownames(spRNA))
  
  integrated = merge(scRNA, spRNA)
  
  # count level normalization
  norm_data = integrated@assays$RNA@counts
  trim_quantil = rowMeans(norm_data[genes_select, ])
  trim_quantil = trim_quantil[trim_quantil < quantile(trim_quantil, prob = 0.95) ] # trim high expression gene, which affect the mean of genes
  
  norm_factor = Matrix::colMeans(norm_data[names(trim_quantil), ]) # based on trimmed genes
  norm_factor = as.numeric(norm_factor)
  
  norm_data@x <- norm_data@x / rep.int(norm_factor, diff(norm_data@p))
  integrated@assays$RNA@data = log1p(norm_data)
  
  # remove cell with number of expressed gene < 95% total
  integrated = integrated[, norm_factor != 0]
  
  rm(norm_data, norm_factor)
  gc()
  
  # check normalization
  #avg_expr = AverageExpression(integrated, group.by="tech", slot="data")
  #boxplot(log1p(avg_expr$RNA[genes_select, ]))
  
  # global level integration
  integrated = run_harmony(integrated, genes_select)
  
  # cluster
  integrated <- FindNeighbors(integrated, reduction="harmony")
  integrated <- FindClusters(integrated, resolution = cluster.resolution)
  
  #DimPlot(integrated, reduction="umap", group.by = "seurat_clusters", raster=FALSE) & NoAxes() 
  
  integrated_cluster = SplitObject(integrated, split.by = "seurat_clusters")
  

  doParallel::registerDoParallel(n.core)
  
  neigbors = foreach::foreach(obj = integrated_cluster, .combine = cbind) %dopar%  {
    # cluster level integration
    obj = run_harmony(obj, genes_select)
    # find neighbor
    obj = FindNeighbors(obj, k.param = 20, reduction="harmony", dims=dims, return.neighbor = T)
    res = sapply(colnames(subset(obj, subset = tech == "scRNA", invert = TRUE)), function(i){
      TopNeighbors(obj@neighbors$RNA.nn, i, n=20)
    })
    res
  }
  neigbors = as.data.frame(neigbors)
  doParallel::stopImplicitCluster()
  
  integrated_scRNA = subset(integrated, subset = tech == "scRNA")
  integrated_merFISH = subset(integrated, subset = tech == "scRNA", invert = TRUE)
  cells_name = colnames(integrated_scRNA)
  
  neigbors = as.list(neigbors)
  neigbors = lapply(neigbors, function(x){
    x = x[x %in% cells_name]
  })
  
  enhancer_expr = integrated_merFISH@assays$RNA@data
  
  # infer expression via scRNA
  enhancer_expr = parallel::mclapply(colnames(enhancer_expr), function(cell){
    cell_neighbors = neigbors[[cell]]
    if (length(cell_neighbors) == 0){
      enhancer_expr[, cell]
    }else if (length(cell_neighbors) == 1){
      (1-RNA.weight) * enhancer_expr[, cell] + RNA.weight * integrated_scRNA@assays$RNA@data[,cell_neighbors]
    }else{
      (1-RNA.weight) * enhancer_expr[, cell] + RNA.weight * rowMeans(integrated_scRNA@assays$RNA@data[,cell_neighbors])
    }
  }, mc.cores = n.core)
  enhancer_expr = do.call(cbind, enhancer_expr)
  colnames(enhancer_expr) = colnames(integrated_merFISH)
  
  # assign inferred expression value to new assay
  integrated_merFISH[[infered.assay]] = CreateAssayObject(data = as(enhancer_expr, "dgCMatrix"))
  
  DefaultAssay(integrated_merFISH) <- infered.assay
  
  # add spatial information
  coord.df = spRNA@images$image@coordinates
  integrated_merFISH@images$image =  new(
    Class = 'SlideSeq',
    assay = infered.assay,
    key = "image_",
    coordinates = coord.df
  )
  
  return(integrated_merFISH)
}

