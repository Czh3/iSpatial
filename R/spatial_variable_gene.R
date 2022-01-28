#' @import Seurat
#' @import SeuratObject
#' @import ggplot2
#'
NULL

#' Find spatial variable genes
#' 
#' @param obj seurat object
#' @param assay select assay
#' @param img select image 
#' @param bin number of bins used to segregate the image
#' @param n.core number of CPU cores used
#' 
#' @return return a matrix wtih spatial variable genes
#' 
spatial_var_gene_slice <- function(obj, assay, img, bin, n.core = 1){
  coord = obj@images[[img]]@coordinates
  
  obj = obj[, rownames(coord)]
  
  bin = as.numeric(bin)
  bin_x = (max(coord[,1]) - min(coord[,1])) / bin
  bin_y = (max(coord[,2]) - min(coord[,2])) / bin
  obj@meta.data$x = floor(coord[,1]/bin_x)
  obj@meta.data$y = floor(coord[,2]/bin_y)
  obj@meta.data$xy = paste(obj@meta.data$x, obj@meta.data$y)
  obj_expr = Seurat::AverageExpression(obj, group.by = "xy", assays = assay)[[1]]
  obj_expr = obj_expr[Matrix::rowSums(obj_expr) != 0, ]
  
  obj_random = obj[rownames(obj_expr), ]
  set.seed(3)
  obj_random@meta.data$xy = sample(obj@meta.data$xy)
  
  obj_random_expr = Seurat::AverageExpression(obj_random, group.by = "xy", assays = assay)[[1]]
  
  SVGs = parallel::mclapply(rownames(obj_expr), function(gene){
    p = wilcox.test(obj_expr[gene, ], obj_random_expr[gene, ])$p.value
    
    obj_expr_sort = sort(obj_expr[gene, ])
    obj_random_expr_sort = sort(obj_random_expr[gene, ])
    
    n_bins = bin*bin * 0.05
    
    #obs_exp = log((mean(obj_expr_sort[(length(obj_expr_sort)-n_bins):length(obj_expr_sort)]) - mean(obj_expr_sort[1:n_bins]))/
    #  (mean(obj_random_expr_sort[(length(obj_random_expr_sort)-n_bins):length(obj_random_expr_sort)]) - mean(obj_random_expr_sort[1:n_bins])))
    #spatial_bias = mean(obj_expr_sort[(length(obj_expr_sort)-n_bins):length(obj_expr_sort)]) - mean(obj_expr_sort[1:n_bins])
    c(gene, img, p)
  }, mc.cores = n.core)
  SVGs = do.call(rbind, SVGs)
  SVGs = as.data.frame(SVGs)
  
  # percent cell expressed
  pct_expr = Matrix::rowSums(obj[[assay]]@data != 0)
  pct_expr = pct_expr/ncol(obj[[assay]]@data)
  
  SVGs$pct_expr = as.numeric(pct_expr[SVGs$V1])
  
  colnames(SVGs) = c("gene", "image", "p.value", "pct.expr")
  SVGs$padj = p.adjust(SVGs$p.value, method = "BH")
  SVGs
}

#' Find spatial variable genes
#' 
#' Bin the images into small chunks. For each gene, detected the numbers of 
#' expression cells in each chunk, comparing the distribution with random. 
#' 
#' @param obj seurat object
#' @param assay select assay
#' @param bin number of bins used to segregate the image
#' @param n.core number of CPU cores used
#' 
#' @return return a matrix wtih spatial variable genes
#' 
#' @export
#' 

spatial_variable_genes = function(obj, assay, bin = 20, n.core = 1){
  if(length(obj@images) == 0){
    stop("Check your spatial seurat object, lost image information.")
  }
  
  images = names(obj@images)
  SVGs = parallel::mclapply(images, function(img) spatial_var_gene_slice(obj, assay, img, bin, n.core)
  , mc.cores = n.core)
  SVGs = do.call(rbind, SVGs)
  SVGs = as.data.frame(SVGs)
  SVGs$pct.expr = as.numeric(SVGs$pct.expr)
  SVGs$p.value = as.numeric(SVGs$p.value)
  rownames(SVGs) = 1:nrow(SVGs)
  SVGs
}

#' get spatial expression features in each image
#' 
#' @param obj seurat object
#' @param assay assay to use
#' @param feature features (genes) to get spatial expression patterns
#' @param img choose image
#' @param bin number of bins used to segregate the image
#' @param n.core number of CPU cores used
#' 
#' @return return a matrix. row: gene, column: spatial expression features
#' 

get_spatial_features_image <- function(obj, assay, feature, img, bin, n.core = 1){
  coord = obj@images[[img]]@coordinates
  
  obj = obj[, rownames(coord)]
  
  if(is.na(feature)){
    feature = rownames(obj)[Matrix::rowSums(obj[[assay]]@data) > 0]
  }
  
  bin = as.numeric(bin)
  bin_x = (max(coord[,1]) - min(coord[,1])) / bin
  bin_y = (max(coord[,2]) - min(coord[,2])) / bin
  obj@meta.data$x = floor(coord[,1]/bin_x)
  obj@meta.data$y = floor(coord[,2]/bin_y)
  obj@meta.data$xy = paste(img, obj@meta.data$x, obj@meta.data$y)
  obj_expr = Seurat::AverageExpression(obj, group.by = "xy", assays = assay, features = feature)[[1]]
  obj_expr = obj_expr[Matrix::rowSums(obj_expr) != 0, ]
  
  obj_expr
}

#' Get spatial expression features 
#' 
#' @param obj seurat object
#' @param assay assay to use
#' @param feature features (genes) to get spatial expression patterns
#' @param bin number of bins used to segregate the image
#' @param n.core number of CPU cores used
#' 
#' @return return a matrix. row: gene, column: spatial expression features
#' 
#' @export
#' 

get_spatial_features <- function(obj, assay, feature, bin, n.core = 1){
  spatial_feature = lapply(names(obj@images), function(img) {
    get_spatial_features_image(obj, assay, feature, img, bin, n.core = n.core)
    })
  spatial_feature = do.call(cbind, spatial_feature)
  
  spatial_feature
}


#' Cluster spatial expression features 
#' 
#' Clustering the genes based on the spatial expression pattern.
#' 
#' @param obj seurat object
#' @param assay assay to use
#' @param feature features (genes) to get spatial expression patterns
#' @param bin number of bins used to segregate the image
#' @param n.cluster number of clusters
#' @param cluster.method hclust method
#' @param cor.method correlation method
#' @param ncol.plot number of columns for plotting
#' @param n.core number of CPU cores used
#' 
#' @return return a vector showing the genes in each cluster
#' 
#' @export
#' 

cluster_spatial_expr_pattern <- function(obj, 
                                         assay = "enhanced",
                                         feature, 
                                         bin = 20, 
                                         n.cluster = 12, 
                                         cluster.method = "ward.D2",
                                         cor.method = "pearson",
                                         ncol.plot = 4, 
                                         n.core = 1){
  spatial_feature = get_spatial_features(obj, assay, feature, bin, n.core = n.core)
  
  SVG_cor = cor(t(spatial_feature), method = cor.method)
  
  h = hclust(as.dist(1-SVG_cor), method = cluster.method)
  
  SVG_cluster = cutree(h, k=n.cluster)
  message("Number of genes in each cluster:")
  message(paste0(capture.output(table(SVG_cluster)), collapse = "\n"))
  
  plots = lapply(1:n.cluster, function(x){
      suppressMessages(spatial_signature_plot(obj, assay, names(SVG_cluster[SVG_cluster==x]), paste("C", as.character(x))))
  })
  
  gridExtra::grid.arrange(grobs = plots, ncol = ncol.plot)
  
  return(SVG_cluster)
}

#' Plot spatial expression
#' 
#' Plot the spatial expression patterns on a group of genes.
#'
#' @param object seurat object
#' @param assay choose assay
#' @param gene.set vector contain the genes to plot
#' @param set.name set the name for the plot
#' @param pt.size.factor point size
#' @param alpha transparency of the points
#' @param color.scale set the color scale. Default:
#' c("gray95", "gray80", "pink", "orange", "red", "red4")
#' @param quantile.cutoff 0 to 1.
#' Set the quantile of expression of the max value. 
#' Default:  0.95% of the max.
#' @param legend whether to plot the legend
#' 
#' @return ggplot object
#' 
#' @export
#'

spatial_signature_plot = function(object, 
                                  assay = "enhanced", 
                                  gene.set, 
                                  set.name, 
                                  pt.size.factor = 2,
                                  alpha = 0.5,
                                  color.scale = c("gray95", "gray80", "wheat", "pink", "orange", "red", "red4"),
                                  quantile.cutoff = 0.95, 
                                  legend = FALSE){
  
  gene.set = gene.set[gene.set %in% rownames(object[[assay]]@data)]
  if(length(gene.set) == 1){
    mean.exp <- object[[assay]]@data[gene.set, ]
  } else{
    mean.exp <- Matrix::colMeans(x = (object[[assay]]@data[gene.set, ]), na.rm = TRUE)
  }
  
  gene.set.name = set.name
  object@meta.data[gene.set.name] <- mean.exp
  
  
  
  max.cutoff = quantile(mean.exp, prob = quantile.cutoff, na.rm=T)
  
  
  if(legend){
    SpatialFeaturePlot(object, features = gene.set.name, pt.size.factor=pt.size.factor, stroke=NA,  alpha = alpha,
                       slot="data", max.cutoff = max.cutoff) +
      theme(legend.position = "right") +
      scale_fill_gradientn(colours = color.scale) 
  }else{
    SpatialFeaturePlot(object, features = gene.set.name, pt.size.factor=pt.size.factor, stroke=NA,  alpha = alpha,
                       slot="data", max.cutoff = max.cutoff) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(set.name)+
      scale_fill_gradientn(colours = color.scale) 
  }
}
