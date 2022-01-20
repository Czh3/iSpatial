
#' Calculate spatially expression correlation
#'
#' Calculate spatially expression correlation between two objects. This function
#' can be used to validation the performance of inferred data.
#'
#' @param obj1 seurat object 1
#' @param assay1 choose assay for obj1
#' @param obj2 seurat object 2
#' @param assay2 choose assay for obj2
#' @param feature select features to calulate 
#' @param bin number of bins used to segregate the image
#' @param n.sample how many times for sampling regions
#'
#' @return Spatial expression correlation of each gene 
#'
#' @export
#'
spatial_cor = function(obj1, 
                       assay1, 
                       obj2, 
                       assay2, 
                       feature = NA,
                       bin = 20,
                       n.core = 10,
                       n.sample = 300){

  if(is.na(feature)){
    feature = Seurat::AverageExpression(obj1)[[1]]
    feature = names(feature[feature[,1] > 0, ])
  }
  
  cor_bin = parallel::mclapply(1:300, function(x) {
    cells = sample_region(obj1)
    if(length(cells) > 2){
      cbind(
        Matrix::rowMeans((obj1[[assay1]]@data[feature, cells]), na.rm = T),
        Matrix::rowMeans((obj2[[assay2]]@data[feature, cells]), na.rm = T)
      )
    }
  }, mc.cores = n.core)
  
  cor_bin = do.call(cbind, cor_bin)
  
  apply(cor_bin, 1, function(x){
    cor(x[seq(1, length(x), by=2)],
        x[seq(2, length(x), by=2)])
  })
}


#' Random sample spatial regions
#'
#' Get the cell names from random sampling spatial regions
#'
#' @param obj seurat object
#' @param bin number of bins used to segregate the image
#'
#' @return A vector including the cell names in a random sampling region
#'
#'

sample_region = function(obj, bin = 20){
  
  bin = 1/bin
  
  image = sample(names(obj@images))
  image.coordinates = obj@images[[image]]@coordinates
  x = sample(image.coordinates$x, 1)
  y = sample(image.coordinates$y, 1)
  
  n_bin_x = (max(image.coordinates$x) - min(image.coordinates$x)) * bin
  n_bin_y = (max(image.coordinates$y) - min(image.coordinates$y)) * bin
  cell_in = rownames(image.coordinates[image.coordinates$x > x & image.coordinates$x < x+n_bin_x  & image.coordinates$y > y & image.coordinates$y < y+n_bin_y,])
  cell_in
}



