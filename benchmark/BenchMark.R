# benchmark
library(ggplot2)
library(Seurat)
library(iSpatial)

Sys.setenv(RETICULATE_PYTHON = "/netscr/chaozhang/miniconda3/bin/python")
library(reticulate)

options(future.globals.maxSize = 30000 * 1024^2)

iSpatial_PCC = function(dataset){
  print(dataset)
  
  scRNA = data.table::fread(paste0("../data/benchmark/", dataset,"/scRNA_count.txt"), data.table = F)
  rownames(scRNA) = scRNA$V1
  scRNA = scRNA[, -1]
  
  scRNA = CreateSeuratObject(counts = scRNA, project = "scRNA",
                             assay = "RNA")
  
  ST = data.table::fread(paste0("../data/benchmark/", dataset,"/Insitu_count.txt"), data.table = F)
  rownames(ST) = paste0("ST", rownames(ST))
  ST = CreateSeuratObject(counts = t(ST), project = "ST",
                          assay = "RNA")
  
  coord = read.table(paste0("../data/benchmark/", dataset,"/Locations.txt"), header = T)
  colnames(coord) = c("x", "y")
  rownames(coord) = colnames(ST)
  
  ST@images$image =  new(
    Class = 'SlideSeq',
    assay = "RNA",
    key = "image_",
    coordinates = coord
  )
  
  ST = NormalizeData(ST)
  
  #10x CV
  np <- import("numpy")
  train_list <- np$load(paste0("../data/benchmark/", dataset,"/train_list.npy"), allow_pickle = TRUE)
  test_list <- np$load(paste0("../data/benchmark/", dataset,"/test_list.npy"), allow_pickle = TRUE)
  
  trainshape = dim(train_list)
  testshape = dim(test_list)
  if (length(trainshape) == 2){
    train_list <- array(train_list, dim = c(trainshape[2],trainshape[1]))
    test_list <- array(test_list, dim = c(testshape[2],testshape[1]))
  }
  if (length(trainshape) == 1){
    train_list <- do.call(cbind, train_list)
    test_list <- do.call(cbind, test_list)
  }
  
  # To speed up, only focus on ST genes.
  scRNA = scRNA[unique(c(train_list)), ]
  scRNA = scRNA[, scRNA$nFeature_RNA > 10]
  scRNA = NormalizeData(scRNA)
  
  iSpatial_genecor = c()
  ##
  for (i in 1:10){
    print(paste0("10X CV: ", i))
    train_gene = train_list[,i]
    validation_gene = test_list[,i]
    
    ST_train = ST[train_gene, ]
    #
    if(length(train_gene) > 50){
      DM = 50
    } else {
      DM = floor(length(train_gene)/2)
    }
    
    iSpatial_group = infer_v2(ST_train,
                             scRNA, 
                             k.neighbor = 30,
                             dims = 1:DM,
                             RNA.weight = 1)
    
    infer_expr = iSpatial_group@assays$enhanced@data[validation_gene, ]
    raw_expr = ST@assays$RNA@counts[validation_gene, colnames(iSpatial_group)]
    
    iSpatial_c = cor(t(as.matrix(expm1(infer_expr))), t(as.matrix((raw_expr))))
    iSpatial_genecor = c(iSpatial_genecor, diag(iSpatial_c))
  }
  
  write.table(as.matrix(iSpatial_genecor),
              paste0("../data/benchmark/", dataset,"/iSpatial_cor.txt"), quote = F, sep = "\t", col.names = F)
  
}

for(i in c(1:5, 7, 9:11)){
  iSpatial_PCC(paste0("Dataset", i))
}

#iSpatial_PCC("Dataset16") # dataset8 does not contain RAW data
#iSpatial_PCC("Dataset17") # dataset8 does not contain RAW data


### plotting
PCC = lapply(c(1:5, 7, 9:11), function(i){
  print(i)
  PCC = read.table(paste0("../data/benchmark/Dataset", i, "/iSpatial_cor.txt"))
  PCC = cbind(PCC$V2, i, "iSpatial")
  return(PCC)
})
PCC = do.call(rbind, PCC)
PCC = as.data.frame(PCC)
colnames(PCC) = c("PCC", "Dataset", "Method")


# From paper
for (i in c(1:5, 7, 9:15)) {
  print(i)
  for (j in c("Tangram", "gimVI", "SpaGE", "Seurat", "novoSpaRc", "SpaOTsc", "LIGER", "stPlus")){
    print(j)
    files = paste0("/nfs4/chaozhang/proj/methods/integrated/SpatialBenchmarking-main/scRNANorm_SpatialNorm/Data", i, "/", j,"_Metrics.txt")
    if(file.exists(files)){
      PCC1 = read.table(files, header = T, sep="\t")
      PCC1 = cbind(PCC1$PCC, i, j)
      PCC1 = as.data.frame(PCC1)
      colnames(PCC1) = c("PCC", "Dataset", "Method")
      
      PCC = rbind(PCC, PCC1)
    }
  }
}

PCC = na.omit(PCC)
PCC$PCC = as.numeric(PCC$PCC)
PCC$Dataset = factor(PCC$Dataset, levels = c(1:5, 7, 9:11))
PCC$Method = factor(PCC$Method, levels = c("iSpatial", "Tangram", "gimVI", "SpaGE", "Seurat", "novoSpaRc", "SpaOTsc", "LIGER", "stPlus"))
ggplot(PCC[PCC$Dataset %in% 1:11,], aes(Dataset, PCC, fill=Method))+
  geom_boxplot(outlier.size = -1,  position = position_dodge(0.8)) +
  ylim(-0.05,0.85) + xlab("Datasets") + ylab("Pearson's correlation") +
  stat_summary(fun=mean, geom="point", shape=18, size=3, position = position_dodge(0.8) ) +
  cowplot::theme_cowplot() + 
  scale_fill_manual(values = as.character(pals::brewer.pastel1(20)))
ggsave("../figure/benchmark_9datasets.pdf", width = 20, height = 5)
  
  
ggplot(PCC[PCC$Dataset %in% 1:11,], aes(Method, PCC, fill=Method))+
  geom_boxplot(outlier.size = -1,  position = position_dodge(0.8)) +
  ylim(-0.05,0.85) + xlab("") + ylab("Pearson's correlation") +
  stat_summary(fun=mean, geom="point", shape=18, size=3, position = position_dodge(0.8) ) +
  cowplot::theme_cowplot() + 
  scale_fill_manual(values = as.character(pals::brewer.pastel1(20)))
ggsave("../figure/benchmark_merge9datasets.pdf", width = 10, height = 5)

wilcox.test(PCC[PCC$Method == "iSpatial", 1], PCC[PCC$Method == "Tangram", 1],) # p-value < 2.2e-16

library(dplyr)
PCC_rank = PCC %>% group_by(Dataset, Method) %>%
  dplyr::summarise(mean(PCC))

for(i in unique(PCC_rank$Dataset)){
  iSpatial_rank = PCC_rank[PCC_rank$Dataset == i, 3]
  print(order(iSpatial_rank, decreasing = T)[1])
}



