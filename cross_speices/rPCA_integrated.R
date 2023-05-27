# Thu Nov  4 10:34:28 2021 ------------------------------
# coding = utf-8
# 进行整合聚类的标准流程

rPCA_integrated <- function(object, split.by, reference, npcs = 40, dims = 1:30, resolution = 0.5){
  
  ifnb.list <- SplitObject(object, split.by = split.by)
  
  ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  features <- SelectIntegrationFeatures(object.list = ifnb.list)
  ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  anchors <- FindIntegrationAnchors(object.list = ifnb.list, reference = reference, reduction = "rpca",
                                    dims = dims)
  immune.combined <- IntegrateData(anchorset = anchors, dims = dims)
  
  immune.combined <- ScaleData(immune.combined, verbose = FALSE)
  immune.combined <- RunPCA(immune.combined, npcs = npcs, verbose = FALSE)
  immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = dims)
  immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = dims)
  immune.combined <- FindClusters(immune.combined, resolution = resolution)
  DefaultAssay(immune.combined) <- "RNA"
  return(immune.combined)
}
