options(warn=-1)
suppressMessages(library(Seurat))
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))

option_list <- list(  # parameters list
  make_option(c("-f", "--file"), type = "character", help = "input rds file path"),
  # make_option(c("-o", "--output"), type = "character", help = "output file path"),
  make_option(c("-b", "--batch"), type = "character", help = "key name of batch"),
  make_option(c("-r", "--resolution"), type = "double", default = 0.5, help = "resolution of louvain clustering")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

seu.rds <- opt$file
batch <- opt$batch
resolution <- opt$resolution

CCA_intergrate <- function(object, split.by, npcs = 40, dims = 1:30, resolution = 0.5, ...){
  
  ifnb.list <- SplitObject(object, split.by = split.by)
  
  ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst")
  })
  features <- SelectIntegrationFeatures(object.list = ifnb.list)
  immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list
                                           , anchor.features = features, ...)
  immune.combined <- IntegrateData(anchorset = immune.anchors)
  DefaultAssay(immune.combined) <- "integrated"
  # Run the standard workflow for visualization and clustering
  immune.combined <- ScaleData(immune.combined, verbose = FALSE)
  immune.combined <- RunPCA(immune.combined, npcs = npcs, verbose = FALSE)
  immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = dims)
  immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = dims)
  immune.combined <- FindClusters(immune.combined, resolution = resolution)
  DefaultAssay(immune.combined) <- "RNA"
  return(immune.combined)
}

seu.obj <- readRDS(seu.rds)
seu.obj <- CCA_intergrate(seu.obj, batch, resolution = resolution)
saveRDS(seu.obj, seu.rds)
