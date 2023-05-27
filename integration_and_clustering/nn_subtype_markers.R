setwd('/data/liangxian/project/MHD_version2/')

library(Seurat)
library(tidyverse)
non.seu <- readRDS('data/non_neuron_named.rds')
Idents(non.seu) <- 'sub_type'
clsmarkers <- FindAllMarkers(non.seu, logfc.threshold = log2(1.5), only.pos = T, test.use = 'MAST')
write.csv(clsmarkers, 'data/non_neuron_subtype_markers.csv')
