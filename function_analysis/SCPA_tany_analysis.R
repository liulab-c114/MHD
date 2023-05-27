library(SCPA)
library(msigdbr)
library(Seurat)
library(tidyverse)
library(patchwork)

rm(list = ls())

setwd('/data/liangxian/project/MHD_version2/')

tany.seu <- readRDS('outputs/tany_seu.rds')

pathways <- 'data/c5.go.bp.v7.5.1.symbols.gmt'

tany.con.seu <- seurat_extract(tany.seu, meta1 = 'diabetes', value_meta1 = 'Control')
tany.db.seu <- seurat_extract(tany.seu, meta1 = 'diabetes', value_meta1 = 'Diabetes')

db_act <- compare_pathways(samples = list(tany.con.seu, tany.db.seu), 
                           pathways = pathways)

write.csv(db_act, 'outputs/tany_SCPA.csv')
