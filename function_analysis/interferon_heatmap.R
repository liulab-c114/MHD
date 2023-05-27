rm(list = ls())

setwd('~/project/singleCell/MHD_version2/')

library(ComplexHeatmap)
library(reshape2)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(patchwork)

inf.genes <- read.csv('CCI/interferon_singnaling.csv'
                      , header = F, check.names = F) %>% pull(V1)

# Neuron_interferon_signaling

cls_of_interest <- c(23, 40, 9, 26, 33)
neuron_list <- list()

for(cls in 1:length(cls_of_interest)){
  neuron_list[[cls]] <- read.csv(paste0("CCI/deg_list/", cls_of_interest[cls], "db_deg.csv"))
  colnames(neuron_list[[cls]])[1] <- "gene"
  neuron_list[[cls]]$Cluster <- cls_of_interest[cls]
  neuron_list[[cls]]$Group <- "db"
}

neuron_df <- do.call("rbind", neuron_list) %>% 
  filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.05 & gene %in% inf.genes) %>% 
  select(gene, Cluster, avg_log2FC)

neuron_df$Cluster <- factor(neuron_df$Cluster, levels = c(33, 26, 9, 40, 23))

neuron_df <- spread(neuron_df, key = Cluster, value = avg_log2FC) %>% 
  column_to_rownames('gene') %>% as.matrix()
neuron_df[is.na(neuron_df)] <- 0

col_fun = colorRamp2(c(0, 2), c("#ecf0f1", "#c0392b"))

pdf('CCI/interferon_neuron.pdf', height = 5, width = 3)
Heatmap(neuron_df, cluster_rows = T
        , cluster_columns = F
        , rect_gp = gpar(col = "#7f8c8d", lwd = 0.5)
        , name = 'logFC'
        , col = col_fun)
dev.off()

# Tanycyte_interferon_signaling

tany_cls <- c('tany_a', 'tany_b')
tany_list <- list()

for(cls in 1:length(tany_cls)){
  tany_list[[cls]] <- read.csv(paste0("CCI/deg_list/", tany_cls[cls], "_db_normal_all_genes.csv"))
  colnames(tany_list[[cls]])[1] <- "gene"
  tany_list[[cls]]$Cluster <- tany_cls[cls]
  tany_list[[cls]]$Group <- "db"
}

tany_df <- do.call("rbind", tany_list) %>% 
  filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.05 & gene %in% inf.genes) %>% 
  select(gene, Cluster, avg_log2FC)

tany_df$Cluster <- factor(tany_df$Cluster, levels = c('tany_a', 'tany_b'))

tany_df <- spread(tany_df, key = Cluster, value = avg_log2FC) %>% 
  column_to_rownames('gene') %>% as.matrix()
tany_df[is.na(tany_df)] <- 0

col_fun = colorRamp2(c(0, 2), c("#ecf0f1", "#2980b9"))
pdf('CCI/interferon_tanycyte.pdf', height = 5, width = 2.5)
Heatmap(tany_df, cluster_rows = T
        , cluster_columns = F
        , name = 'logFC'
        , rect_gp = gpar(col = "#7f8c8d", lwd = 0.5)
        , col = col_fun)
dev.off()

# Microglia_interferon_signaling

micro_cls <- c('micro-1', 'micro-2')
micro_list <- list()

for(cls in 1:length(micro_cls)){
  micro_list[[cls]] <- read.csv(paste0("CCI/deg_list/", micro_cls[cls], "_db_normal_all_genes.csv"))
  colnames(micro_list[[cls]])[1] <- "gene"
  micro_list[[cls]]$Cluster <- micro_cls[cls]
  micro_list[[cls]]$Group <- "db"
}

micro_df <- do.call("rbind", micro_list) %>% 
  filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.05 & gene %in% inf.genes) %>% 
  select(gene, Cluster, avg_log2FC)

micro_df$Cluster <- factor(micro_df$Cluster, levels = c('micro-1', 'micro-2'))

micro_df <- spread(micro_df, key = Cluster, value = avg_log2FC) %>% 
  column_to_rownames('gene') %>% as.matrix()
micro_df[is.na(micro_df)] <- 0

col_fun = colorRamp2(c(0, 2), c("#ecf0f1", "#16a085"))
pdf('CCI/interferon_microglia.pdf', height = 5, width = 2.5)
Heatmap(micro_df, cluster_rows = T
        , cluster_columns = F
        , name = 'logFC'
        , rect_gp = gpar(col = "#7f8c8d", lwd = 0.5)
        , col = col_fun)
dev.off()
