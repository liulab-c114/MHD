rm(list = ls())

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(patchwork)

plot_deg_heat <- function(subtype_name){
  # Obesity
  ob_list <- list()
  for(cls in 1:length(subtype_name)){
    ob_list[[cls]] <- read.csv(paste0("data/nn_DEG/", subtype_name[cls], "_ob_normal_all_genes.csv"))
    colnames(ob_list[[cls]])[1] <- "gene"
    ob_list[[cls]]$Cluster <- subtype_name[cls]
    ob_list[[cls]]$Group <- "OB"
  }
  
  ob_df <- do.call("rbind", ob_list) %>% 
    filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.05)
  
  # Diabetes
  
  db_list <- list()
  for(cls in 1:length(subtype_name)){
    db_list[[cls]] <- read.csv(paste0("data/nn_DEG/", subtype_name[cls], "_db_normal_all_genes.csv"))
    colnames(db_list[[cls]])[1] <- "gene"
    db_list[[cls]]$Cluster <- subtype_name[cls]
    db_list[[cls]]$Group <- "DB"
  }
  
  db_df <- do.call("rbind", db_list) %>% 
    filter(abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.05)
  
  deg_mt <- rbind(ob_df, db_df) %>%
    mutate(Cls_group = paste(Cluster, Group, sep = "_")) %>%
    select(gene, Cls_group, avg_log2FC) %>% 
    spread(key = Cls_group, value = avg_log2FC) %>%
    column_to_rownames('gene') %>% as.matrix()
  deg_mt[is.na(deg_mt)] <- 0
  
  col_fun = colorRamp2(c(-2, 0, 2), c("#0097e6", "#ecf0f1", "#c23616"))
  
  deg_heatmap <- Heatmap(deg_mt, col = col_fun
                         , rect_gp = gpar(col = "#7f8c8d", lwd = 0.5)
                         , name = 'logFC')
  deg_heatmap
}

ast <- c('Ast-1', 'Ast-2')
fibro <- c('Fibro-1', 'Fibro-2')
# c('tany_a', 'tany_b', 'micro-1', 'micro-2')
endo <- c('Endothelial')
ependy <- c('Ependy-1', 'Ependy-2')
oligo <- c('Oligo-1', 'Oligo-2')
opc <- c('OPC')

ast_heatmap <- plot_deg_heat(ast)
fibro_heatmap <- plot_deg_heat(fibro)
endo_heatmap <- plot_deg_heat(endo)
ependy_heatmap <- plot_deg_heat(ependy)
oligo_heatmap <- plot_deg_heat(oligo)
opc_heatmap <- plot_deg_heat(opc)

pdf('outputs/nn_deg_heatmap/ast_heatmap.pdf', height = 4, width = 4)
draw(ast_heatmap)
dev.off()

pdf('outputs/nn_deg_heatmap/fibro_heatmap.pdf', height = 20, width = 4)
draw(fibro_heatmap)
dev.off()

pdf('outputs/nn_deg_heatmap/endo_heatmap.pdf', height = 5, width = 3)
draw(endo_heatmap)
dev.off()

pdf('outputs/nn_deg_heatmap/ependy_heatmap.pdf', height = 8, width = 4)
draw(ependy_heatmap)
dev.off()

pdf('outputs/nn_deg_heatmap/oligo_heatmap.pdf', height = 6, width = 4)
draw(oligo_heatmap)
dev.off()

pdf('outputs/nn_deg_heatmap/opc_heatmap.pdf', height = 4, width = 3)
draw(opc_heatmap)
dev.off()
