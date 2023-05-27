rm(list = ls())
setwd('~/project/singleCell/MHD_version2/')
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(patchwork)

inf_gsea <- read.csv('data/gsea_cls_INF.csv', row.names = 1)

immune_pathway <- c('defense response'
                    # , 'response to virus'
                    # , 'regulation of type I interferon production'
                    , 'type I interferon production'
                    , 'response to cytokine'
                    , 'response to type I interferon'
                    , 'innate immune response')

synapse_pathway <- c('chemical synaptic transmission', 'synapse assembly'
                     # , 'central nervous system development'
                     # , 'regulation of neuron projection development'
                     # , 'dendrite development'
                     , 'regulation of synaptic plasticity'
                     , 'synapse organization', 'synaptic signaling')

immune_df <- filter(inf_gsea, Description %in% immune_pathway) %>% 
  mutate(cluster = 'immune'
         , group_cls = paste(group, cls, sep = '_')) %>% 
  select(Description, NES, group_cls)

synapse_df <- filter(inf_gsea, Description %in% synapse_pathway) %>% 
  mutate(cluster = 'synapse'
  , group_cls = paste(group, cls, sep = '_')) %>% 
  select(Description, NES, group_cls)

selected_df <- rbind(immune_df, synapse_df)

selected_df <- spread(selected_df, key = group_cls, value = NES) %>% 
  column_to_rownames('Description') %>% 
  t()
col_fun = colorRamp2(c(-3, 0, 3), c("#0097e6", "#ecf0f1", "#c23616"))
col_fun(seq(-3, 3))

Heatmap(selected_df, cluster_rows = F, name = 'NES', col = col_fun
        , na_col = '#ecf0f1'
        , rect_gp = gpar(col = "#2c3e50", lwd = 0.5))

pdf('outputs/inf_pathway_heatmap.pdf', width = 5, height = 3)
Heatmap(selected_df, cluster_rows = F, name = 'NES', col = col_fun
        , show_row_names = F, show_column_names = F, na_col = '#ecf0f1'
          , rect_gp = gpar(col = "#2c3e50", lwd = 0.5))
dev.off()



# install.packages("devtools")
devtools::install_version("crossmatch", version = "1.3.1", repos = "http://cran.us.r-project.org")
devtools::install_version("multicross", version = "2.1.0", repos = "http://cran.us.r-project.org")
devtools::install_github("jackbibby1/SCPA")

