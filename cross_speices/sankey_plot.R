rm(list = ls())
setwd('~/project/singleCell/MHD_version2/MHD_hypo/')
library(ggsankey)
library(tidyverse)
library(patchwork)

pomc.meta <- read.csv('pomc/pomc_metadata.csv', row.names = 1)
monkey.meta <- filter(pomc.meta, species == '0')


df <- monkey.meta %>%
  make_long(subcls_name, predicted_C185_named)

pl <- ggplot(df, aes(x = x,                        
                     next_x = next_x,                                     
                     node = node,
                     next_node = next_node,        
                     fill = factor(node)))

pl <- pl +geom_sankey(flow.alpha = 0.5,          #This Creates the transparency of your node 
                      node.color = "black",     # This is your node color        
                      show.legend = TRUE)        # This determines if you want your legend to show

pl <- pl + theme_bw()
pl <- pl + theme(legend.position = 'none')
pl <- pl + theme(axis.title = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),
                 panel.grid = element_blank())

# pl <- pl + scale_fill_viridis_d(option = "inferno")
# pl <- pl + scale_fill_manual(values = c('POMCh_LEPRh' = "#c0392b",
#                                         'POMCh_LEPRl' = "#f39c12",
#                                         'POMCl_LEPRl' = "#2980b9",
#                                         'C185-48: Anxa2.Pomc.GLU-5' = "#d35400",
#                                         'C185-49: Glipr1.Pomc.GLU-5' = "#8e44ad",
#                                         'C185-50: Ttr.Pomc.GLU-5' = "#27ae60"))

pl
ggsave('pomc/sankey_plot.pdf', pl, height = 4, width = 4)
