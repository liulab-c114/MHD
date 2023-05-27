library(SCPA)
library(msigdbr)
library(Seurat)
library(tidyverse)
library(patchwork)

rm(list = ls())

setwd('/data/liangxian/project/MHD_version2/')

pvn.seu <- readRDS('data/pvn_seu.rds')

pathways <- 'data/metabolic_pathways.csv'

pvn.con.seu <- seurat_extract(pvn.seu, meta1 = 'diabetes', value_meta1 = 'normal')
pvn.db.seu <- seurat_extract(pvn.seu, meta1 = 'diabetes', value_meta1 = 'diabetes')

db_act <- compare_pathways(samples = list(pvn.con.seu, pvn.db.seu), 
                             pathways = pathways)

db_act <- db_act %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

aa_path <- db_act %>% 
  filter(grepl(pattern = "reactome_arachi", ignore.case = T, x = Pathway))

ggplot(db_act, aes(FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = db_act$color, stroke = 0.3) +
  # geom_point(data = aa_path, shape = 21, cex = 2.8, fill = "orangered2", color = "black", stroke = 0.3) +
  xlim(-10, 45) +
  ylim(0, 5.1) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)

ggsave('outputs/pvn_SCPA_metab.pdf', width = 3, height = 4)

write.csv(db_act, 'outputs/pvn_SCPA.csv')
