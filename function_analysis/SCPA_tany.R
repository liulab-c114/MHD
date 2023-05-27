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

db_act <- read.csv('outputs/tany_SCPA.csv', row.names = 1)

plot_rank(scpa_out = db_act, 
          pathway = c("dna_damage"), 
          base_point_size = 2, 
          highlight_point_color = '#2980b9', 
          highlight_point_size = 3, 
          label_pathway = F) + theme_classic()
ggsave('tany_dna_damage_rank.pdf', width = 4, height = 3)

db_act <- db_act %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

ggplot(db_act, aes(FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_jitter(color = db_act$color, size = 0.5
              # , cex = 0.1, shape = 21, , stroke = 0.3
              ) +
  # geom_point(data = aa_path, shape = 21, cex = 2.8, fill = "orangered2", color = "black", stroke = 0.3) +
  xlim(-80, 10) +
  ylim(0, 7) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)
