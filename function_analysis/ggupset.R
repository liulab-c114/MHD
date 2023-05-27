rm(list = ls())
setwd('~/project/singleCell/MHD_version2/MHD_hypo/')

library(ggupset)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(tidyverse)
library(patchwork)

ob_deg <- read.csv('agrp/23ob_deg.csv', row.names = 1) %>% filter(abs(avg_log2FC) > log2(1.5))
db_deg <- read.csv('agrp/23db_deg.csv', row.names = 1) %>% filter(abs(avg_log2FC) > log2(1.5))

ob_deg$gene <- rownames(ob_deg)
ob_deg %>% mutate(type = ifelse(avg_log2FC > 0, "Obesity upregulate", "Obesity downregulate")) %>% 
  dplyr::select(c(gene, type)) -> ob.df

db_deg$gene <- rownames(db_deg)
db_deg %>% mutate(type = ifelse(avg_log2FC > 0, "Diabetes upregulate", "Diabetes downregulate")) %>% 
  dplyr::select(c(gene, type)) -> db.df

merge.df <- rbind(ob.df, db.df) %>% 
  dplyr::group_by(gene) %>%
  distinct(type, .keep_all = T) %>%
  summarise(cls = list(type))

# 绘制DEG的Overlap
bar_deg <- ggplot(merge.df, aes(x = cls)) + 
  geom_bar(
    fill = c(rep("black", 2), rep("#f28073", 2), rep("black", 2))
    # fill = c("#f28073", "black", "#f28073", rep("black", 3))
  ) + 
  geom_text(stat='count', aes(label=..count..)
            , vjust=c(rep(1.5, 6))
            , colour = c(rep("white", 6))
            , size = 3.2)+
  theme_bw() + 
  labs(title = "Number of DEG", x = "", y = "") +
  scale_x_upset() + 
  theme(panel.grid = element_blank())
ggsave(filename = "agrp/agrp_degBar.pdf", plot = bar_deg, width = 6, height = 6)
