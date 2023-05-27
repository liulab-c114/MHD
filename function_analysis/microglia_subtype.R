# perform microglia analysis

rm(list = ls())

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(patchwork)

plotClsRatio <- function(cls, info, manual_col = F, coltype = "fill", colornum = 8, colortype = "Paired", man_cols = vector()){
  require(RColorBrewer)
  colourCount = length(unique(info))
  getPalette = colorRampPalette(brewer.pal(colornum, colortype))
  dataFrame <- as.data.frame(table(cls, info))
  colnames(dataFrame) <- c("clusters", "Group", "Number")
  
  if(manual_col){
    fill_color <- man_cols
  }else{
    fill_color <- getPalette(colourCount)
  }
  
  ggplot(dataFrame, aes(x = clusters, y = Number, fill = Group)) + 
    geom_col(position = coltype, color = "black") + 
    scale_fill_manual(values = fill_color) + 
    coord_flip() + 
    theme(panel.background = element_blank()
          , axis.text = element_text(face = "bold")
          , axis.text.x = element_text(angle = 45)
    )
}


micro.seu <- readRDS('outputs/microglia/micro_seu.rds')

DimPlot(micro.seu, group.by = 'sub_type')

clsMarkers <- FindMarkers(micro.seu, logfc.threshold = 0
                          , ident.1 = 'micro-1', ident.2 = 'micro-2', group.by = 'sub_type'
                          , test.use = 'MAST')

write.csv(clsMarkers, 'outputs/micro_clsmarkers_mast.csv')

clsMarkers_1 <- FindMarkers(micro.seu, logfc.threshold = 0
                          , ident.1 = 'micro-1', ident.2 = 'micro-2', group.by = 'sub_type')

write.csv(clsMarkers, 'outputs/micro_clsmarkers_wilcoxon.csv')

micro.seu$nuclei <- factor(micro.seu$nuclei, levels = c('PVN', 'DMH', "LAH", 'VMH', 'IN')
                          , labels = c('PVN', 'DMH', 'LHA', 'VMH', 'INF'))
table(micro.seu$nuclei)
plotClsRatio(micro.seu$nuclei, micro.seu$sub_type
             , manual_col = T
             , man_cols = c('#d35400', '#ff793f')) + NoLegend()
ggsave('outputs/microglia/subtype_ratio.pdf', width = 6, height = 2)

rm(list=ls())
library(tidyverse)
library(patchwork)

clsMarker <- read.csv('outputs/micro_clsmarkers_mast.csv', row.names = 1)
clsMarker_1 <- read.csv('outputs/micro_clsmarkers_wilcoxon.csv', row.names = 1)

volcanoPlot <- function(data, logFC_thre = 1, pVal_thre = 0.05, plot.title = "", save.dir = ""){
  ## 用于绘制火山图的函数，data为需要传入的数据框；
  # logFC_thre为logfoldchange的阈值；
  # pVal_thre为p值(校正后)的阈值；
  # title为所绘制的火山图的标题。
  
  #获得表达上调和下调的基因
  data <- transform(data, p_val_adj = -1*log10(data$p_val_adj))
  data <- data[sort(data$avg_log2FC,index.return = T)$ix,]
  data$label <- ''
  data$label[data$avg_log2FC<= -logFC_thre & data$p_val_adj > -1*log10(pVal_thre)] <- 'down'
  
  data$label[data$avg_log2FC >= logFC_thre & data$p_val_adj > -1*log10(pVal_thre)] <- 'up'
  
  data$label[data$label == ''] <- 'no'
  
  #显示表达上调和下调的代表性几个基因
  
  # data$name <- ''
  # 
  # up <- row.names(tail(data[data$label=='up', ], 50))
  # down <- row.names(head(data[data$label=='down', ], 50))
  # labeled_gene <- c(up,down)
  # for(j in 1:nrow(data))
  # {
  #   for(i in labeled_gene)
  #   {
  #     if(row.names(data[j,]) == i)
  #     {
  #       data[j,'name'] <- row.names(data[j,])
  #     }
  #   }
  # }
  
  #绘制火山图
  
  set_theme <-  theme(panel.background = element_blank()
                      ,plot.title = element_text(hjust = 0.5, size = 18)
                      ,axis.title.x = element_text(vjust=-0.45, size = 16)
                      ,axis.title.y = element_text(vjust=1.2, size = 16)
                      ,axis.text = element_text(size = 8)
                      ,axis.ticks = element_line(colour="black")
                      ,axis.line = element_line())
  
  volcanoP <- ggplot(data = data, aes(avg_log2FC, p_val_adj))+
    geom_point(size = 1.5, aes(color = label))+
    scale_color_manual(values = c('#2980b9', '#7f8c8d', '#c0392b'))+
    geom_vline(xintercept = c(-logFC_thre,logFC_thre), lty = 4, col = '#2c3e50', lwd = 0.5)+
    geom_hline(yintercept = -log10(pVal_thre), lty = 4, col = '#2c3e50', lwd = 0.5)+
    labs(x = 'log2(FoldChange)', y = '-log10(P_adj)', title = plot.title)+
    # geom_text_repel(aes(label = name), max.overlaps = 10)+
    set_theme
  savefile <- paste(save.dir, plot.title, pVal_thre,'.png',sep = '')
  ggsave(volcanoP, filename = savefile, dpi = 300, height = 15, width = 15)
  return(volcanoP)
}

p_mast <- volcanoPlot(clsMarker, logFC_thre = 0.3, pVal_thre = 0.05, plot.title = 'microglia_subtype_mast', save.dir = 'outputs/')
p_wilcox <-  volcanoPlot(clsMarker, logFC_thre = 0.3, pVal_thre = 0.05, plot.title = 'microglia_subtype_wilcox', save.dir = 'outputs/')

ggsave(filename = 'outputs/volcano_p_mast.pdf', width = 6, height = 4)

