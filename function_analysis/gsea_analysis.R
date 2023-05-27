setwd("/data/liangxian/project/singleCell/GSEA")
Sys.time()
rm(list = ls())

# 加载数据集

library(Seurat)
# library(enrichplot)
library(patchwork)
# library(simplifyEnrichment)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

cls_of_interest <- c(23, 40, 9, 26, 33)

# 读取DEG
ob_deg <- list()
for(cls in 1:length(cls_of_interest)){
  ob_deg[[cls]] <- read.csv(paste0("data/", cls_of_interest[cls], "ob_deg.csv"))
  ob_deg[[cls]] <- arrange(ob_deg[[cls]], desc(avg_log2FC))
  ob_deg[[cls]]$cls <- cls_of_interest[cls]
}
# names(ob_deg) <- cls_of_interest

# 读取 gmt文件
cp.gmt <- read.gmt("data/c2.cp.v7.5.1.symbols.gmt")

# GSEA富集分析 cp enrich
gse.CP.list <- list()
cp.ob.res <- list()

for(cls in 1:length(cls_of_interest)){
  data_all_sort <- ob_deg[[cls]] %>%
    arrange(desc(avg_log2FC))
  head(data_all_sort)

  geneList = data_all_sort$avg_log2FC #把foldchange按照从大到小提取出来
  names(geneList) <- data_all_sort$X #给上面提取的foldchange加上对应上ENTREZID
  gse.CP <- GSEA(geneList, TERM2GENE = cp.gmt, pvalueCutoff = 1)
  # saveRDS(gse.GO, file = paste0("outputs/gsea/", cls_of_interest[cls], "gse_go.rds"))
  gse.CP.list[[cls]] <- gse.CP
  cp.ob.res[[cls]] <- gse.CP@result
}

for(cls in 1:length(cls_of_interest)){
  cp.ob.res[[cls]]$cls <- cls_of_interest[cls]
  cp.ob.res[[cls]]$group <- "obesity"
}
cp.ob.df <- do.call(rbind, cp.ob.res)
cp.ob.df$Description <- tolower(cp.ob.df$Description)
write.csv(cp.ob.df, file = "outputs/gsea_ob_cp.csv")

# go enrichment analysis

# GSEA富集分析
gse.GO.list <- list()
go.ob.res <- list()

for(cls in 1:length(cls_of_interest)){
  colnames(ob_deg[[cls]])[1] <- "SYMBOL"
  gene <- ob_deg[[cls]]$SYMBOL
  #开始ID转换，会有丢失
  gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  #去重
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

  data_all <- ob_deg[[cls]] %>%
    inner_join(gene,by="SYMBOL")
  dim(data_all)
  head(data_all)

  data_all_sort <- data_all %>%
    arrange(desc(avg_log2FC))
  head(data_all_sort)

  geneList = data_all_sort$avg_log2FC #把foldchange按照从大到小提取出来
  names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
  gse.GO <- gseGO(
    geneList, #geneList
    ont = "BP",  # 可选"BP"、"MF"和"CC"或"ALL"
    OrgDb = org.Hs.eg.db, #人 注释基因
    keyType = "ENTREZID",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",#p值校正方法
  )
  # saveRDS(gse.GO, file = paste0("week2/gsea/", cls_of_interest[cls], "ob_gse_go.rds"))
  gse.GO.list[[cls]] <- gse.GO
  go.ob.res[[cls]] <- gse.GO@result
}


for(cls in 1:length(cls_of_interest)){
  go.ob.res[[cls]]$cls <- cls_of_interest[cls]
  go.ob.res[[cls]]$group <- "obesity"
}
go.ob.df <- do.call(rbind, go.ob.res)

go.df <- rbind(go.ob.df, go.ob.df)
write.csv(go.ob.df, file = "outputs/gsea_ob_go.csv")

# ==================== Diabetes ========================== #

db_deg <- list()
for(cls in 1:length(cls_of_interest)){
  db_deg[[cls]] <- read.csv(paste0("data/", cls_of_interest[cls], "db_deg.csv"))
  db_deg[[cls]] <- arrange(db_deg[[cls]], desc(avg_log2FC))
  db_deg[[cls]]$cls <- cls_of_interest[cls]
}

# GSEA富集分析 cp enrich
gse.CP.list <- list()
cp.db.res <- list()

for(cls in 1:length(cls_of_interest)){
  data_all_sort <- db_deg[[cls]] %>%
    arrange(desc(avg_log2FC))
  head(data_all_sort)

  geneList = data_all_sort$avg_log2FC #把foldchange按照从大到小提取出来
  names(geneList) <- data_all_sort$X #给上面提取的foldchange加上对应上ENTREZID
  gse.CP <- GSEA(geneList, TERM2GENE = cp.gmt, pvalueCutoff = 1)
  # saveRDS(gse.GO, file = paste0("outputs/gsea/", cls_of_interest[cls], "gse_go.rds"))
  gse.CP.list[[cls]] <- gse.CP
  cp.db.res[[cls]] <- gse.CP@result
}

for(cls in 1:length(cls_of_interest)){
  cp.db.res[[cls]]$cls <- cls_of_interest[cls]
  cp.db.res[[cls]]$group <- "diabetes"
}
cp.db.df <- do.call(rbind, cp.db.res)
cp.db.df$Description <- tolower(cp.db.df$Description)
write.csv(cp.db.df, file = "outputs/gsea_db_cp.csv")

# go enrichment analysis

# GSEA富集分析
gse.GO.list <- list()
go.db.res <- list()

for(cls in 1:length(cls_of_interest)){
  colnames(db_deg[[cls]])[1] <- "SYMBOL"
  gene <- db_deg[[cls]]$SYMBOL
  #开始ID转换，会有丢失
  gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  #去重
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

  data_all <- db_deg[[cls]] %>%
    inner_join(gene,by="SYMBOL")
  dim(data_all)
  head(data_all)

  data_all_sort <- data_all %>%
    arrange(desc(avg_log2FC))
  head(data_all_sort)

  geneList = data_all_sort$avg_log2FC #把foldchange按照从大到小提取出来
  names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
  gse.GO <- gseGO(
    geneList, #geneList
    ont = "BP",  # 可选"BP"、"MF"和"CC"或"ALL"
    OrgDb = org.Hs.eg.db, #人 注释基因
    keyType = "ENTREZID",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",#p值校正方法
  )
  # saveRDS(gse.GO, file = paste0("week2/gsea/", cls_of_interest[cls], "db_gse_go.rds"))
  gse.GO.list[[cls]] <- gse.GO
  go.db.res[[cls]] <- gse.GO@result
}


for(cls in 1:length(cls_of_interest)){
  go.db.res[[cls]]$cls <- cls_of_interest[cls]
  go.db.res[[cls]]$group <- "diabetes"
}
go.db.df <- do.call(rbind, go.db.res)

go.df <- rbind(go.db.df, go.db.df)
write.csv(go.db.df, file = "outputs/gsea_db_go.csv")

Sys.time()
