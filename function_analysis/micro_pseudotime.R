library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)

load('outputs/micro_cds.Rdata')

pData(cds)$Cluster=pData(cds)$diabetes

Sys.time()
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Cluster")
Sys.time()

sig_genes <- subset(diff_test_res, qval < 1e-10)
sig_genes <- sig_genes[order(sig_genes$pval),]
# head(sig_genes[,c("gene_short_name", "pval", "qval")] ) 
cg <- as.character(head(sig_genes$gene_short_name, 25)) 

# 第一步: 挑选合适的基因. 有多个方法，例如提供已知的基因集，这里选取统计学显著的差异基因列表
ordering_genes <- row.names(subset(diff_test_res, qval < 1e-5))
cds <- setOrderingFilter(cds, ordering_genes)

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)

save(cds, file = 'outputs/micro_processed_cds.Rdata')

pdf('outputs/micro_pseudotime.pdf')
plot_cell_trajectory(cds, color_by = "Pseudotime")
dev.off()

pdf('outputs/micro_cluster.pdf')
plot_cell_trajectory(cds, color_by = "Cluster")
dev.off()

pdf('outputs/micro_state.pdf')
plot_cell_trajectory(cds, color_by = "State")
dev.off()

