library(Seurat)
library(CellChat)
library(tidyverse)
library(patchwork)

seu.obj <- readRDS('outputs/inf_neu_tany_micro_seuobj.rds')

group <- c("Control","Obesity","Diabetes")
data.obj = seu.obj@assays$RNA@data # normalized data matrix
objmeta = seu.obj@meta.data # a dataframe with rownames containing cell mata data

object.list <- list()

for (i in 1:3) {
  g <- group[i]
  setwd(paste0("/data/liangxian/project/MHD_version2/outputs/cellchat/", group[i],"/"))
  cell.use = rownames(objmeta)[objmeta$diabetes == g] # extract the cell names from disease data
  # Prepare input data for CelChat analysis
  data.input = data.obj[, cell.use]
  meta = objmeta[cell.use, ]
  unique(meta$subtype) # check the cell labels
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "subtype4cc")
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "subtype4cc") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  showDatabaseCategory(CellChatDB)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multiprocess", workers = 8) # do parallel
  #> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
  #> explicitly specify either 'multisession' or 'multicore'. In the current R
  #> session, 'multiprocess' equals 'multisession'.
  #> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
  #> processing ('multicore') is not supported when running R from RStudio
  #> because it is considered unstable. For more details, how to control forked
  #> processing or not, and how to silence this warning in future R sessions, see ?
  #> parallelly::supportsMulticore
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  cellchat <- projectData(cellchat, PPI.human)
  
  cellchat <- computeCommunProb(cellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  
  saveRDS(cellchat, paste0('cellchat_', g, '.rds'))
  
  pdf("1.netVisual_circle.pdf",width = 10)
  par(mfrow = c(1,2), xpd=TRUE)
  p1 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  p2 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  print(p2)
  dev.off()
  
  
  mat <- cellchat@net$weight
  pdf("2.celltype_weight.pdf",width = 15,height = 25)
  par(mfrow = c(5,3), xpd=TRUE)
  for (j in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[j, ] <- mat[j, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[j])
  }
  dev.off()
  
  #Part III: Visualization of cell-cell communication network
  print(cellchat@netP$pathways)
  pathways.show <- c("PSAP") 
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  vertex.receiver = seq(1,4) # a numeric vector. 
  
  netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  
  
  # Circle plot
  pdf("3.PSAP_pathway_network.pdf")
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  
  # Chord diagram
  pdf("4.PSAP_pathway_network_chord.pdf")
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  dev.off()
  
  
  # Heatmap
  pdf("5.PSAP_pathway_heat.pdf")
  par(mfrow=c(1,1))
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  dev.off()
  
  #> Do heatmap based on a single object
  
  # Chord diagram
  #group.cellType <- c(rep("MIC_3", 4), rep("MIC_4", 4), rep("n12:AGRP/MCTP2", 4)) # grouping cell clusters into fibroblast, DC and TC cells
  #names(group.cellType) <- levels(cellchat@idents)
  #netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
  #> Plot the aggregated cell-cell communication network at the signaling pathway level
  pdf("6.LR_contriubtion.pdf",height = 3)
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  dev.off()
  
  
  pairLR.PSAP <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
  LR.show <- pairLR.PSAP[1,] # show one ligand-receptor pair
  # Hierarchy plot
  vertex.receiver = seq(1,4) # a numeric vector
  pdf("7.PSAP_net.pdf")
  netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
  dev.off()
  #> [[1]]
  # Circle plot
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
  
  # Chord diagram
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
  
  #Automatically save the plots of the all inferred network for quick exploration
  pathways.show.all <- cellchat@netP$pathways
  # check the order of cell identity to set suitable vertex.receiver
  levels(cellchat@idents)
  vertex.receiver = seq(1,4)
  dir.create("8.pathway_net_analysis/")
  for (j in 1:length(pathways.show.all)) {
    # Visualize communication network associated with both signaling pathway and individual L-R pairs
    net <- netVisual(cellchat, signaling = pathways.show.all[j], vertex.receiver = vertex.receiver, layout = "hierarchy")
    # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
    gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[j])
    ggsave(filename=paste0("8.pathway_net_analysis/",pathways.show.all[j], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
  }
  
  #Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
  # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
  netVisual_bubble(cellchat, sources.use = c(3:7), targets.use = c(1:2), remove.isolate = FALSE)
  #> Comparing communications on a single object
  
  # show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
  pdf("9.allpathway_LR_bubble.pdf",width = 5,height = 5)
  pairLR.use <- extractEnrichedLR(cellchat, signaling = pathways.show.all)
  netVisual_bubble(cellchat, sources.use = c(3:7), targets.use = c(1:2), pairLR.use = pairLR.use, remove.isolate = TRUE)
  dev.off()
  #> Comparing communications on a single object
  
  # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
  # show all the interactions sending from Inflam.FIB
  pdf("10.source_mic_target_neu.pdf",width = 10)
  netVisual_chord_gene(cellchat, sources.use = c(3:7), targets.use = c(1:2), lab.cex = 0.5,legend.pos.y = 30)
  dev.off()
  
  # show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
  #netVisual_chord_gene(cellchat, sources.use = c(3:7), targets.use = c(1:2), slot.name = "netP", legend.pos.x = 10)
  
  #Plot the signaling gene expression distribution using violin/dot plot
  pdf("11.PSAP_gene_expression.pdf")
  #plotGeneExpression(cellchat, signaling = "PSAP")
  plotGeneExpression(cellchat, signaling = "PSAP", enriched.only = FALSE)
  dev.off()
  
  
  #> Registered S3 method overwritten by 'spatstat.geom':
  #>   method     from
  #>   print.boxx cli
  #> Scale for 'y' is already present. Adding another scale for 'y', which will
  #> replace the existing scale.
  #> Scale for 'y' is already present. Adding another scale for 'y', which will
  #> replace the existing scale.
  #> Scale for 'y' is already present. Adding another scale for 'y', which will
  #> replace the existing scale.
  
  #Part IV: Systems analysis of cell-cell communication network
  # Compute the network centrality scores
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
  pdf("12.PSAP_netAnalysis_signalingRole_network.pdf",height = 3,width = 5)
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  dev.off()
  dir.create("12.all_pathway_netAnalysis/")
  for (j in 1:length(pathways.show.all)) {
    pdf(paste0("12.all_pathway_netAnalysis/",pathways.show.all[j],".pdf"),height = 3,width = 5)
    netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all[j], width = 8, height = 2.5, font.size = 10)
    dev.off()
  }
  
  #Visualize the dominant senders (sources) and receivers (targets) in a 2D space
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  gg1 <- netAnalysis_signalingRole_scatter(cellchat)
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  # Signaling role analysis on the cell-cell communication networks of interest
  gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show.all)
  #> Signaling role analysis on the cell-cell communication network from user's input
  pdf("13.dominant_senders_receivers.pdf",width = 14)
  print(gg1 + gg2)
  dev.off()
  
  #Identify signals contributing most to outgoing or incoming signaling of certain cell groups
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
  
  pdf("14.out_in_coming_signaling.pdf",width = 14)
  print(ht1 + ht2)
  dev.off()
  
  #ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("PSAP", "NRG"))
  #Identify and visualize outgoing communication pattern of secreting cells
  library(NMF)
  #> Loading required package: pkgmaker
  #> Loading required package: registry
  #> Loading required package: rngtools
  #> Loading required package: cluster
  #> NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 15/16
  #>   To enable shared memory capabilities, try: install.extras('
  #> NMF
  #> ')
  #> 
  #> Attaching package: 'NMF'
  #> The following objects are masked from 'package:igraph':
  #> 
  #>     algorithm, compare
  library(ggalluvial)
  #selectK(cellchat, pattern = "outgoing")
  
  
  #Manifold and classification learning analysis of signaling networks
  cellchat <- computeNetSimilarity(cellchat, type = "functional")
  cellchat <- netEmbedding(cellchat, type = "functional")
  #> Manifold learning of the signaling networks for a single dataset
  cellchat <- netClustering(cellchat, type = "functional")
  #> Classification learning of the signaling networks for a single dataset
  # Visualization in 2D-space
  pdf("15.netcluster.pdf")
  netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
  dev.off()
  
  cellchat <- computeNetSimilarity(cellchat, type = "structural")
  cellchat <- netEmbedding(cellchat, type = "structural")
  #> Manifold learning of the signaling networks for a single dataset
  cellchat <- netClustering(cellchat, type = "structural")
  #> Classification learning of the signaling networks for a single dataset
  # Visualization in 2D-space
  pdf("15.netvisual.pdf")
  netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
  dev.off()
  
  saveRDS(cellchat,paste0(group[i],"_cellchat.rds"))
  
}

names(object.list) <- c("Control","Obesity","Diabetes")


cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat
#> An object of class CellChat created from a merged object with multiple datasets 
#>  555 signaling genes.
#>  7563 cells.


gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")

pdf("compare_interaction.pdf",width = 8)
print(gg1 + gg2)
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
p1 <- netVisual_diffInteraction(cellchat, weight.scale = T)
p2 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
