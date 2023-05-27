#' Extract genes in pathway
#' 
#' @description extract genes of selected pathway in GSEA result
#' 
#' @param select.pathway A vector contains pathways of interest. 
#' Use ID like "GO:0030001"
#' 
#' @param all.pathway A data.frame contains results of GSEA analysis. 
#' Column names must contain "ID" and "core_enrichment". 
#' Core_enrichment is the subset of genes that contributes most to the enrichment results
#' 
#' @param savefile A directory to save genes. use ".csv" as suffix
#' 
#' @return A vector of all selected genes

extract_genes_in_pathway <- function(select.pathway, all.pathway, savefile){
  
  require(clusterProfiler)
  require(tidyverse)
  
  POI <- filter(all.pathway, ID %in% select.pathway)
  
  # extract gene and transfer to SYMBOL
  
  gene <- unlist(strsplit(POI$core_enrichment, split = "/")) %>% unique() %>% 
    bitr(fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db") 
  
  colnames(gene)[2] <- "gene"
  write.csv(gene$gene, file = savefile)
  gene$gene
}


#' Preform default NicheNet analysis
#' 
#' @description 
#' preform NicheNet analysis by a given Seurat Object and target genes
#' genes are usually selected from some GO or KEGG pathways 
#' 
#' @param seu.obj A Seurat Object containing a column in meta.data "subtype" 
#' gives the information of receiver cell cluster and sender cell clusters
#' 
#' @param target.gene A vector containing genes of intereset
#' 
#' @param receiver A character of receiver cell cluster's name in "subtype"
#' 
#' @param sender A character or vector containing sender cell clusters's names in "subtype"
#' 
#' @param group.by A character of variables to group cells by;
#' 
#' @param ident.1 Cell class identity 1.
#' @param ident.2 Cell class identity 2.
#' 
#' @param assay_oi which expression matrix to use. default is "RNA"
#' 
#' @param ligand_target_matrix NicheNet’s ligand-target prior model
#' @param lr_network NicheNet’s ligand-target prior model
#' @param weighted_networks NicheNet’s ligand-target prior model
#' 
#' @return  A list containing main results of NicheNet analysis

preform_nn_analysis <- function(seu.obj
                                , target.gene
                                
                                , ligand_target_matrix
                                , lr_network
                                , weighted_networks
                                
                                , receiver
                                , sender
                                
                                , group.by
                                , ident.1
                                , ident.2
                                
                                , assay_oi = "RNA"
                                , logfc.threshold = 0.25
                                , test.use = "wilcox"
                                , min.pct = 0.1
                                ){
  ## load packages
  
  require(nichenetr)
  require(Seurat)
  require(tidyverse)
  
  # ============= define Receiver and Sender cell type =========================
  
  
  expressed_genes_receiver = get_expressed_genes(receiver
                                                 , seu.obj
                                                 , pct = 0.10
                                                 , assay_oi = assay_oi) # 
  background_expressed_genes = expressed_genes_receiver %>% 
    .[. %in% rownames(ligand_target_matrix)]
  
  
  sender_celltypes = sender
  list_expressed_genes_sender = sender_celltypes %>% 
    unique() %>% 
    lapply(get_expressed_genes, seu.obj, 0.10, assay_oi = assay_oi) # lapply to get the expressed genes of every sender cell type separately here
  expressed_genes_sender = list_expressed_genes_sender %>% 
    unlist() %>% 
    unique()
  
  
  # ================= Define a gene set of interest =====================
  # these are the genes in the “receiver/target” cell population that are \
  # potentially affected by ligands expressed by interacting cells 
  
  receiver.seu <- subset(seu.obj, subtype == receiver)
  receiver.deg <- FindMarkers(receiver.seu
                           , ident.1 = ident.1
                           , ident.2 = ident.2
                           , group.by = group.by
                           , logfc.threshold = logfc.threshold
                           , min.pct = min.pct
                           , test.use = test.use
                           ) %>% rownames()
  
  geneset_oi <- intersect(target.gene, receiver.deg)
  
  
  # ================ Define a set of potential ligands ====================
  # these are ligands that are expressed by the “sender/niche” cell population \
  # and bind a (putative) receptor expressed by the “receiver/target” population
  
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  
  potential_ligands = lr_network %>% 
    filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
    pull(from) %>% unique()
  
  
  # ================ Perform NicheNet ligand activity analysis ================
  # rank the potential ligands based on the presence of their target genes \
  # in the gene set of interest (compared to the background set of genes)
  
  ligand_activities = predict_ligand_activities(geneset = geneset_oi
                                                , background_expressed_genes = background_expressed_genes
                                                , ligand_target_matrix = ligand_target_matrix
                                                , potential_ligands = potential_ligands)
  
  ligand_activities = ligand_activities %>% arrange(-pearson) %>% 
    mutate(rank = rank(desc(pearson)))
  
  best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% 
    arrange(-pearson) %>% pull(test_ligand) %>% unique()
  
  
  # ===== Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis ===========
  
  active_ligand_target_links_df = best_upstream_ligands %>% 
    lapply(get_weighted_ligand_target_links
           , geneset = geneset_oi
           , ligand_target_matrix = ligand_target_matrix
           , n = 200) %>% bind_rows() %>% drop_na()
  
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df
                                                                   , ligand_target_matrix = ligand_target_matrix
                                                                   , cutoff = 0.33)
  
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% 
    rev() %>% make.names()
  
  order_targets = active_ligand_target_links_df$target %>% unique() %>% 
    intersect(rownames(active_ligand_target_links)) %>% make.names()
  
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  
  p_ligand_target_network = vis_ligand_target %>% 
    make_heatmap_ggplot("Prioritized ligands","Predicted target genes"
                        , color = "purple",legend_position = "top"
                        , x_axis_position = "top",legend_title = "Regulatory potential") + 
    theme(axis.text.x = element_text(face = "italic")) + 
    scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
  
  # ================= Receptors of top-ranked ligands ======================
  
  lr_network_top = lr_network %>% 
    filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% 
    distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  
  lr_network_top_df_large = weighted_networks_lr %>% 
    filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  lr_network_top_df = lr_network_top_df_large %>% 
    spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% 
    select(-to) %>% as.matrix() %>% 
    magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  
  p_ligand_receptor_network = vis_ligand_receptor_network %>% 
    t() %>% 
    make_heatmap_ggplot("Ligands","Receptors"
                        , color = "mediumvioletred"
                        , x_axis_position = "top"
                        , legend_title = "Prior interaction potential")
  
  nn.res <- list(ligand_activities = ligand_activities
                 , top_ligands = best_upstream_ligands
                 , active_ligand_target_links_df = active_ligand_target_links_df
                 , p_ligand_target_network = p_ligand_target_network
                 , vis_ligand_receptor_network = vis_ligand_receptor_network
                 , p_ligand_receptor_network = p_ligand_receptor_network)
  nn.res
}

