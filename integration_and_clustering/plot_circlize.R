rm(list = ls())

library(circlize)
library(tidyverse)
library(patchwork)

# plot1cell function ===============================

transform_coordinates <- function(
  coord_data, 
  zoom = 0.8
){
  center_data<-coord_data-mean(c(min(coord_data),max(coord_data)))
  max_data<-max(center_data)
  new_data<-center_data * zoom / max_data
  new_data
}

get_segment <- function(
  dat, 
  group
){
  dat<-dat[order(dat[,group],decreasing = F), ]
  rownames(dat)<-1:nrow(dat)
  dat<-dat[!duplicated(dat[,group]),]
  dat_seg<-as.integer(rownames(dat))
  dat_seg
}

# Setting colors
ClsName <- c("EX", "IN", "OXT", "HDC", "TANY", "Ependy", "AST", "ENDO", "MIC", "OLI", "OPC")
ClsColor <- c("#0073a8", "#434da2", "#00afcc", "#00a497", "#d70035", "#ea5549", "#00984f"
              , "#d9e367", "#7f1184", "#fcc800", "#e5a323")

NucName <- c("PVN", "DMH", "LHA", "VMH", "INF")
NucColor <- rev(c("#ea5550", "#fdd35c", "#00ac97", "#00a1e9", "#4d4398"))

groupName <- c("Control", "Obesity", "Diabetes")
groupColor <- c("#2980b9", "#f39c12", "#c0392b")

SubtypeColor <- c()

# initiation =======================================

AllCell.meta <- readRDS("data/all_umap0.01_meta.rds")
AllCell.meta <- AllCell.meta %>% select(-contains("pANN"))

table(AllCell.meta$celltype)
table(AllCell.meta$subtype)

AllCell.meta$celltype <- as.character(AllCell.meta$celltype)
AllCell.meta[AllCell.meta$subtype == "HDC" | AllCell.meta$subtype == "HDC/NEFL", "celltype"] <- "HDC"
AllCell.meta$celltype <- factor(AllCell.meta$celltype, levels = c("EX", "IN", "OXT", "HDC", "TANY", "Ependy"
                                                                  , "AST", "ENDO", "MIC", "OLI", "OPC"))

# get embedding coordinate

AllCell.umap <- read.csv("week2/allcell_umap_embedding.csv"
                         , row.names = 1)

AllCell.tsne <- read.csv("week2/allcell_tsne_embedding.csv"
                         , row.names = 1)

AllCell.meta <- cbind(AllCell.meta, AllCell.tsne)
AllCell.meta <- cbind(AllCell.meta, AllCell.umap)

AllCell.meta$p_UMAP1 <- transform_coordinates(AllCell.meta$UMAP1)
AllCell.meta$p_UMAP2 <- transform_coordinates(AllCell.meta$UMAP2)

cellnames<-rownames(AllCell.meta)
AllCell.meta$cells<-rownames(AllCell.meta)

# order cells

celltypes <- c("EX", "IN", "OXT", "HDC", "TANY", "Ependy", "AST", "ENDO", "MIC", "OLI", "OPC")
AllCell.meta$celltype<-as.character(AllCell.meta$celltype)

new_dat <- list()
for (i in 1:length(celltypes)){
  dat1<-AllCell.meta[AllCell.meta$celltype == celltypes[i],]
  dat1$x_polar<-1:nrow(dat1)
  new_dat[[i]]<-dat1
}
new_dat<-do.call('rbind', new_dat)
new_dat$x_polar2 <- log10(new_dat$x_polar)
AllCell.meta <- new_dat

# Plot circlize
pdf("week3/plot_circlize.pdf", width = 12, height = 12)
circos.clear()
par(bg = "#FFFFFF")
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0.01,0),"track.height" = 0.01, gap.degree =c(rep(2, (length(celltypes)-1)),12),points.overflow.warning=FALSE)
circos.initialize(sectors =  AllCell.meta$celltype, x = AllCell.meta$x_polar2)

circos.track(AllCell.meta$celltype, AllCell.meta$x_polar2, y=AllCell.meta$UMAP2, bg.border=NA,panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter,
              CELL_META$cell.ylim[2]+ mm_y(4),
              CELL_META$sector.index,
              cex=1.5, col = 'black', facing = "bending.inside", niceFacing = T)
  circos.axis(labels.cex = 0.5, col = 'black', labels.col =  'black')
})

for(i in 1:length(celltypes)){
  dd<-AllCell.meta[AllCell.meta$celltype==celltypes[i],]
  circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), y1 = 0, col = ClsColor[i],  lwd=3, sector.index = celltypes[i])
}
text(x = 1, y=0.1, labels = "Cluster", cex = 1, col = 'black',srt=-90)

# # add subtype track
# circos.track(AllCell.meta$celltype, AllCell.meta$x_polar2, y=AllCell.meta$UMAP2, bg.border=NA)
# for(i in 1:length(celltypes)) {
#   AllCell.meta_cl<-AllCell.meta[AllCell.meta$celltype==celltypes[i],]
#   dat_seg<-get_segment(AllCell.meta_cl, group = "subtype")
#   dat_seg2<-c(dat_seg[-1]-1, nrow(AllCell.meta_cl))
#   scale_factor<-max(AllCell.meta_cl$x_polar2)/nrow(AllCell.meta_cl)
#   dat_seg<-scale_factor*dat_seg
#   dat_seg2<-scale_factor*dat_seg2
#   circos.segments(x0 = dat_seg, y0 = 0, x1 = dat_seg2, y1 = 0, col = NucColor, sector.index = celltypes[i], lwd=3)
# }
# text(x = (1-0.03*(2-1)), y=0.1, labels = "Subtype", cex = 0.4, col = 'black',srt=-90)
# 

# add nuclei track
circos.track(AllCell.meta$celltype, AllCell.meta$x_polar2, y=AllCell.meta$UMAP2, bg.border=NA)
for(i in 1:length(celltypes)) {
  AllCell.meta_cl<-AllCell.meta[AllCell.meta$celltype==celltypes[i],]
  dat_seg<-get_segment(AllCell.meta_cl, group = "nuclei")
  dat_seg2<-c(dat_seg[-1]-1, nrow(AllCell.meta_cl))
  scale_factor<-max(AllCell.meta_cl$x_polar2)/nrow(AllCell.meta_cl)
  dat_seg<-scale_factor*dat_seg
  dat_seg2<-scale_factor*dat_seg2
  circos.segments(x0 = dat_seg, y0 = 0, x1 = dat_seg2, y1 = 0, col = NucColor, sector.index = celltypes[i], lwd=3)
}
text(x = (1-0.03*(2-1)), y=0.1, labels = "Nuclei", cex = 1, col = 'black',srt=-90)


# add group track
circos.track(AllCell.meta$celltype, AllCell.meta$x_polar2, y=AllCell.meta$UMAP2, bg.border=NA)
for(i in 1:length(celltypes)) {
  AllCell.meta_cl<-AllCell.meta[AllCell.meta$celltype==celltypes[i],]
  dat_seg<-get_segment(AllCell.meta_cl, group = "diabetes")
  dat_seg2<-c(dat_seg[-1]-1, nrow(AllCell.meta_cl))
  scale_factor<-max(AllCell.meta_cl$x_polar2)/nrow(AllCell.meta_cl)
  dat_seg<-scale_factor*dat_seg
  dat_seg2<-scale_factor*dat_seg2
  circos.segments(x0 = dat_seg, y0 = 0, x1 = dat_seg2, y1 = 0, col = groupColor, sector.index = celltypes[i], lwd=3)
}
text(x = (1-0.03*(3-1)), y=0.1, labels = "Group", cex = 1, col = 'black',srt=-90)

dev.off()

