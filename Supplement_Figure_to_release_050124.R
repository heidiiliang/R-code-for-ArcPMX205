library(Matrix)
library(stringr)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(scales))
library(DESeq2)
library(ggforce)
library(RColorBrewer)
library(Seurat)
options("stringsAsFactors" = FALSE)
library(ComplexHeatmap)
library(ggplot2)
library(devtools::install_github("sqjin/CellChat"))
library(patchwork)
library(NMF)
library(ggalluvial)
library(CellChat)

##################################
########## FIGURE S2 #############
##################################
ad.combined.annotate <- readRDS("sc_sn_postSCT_20kanchors_RPCA_res1_annotated_050522.rds")
####Figure S2A
palette7 <- colorRampPalette(c( "#49796b","#409cb3", "#c6588f", "#928678", "#a92a3f", "#bc8400","#a99fbf"))(n=7)
ad.combined.annotate <-subset(ad.combined.annotate, subset = cellType == "Mixed", invert = T)
all_RNA_counts_vln <- VlnPlot(object=ad.combined.annotate, features = c("C1qa","C1qb","C1qc","C2","C4a","C4b","C3","Hc"), assay ="RNA", slot="counts", group.by = 'cellType', cols = palette7)
all_RNA_counts_vln

####Figure S2B
mic_sub <- readRDS(file = "Splitseq_micsub_042423.rds")
ms_clusters_color <- colorRampPalette(c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99"))(n=11)
micsub_RNA_counts_vln <- VlnPlot(object=mic_sub, features = c('Selplg', 'Csf1r', 'Tgfbr1', 'Hexb','Cx3cr1','P2ry12'), assay ="RNA", slot="counts", group.by = 'seurat_clusters', cols = ms_clusters_color)
micsub_RNA_counts_vln

####Figure S2C
ast <- readRDS("Splitseq_ast_042423.rds")
ast_clusters_color <- colorRampPalette(c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928"))(n=12)
ast_RNA_counts_vln <- VlnPlot(object=ast, features = c('Slc1a2','Slc1a3','Aqp4','Apoe','Clu'), assay ="RNA", slot="counts", group.by = 'seurat_clusters', cols = ast_clusters_color)
ast_RNA_counts_vln

##################################
########## FIGURE S3 #############
##################################
####Figure S3
mic_sub <- readRDS(file = "Splitseq_micsub_042423.rds")
ms_clusters_color <- colorRampPalette(c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99"))(n=11)
microglia_markers <-c("Spi1",'Mef2c','Tmem119', 'Selplg', 'Csf1r', 'Tgfbr1','Tgfbr2', 'Hexb','Cx3cr1','P2ry12'
                      , 'Clec7a','Tyrobp', 'Apoe', 'Trem2', 'Fth1','B2m', 'C1qa','C1qb','C1qc','Lyz2' #DAM 1
                      , 'Lpl', 'Cst7','Fgf13','Spp1', 'Csf1', 'Axl', 'Itgax', 'Cd9', #DAM 2
                      'C5ar1',
                      "Cd44","Cdk8","Gphn",#Cluster6
                      "Inpp5d","Mertk") #Cluster9

# average expression across clusters
DefaultAssay(mic_sub) = "integrated"
matrix = AverageExpression(mic_sub, group.by = "exp_gen_treat_sex", return.seurat = F)
integrated_matrix = matrix$integrated # get integrated average expression, other options are raw RNA and normalized SCT
microglia_markers = microglia_markers[microglia_markers %in% rownames(integrated_matrix)] # subset by genes of interest
col.order = as.character(c("M-WT-H2O-F","M-WT-H2O-M","M-WT-PMX-F","M-WT-PMX-M",
                           "N-WT-H2O-F","N-WT-H2O-M","N-WT-PMX-F","N-WT-PMX-M",
                           "M-Arc-H2O-F","M-Arc-H2O-M","M-Arc-PMX-F","M-Arc-PMX-M",
                           "N-Arc-H2O-F","N-Arc-H2O-M","N-Arc-PMX-F","N-Arc-PMX-M"))
integrated_matrix = integrated_matrix[microglia_markers,col.order] # make sure rows and columns are in order using gene list and col.order
integrated_matrix = t(apply(integrated_matrix, 1, function(x) rescale(x, to = c(-2, 2)))) # normalize matrix between -2 and 2
cluster_celltype = unique(mic_sub@meta.data[,c("exp_gen_treat_sex","cellType")]) # get cluster-celltype table
cluster_celltype = cluster_celltype[match(col.order,cluster_celltype$seurat_clusters),] # put table in same order as matrix
cluster_celltype$cellType = factor(cluster_celltype$cellType, levels = "Microglia")

# make top color bar annotation                                                                   
column_ha = HeatmapAnnotation( 
  exp_gen_treat_sex = cluster_celltype$exp_gen_treat_sex,
  col = list(exp_gen_treat_sex = c("1"="#a6cee3","2"="#1f78b4","3"="#b2df8a","4"="#33a02c","5"="#fb9a99","6"="#e31a1c","7"="#fdbf6f","8"="#ff7f00","9"="#cab2d6","10"="#6a3d9a","11"="#ffff99")))
# complex heatmap
micsub_expgentreatsex_all_heatmap <-Heatmap(integrated_matrix, name = "Expression", show_row_dend = F, show_column_dend = T, 
                                            cluster_columns = F, cluster_rows = F, col = viridis(100), column_names_side = "top",
                                            top_annotation = column_ha)
micsub_expgentreatsex_all_heatmap

##################################
########## FIGURE S7 #############
##################################
Arc.PMX.cellchat <- readRDS("ArcPMX_cellchat_trim02_mincells10_nosubsetDB_100422.rds")
Arc.V.cellchat <- readRDS("ArcV_cellchat_trim02_mincells10_nosubsetDB_100422.rds")
WT.PMX.cellchat <- readRDS("WTPMX_cellchat_trim02_mincells10_nosubsetDB_100322.rds")
WT.V.cellchat <- readRDS("WTV_cellchat_trim02_mincells10_nosubsetDB_100322.rds")

####Figure S7A
netVisual_circle(WT.V.cellchat@net$count, weight.scale = T, label.edge= T,
                 edge.label.cex = 1.2,arrow.size =0.1, arrow.width = 1, vertex.label.cex = 1,vertex.size.max = 20,
                 color.use = c( "#49796b" , "#409cb3" , "#c6588f", "#e4994b", "#a92a3f", "#bc8400" , "#a99fbf"),
                 title.name = "WT.Veh Number of interactions")

####Figure S7B
netVisual_circle(WT.PMX.cellchat@net$count, weight.scale = T, label.edge= T,
                 edge.label.cex = 1.2,arrow.size =0.1, arrow.width = 1, vertex.label.cex = 1,vertex.size.max = 20,
                 color.use = c( "#49796b" , "#409cb3" , "#c6588f", "#e4994b", "#a92a3f", "#bc8400" , "#a99fbf"),
                 title.name = "WT.PMX Number of interactions")

####Figure S7C
netVisual_circle(Arc.V.cellchat@net$count, weight.scale = T, label.edge= T,
                 edge.label.cex = 1.2,arrow.size =0.1, arrow.width = 1, vertex.label.cex = 1,vertex.size.max = 20,
                 color.use = c( "#49796b" , "#409cb3" , "#c6588f", "#e4994b", "#a92a3f", "#bc8400" , "#a99fbf"),
                 title.name = "Arc.Veh Number of interactions")

####Figure S7D
netVisual_circle(Arc.PMX.cellchat@net$count, weight.scale = T, label.edge= T,
                 edge.label.cex = 1.2,arrow.size =0.1, arrow.width = 1, vertex.label.cex = 1,vertex.size.max = 20,
                 color.use = c( "#49796b" , "#409cb3" , "#c6588f", "#e4994b", "#a92a3f", "#bc8400" , "#a99fbf"),
                 title.name = "Arc.PMX Number of interactions")


####Figure S7E, F
WT.PMX.cellchat_a1 <- netAnalysis_computeCentrality(WT.PMX.cellchat, slot.name = "netP")
WT.V.cellchat_a1 <- netAnalysis_computeCentrality(WT.V.cellchat, slot.name = "netP")
object.list <- list(WTVeh = WT.V.cellchat_a1,WTPMX = WT.PMX.cellchat_a1)
WT.cellchat <- mergeCellChat(object.list, add.names = names(object.list))

WT.V.cellchat_gg1 <- netAnalysis_signalingRole_scatter(WT.V.cellchat_a1, title = "WT.Vehicle")
WT.PMX.cellchat_gg1 <- netAnalysis_signalingRole_scatter(WT.PMX.cellchat_a1, title = "WT.PMX")
Veh_num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(Veh_num.link), max(Veh_num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, dot.size = c(1, 6),
                                               color.use = c( "#49796b", "#409cb3", "#c6588f", "#e4994b", "#a92a3f", "#bc8400", "#a99fbf"))
}
patchwork::wrap_plots(plots = gg)

####Figure S7G, H
Arc.PMX.cellchat_a1 <- netAnalysis_computeCentrality(Arc.PMX.cellchat, slot.name = "netP")
Arc.V.cellchat_a1 <- netAnalysis_computeCentrality(Arc.V.cellchat, slot.name = "netP")
object.list <- list(ArcVeh = Arc.V.cellchat_a1, ArcPMX = Arc.PMX.cellchat_a1)
Arc.cellchat <- mergeCellChat(object.list, add.names = names(object.list))

Arc.PMX.cellchat_gg1 <- netAnalysis_signalingRole_scatter(Arc.PMX.cellchat_a1, title = "Arc.PMX205")
Arc.V.cellchat_gg1 <- netAnalysis_signalingRole_scatter(Arc.V.cellchat_a1, title = "Arc.Vehicle")
Arc_num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(Arc_num.link), max(Arc_num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, dot.size = c(1, 6),
                                               color.use = c( "#49796b", "#409cb3", "#c6588f", "#e4994b", "#a92a3f", "#bc8400", "#a99fbf"))
}
patchwork::wrap_plots(plots = gg)

####Figure S7I
Arc.V.cellchat_a1 <- netAnalysis_computeCentrality(Arc.V.cellchat, slot.name = "netP")
WT.V.cellchat_a1 <- netAnalysis_computeCentrality(WT.V.cellchat, slot.name = "netP")
object.list <- list(WTVeh = WT.V.cellchat_a1, ArcVeh = Arc.V.cellchat_a1)
Veh.cellchat <- mergeCellChat(object.list, add.names = names(object.list))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(Veh.cellchat, 
                          weight.scale = T, label.edge = T,
                          edge.label.cex = 1.2,arrow.size =0.1, arrow.width = 1, vertex.label.cex = 1,vertex.size.max = 20,
                          color.use = c( "#49796b", "#409cb3", "#c6588f", "#e4994b", "#a92a3f", "#bc8400", "#a99fbf"),
                          title.name = "ArcVeh vs WTVeh Differential number of interactions")
netVisual_diffInteraction(Veh.cellchat, weight.scale = T, measure = "weight",label.edge = T, 
                          edge.label.cex = 1.2,arrow.size =0.1, arrow.width = 1, vertex.label.cex = 1,vertex.size.max = 20,
                          color.use = c( "#49796b", "#409cb3", "#c6588f", "#e4994b", "#a92a3f", "#bc8400", "#a99fbf"),
                          title.name = "ArcVeh vs WTVeh Differential interaction Strength")

####Figure S7J
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(Arc.cellchat, 
                          weight.scale = T, label.edge = T,
                          edge.label.cex = 1.2,arrow.size =0.1, arrow.width = 1, vertex.label.cex = 1,vertex.size.max = 20,
                          color.use = c( "#49796b", "#409cb3", "#c6588f", "#e4994b", "#a92a3f", "#bc8400", "#a99fbf"),
                          title.name = "ArcPMX vs ArcVeh Differential number of interactions")
netVisual_diffInteraction(Arc.cellchat, weight.scale = T, measure = "weight",label.edge = T, 
                          edge.label.cex = 1.2,arrow.size =0.1, arrow.width = 1, vertex.label.cex = 1,vertex.size.max = 20,
                          color.use = c( "#49796b", "#409cb3", "#c6588f", "#e4994b", "#a92a3f", "#bc8400", "#a99fbf"),
                          title.name = "ArcPMX vs ArcVeh Differential interaction Strength")

#Differential number of interactions or interaction strength among different cell populations
#The differential number of interactions or interaction strength in the cell-cell communication network between two datasets can be visualized using circle plot, where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
palette8 <- colorRampPalette(c( "#49796b" , "#409cb3", "#c6588f", "#e4994b", "#a92a3f", "#bc8400", "#928678","#a99fbf"))(n=8)
#"Neurons" = "#49796b", "Astrocytes" = "#409cb3", "Oligodendrocytes" = "#c6588f","OPCs" = "#a92a3f", "Endothelial" = "#bc8400", "Pericytes" = "#a99fbf","Microglia" = "#e4994b","Mixed" = "#928678")))

netVisual_diffInteraction(Arc.cellchat, 
                          weight.scale = T, label.edge = T,
                          edge.label.cex = 1.2,arrow.size =0.1, arrow.width = 1, vertex.label.cex = 1,vertex.size.max = 20,
                          color.use = c( "#49796b", "#409cb3", "#c6588f", "#e4994b", "#a92a3f", "#bc8400", "#a99fbf"),
                          title.name = "ArcPMX vs ArcVeh Differential number of interactions")

####Figure S8A, B
# relative information flow and information flow for ArcVeh and ArcPMX
Arc_r1 <- rankNet(Arc.cellchat, mode = "comparison", stacked = T, do.stat = TRUE)

####Figure S8B
Arc_r2 <- rankNet(Arc.cellchat, mode = "comparison", stacked = F, do.stat = TRUE)

####Figure S9
#> Comparing communications on a merged object
Arc_netvisualbubble1 <- netVisual_bubble(Arc.cellchat, sources.use = c(1:7), targets.use = c(1:7),  comparison = c(1, 2), max.dataset = 2, 
                                         signaling=c("CLDN","CDH","SEMA5","FGF","MPZ","CDH5","SEMA7","TGFb","BMP", #significant in ArcPMX
                                                     "PECAM1","TENASCIN","FN1","EPHA","COLLAGEN","PROS","AGRN","LAMININ","ANGPT","CSF","MAG","VTN"),#significant in ArcVeh
                                         title.name = "Increased signaling in ArcPMX (ArcPMX vs ArcVeh)", color.heatmap = "viridis",
                                         angle.x = 45, remove.isolate = T, show.legend = TRUE, font.size = 13, font.size.title = 24)
Arc_netvisualbubble1 