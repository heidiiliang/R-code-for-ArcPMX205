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
library(randomcoloR)

#######
#Identify clusters:
#The top 50 most variable genes included microglial identity genes (Crybb1, Alox5ap, Maf), 
#Hif1a and its target genes (Hif1a, Igf1, Spp1, Pkm, Gapdh), α-defensin genes (Defa20, Defa21, Defa24, Gm15284), 
#chemokines (Ccl3, Ccl4, Ccl6) and lysosomal genes (Lyz2, Cst7, Ctsa)
#https://www.nature.com/articles/s41467-021-23111-1#Sec2
#AQP4 in astrocytes28, APOLD1 in endothelial cells29, CCL4 in microglia30, RELN in neurons31, PLP1 in oligodendrocytes32, and PDGFRA in OPCs33
#https://www.nature.com/articles/s41598-018-27293-5#MOESM1
#Astrocytes: Aqp4, Gfap, Slc1a2, Gja1, Gjb6, Clu, Aldoc, Slc1a3, Slc1a2
#Endothelial:Apold1, Ftl1, Cd34, Vwf, Abcb1, Rgs5, Apold1, Ifitm1, Itm2a, Bsg
#Olygodendrocytes: Plp1, Mbp, Mag, Mog, Cnp, Cldn11, Cnp, Mal, Mog, Cryab
#Microglia: Cxcr3, Itgam, Ccl3, Ccl4, Csf1r, P2ry12, C1qb, Ctss, Tyrobp, Cd83, Crybb1, Alox5ap, Maf
##Neuron: Rbfox3, Npy, Bi, Scg2, Gad2, Reln, Npy, Snap25, 
#Perycytes: Pdgfrb, Cd248, Rgs5, Foxf2, Cd19
#ependymal:
#https://www.panglaodb.se/markers.html?cell_type=%27Ependymal%20cells%27
#OPCs: Pdgfra, Vcan, Gpr17, Sox10, Olig1, Olig2, Tnr, Cntn1

#create a pallete with random distinct colors
n <- 42
palette42 <- distinctColorPalette(n)
#visualize the colors in a pie plot
pie(rep(1, n), col=palette42)
print(palette42) 

#colors for 42 clusters
palette42 <- colorRampPalette(c("#88E766", "#8C64E6", "#77F2E4", "#D673DE", "#76A1DA", "#572878", "#A7E490", "#DB5061", "#BEEAE5" 
                                , "#B58CD9", "#E48443", "#A7BAA8", "#D6B6E4", "#91A63F", "#E853A7", "#fcafac", "#DBEA8C", "#62DAE4" 
                                , "#7738E6", "#E3C070", "#cc324c", "#E78CC3", "#4DBF57", "#E2B4C2", "#66A275", "#60EC3B", "#53C4AE" 
                                , "#DEB798", "#71C4E5", "#71A2AA", "#A7EEC9", "#748BE6", "#E2E9BB", "#61E3A2", "#EBDDD7", "#C0CEE3"
                                , "#D1F04A", "#E8D03F", "#E28F83", "#89828E", "#87683D", "#A06280"))(n=42) 

# Read combined object
ad.combined <- readRDS('sc_sn_postSCT_20kanchors_RPCA_res1.rds')
levels(ad.combined)

DimPlot(ad.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE,repel = TRUE, cols = palette42)
ad.combined.markers <- FindAllMarkers(ad.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(ad.combined.markers, "ad.combined_clusters_allgene_081122.txt",sep="\t",row.names = F)

ad.combined.markers <- FindAllMarkers(ad.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ad.combined.markers %>%
  group_by(cluster) %>%
  top_n(n =20, wt = avg_log2FC) -> ad.combined.markers_20

ad.combined.markers_20

DoHeatmap(ad.combined.markers,features = ad.combined.markers_20$gene,draw.lines = T,lines.width = 1 
          ,group.by = 'ident') + #NoLegend() + group.colors = new_palette_micro
  scale_fill_viridis( na.value = "white") +
  theme(text = element_text(size = 20, color = 'black'))

write.table(ad.combined.markers_20, "ad.combined_clusters_top20_120122.txt",sep="\t",row.names = F)

ad.combined.sample <- as.data.frame(table(ad.combined@meta.data$sample))

Idents(ad.combined)
levels(ad.combined)
ad.combined <- RenameIdents(object = ad.combined, '0'="Neurons", '1'="Neurons", '2'="Neurons", '3'="Neurons", '4'= "Astrocytes", '5'="Neurons"
                            , '6'="Oligodendrocytes", '7'="Oligodendrocytes", '8'="Neurons", '9'="Neurons"
                            , '10'="Astrocytes", '11'="Neurons", '12'="Neurons"
                            ,'13'="Microglia", '14'="Neurons", '15'="Microglia", '16'="Microglia", '17'="Neurons", '18'="OPCs", '19'="Endothelial"
                            ,'20'="Neurons", '21'="Neurons", '22'="Neurons", '23'="Microglia", '24'="Neurons", '25'="Mixed", '26'="Microglia"
                            ,'27'="Neurons", '28'="Neurons", '29'="Neurons", '30'="Pericytes", '31'="Neurons", '32'="Neurons", '33'="Mixed"
                            ,'34'="Neurons", '35'="Neurons", '36'="Neurons", '37'="OPCs", '38'="Astrocytes", '39'="Neurons", '40'="Neurons", '41'="Neurons")
levels(ad.combined)
Idents(ad.combined)
ad.combined[["cellType"]] <- Idents(object = ad.combined)
ad.combined@meta.data$cl_cellType = paste0(ad.combined@meta.data$cellType, "_", ad.combined@meta.data$seurat_clusters)
saveRDS(ad.combined, file ="sc_sn_postSCT_20kanchors_RPCA_res1_annotated_050522.rds")




ad.combined.annotate <- readRDS("sc_sn_postSCT_20kanchors_RPCA_res1_annotated_050522.rds")
DefaultAssay(ad.combined.annotate)

cellType_levels <- c(0,1,2,3,5,8,9,11,12,14,17,20,21,22,24,27,28,29,31,32,34,35,36,39,40,41,
                     4,10,38,
                     6,7,18,37,
                     19,
                     30,
                     13,15,16,23,26,
                     25,33)

palette8 <- colorRampPalette(c( "#49796b", "#409cb3", "#c6588f", "#a92a3f", "#bc8400", "#a99fbf", "#e4994b","#928678"))(n=8)
#"Neurons" = "#49796b", "Astrocytes" = "#409cb3", "Oligodendrocytes" = "#c6588f",#"OPCs" = "#a92a3f", "Endothelial" = "#bc8400", "Pericytes" = "#a99fbf",#"Microglia" = "#e4994b","Mixed" = "#928678")))

####Figure 1B
all_dimplot_clusters <- DimPlot(ad.combined.annotate, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, cols = palette42)

####Figure 1C
all_dimplt_exp <-DimPlot(ad.combined.annotate, reduction = "umap", group.by = "experiment",label = FALSE,repel = TRUE, label.size = 8, cols=c("#000080","#F19CBB"))

####Figure 1D
ad.combined.annotate_celltype <- ad.combined.annotate
ad.combined.annotate_celltype@meta.data$seurat_clusters <- factor(ad.combined.annotate_celltype@meta.data$seurat_clusters, # Change ordering manually
                                                                  levels = c("0","1","2","3","5","8","9","11","12","14","17","20","21","22","24","27","28","29","31","32","34","35","36","39","40","41",
                                                                             "4","10","38",
                                                                             "6","7","18","37",
                                                                             "19",
                                                                             "30",
                                                                             "13","15","16","23","26",
                                                                             "25","33"))

all_barplot_exp <- ggplot(ad.combined.annotate_celltype@meta.data, aes(fill=experiment, x=seurat_clusters, width=.5)) + geom_bar(position = "fill")+ 
  scale_fill_manual(values = c("#000080","#F19CBB"),labels = c("SC","SN")) + 
  theme_bw() + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 16))+ 
  theme(legend.position="right", legend.title =element_text(size=16),legend.text = element_text(size=14))+
  theme(axis.text = element_text(size = 18), axis.text.x = element_text(size = 14, color = 'black'), axis.text.y = element_text(size = 14, color = 'black')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank()) +
  ggtitle("SN+SC Experiment Porportion")

all_barplot_exp

####Figure 1E
all_barplot_gentreat <- ggplot(ad.combined.annotate_celltype@meta.data, aes(fill=gen_treat, x=seurat_clusters, width=.5)) + geom_bar(position = "fill")+ 
  theme_bw() + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 16))+ 
  theme(legend.position="right", legend.title =element_text(size=16),legend.text = element_text(size=14))+
  theme(axis.text = element_text(size = 18), axis.text.x = element_text(size = 14, color = 'black'), axis.text.y = element_text(size = 14, color = 'black')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank()) +
  ggtitle("SN+SC Treatment Porportion")

####Figure1F
ad_markers_complement <- c("C1qa","C1qb","C1qc","C1r","C1s",
                           "C2","Cfb","C4a","C4b","C3","C5",
                           "C6","C7","C8a","C8b","C8g","C9",
                           "Cfd","Cfp","Mbl2","Masp1","Masp2",
                           "Fcn1","Fcn2","Fcn3",
                           "C4bpa","C4bpb",
                           "Cr1","Cr2",
                           "Cd55", "Cd46",
                           "Cfh","Cfhr1","Cfhr2","Cfhr3","Cfhr4","Cfhr5","Hc",
                           "Csmd1","Csmd2","Csmd3",
                           "Cd59","Cfi","Serping1","Masp1","Masp2",
                           "Itgam","Itgb2","Itgax","C5ar1","C5ar2","C3ar1", "Cd93")

# average expression across clusters
matrix = AverageExpression(ad.combined.annotate, group.by = "seurat_clusters", return.seurat = F)
integrated_matrix = matrix$integrated # get integrated average expression, other options are raw RNA and normalized SCT
ad_markers_complement = ad_markers_complement[ad_markers_complement %in% rownames(integrated_matrix)] # subset by genes of interest
# set column order
col.order = as.character(c(0,1,2,3,5,8,9,11,12,14,17,20,21,22,24,27,28,29,31,32,34,35,36,39,40,41,
                           4,10,38,
                           6,7,18,37,
                           19,
                           30,
                           13,15,16,23,26,
                           25,33))
integrated_matrix = integrated_matrix[ad_markers_complement,col.order] # make sure rows and columns are in order using gene list and col.order
integrated_matrix = t(apply(integrated_matrix, 1, function(x) rescale(x, to = c(-2, 2)))) # normalize matrix between -2 and 2
cluster_celltype = unique(ad.combined.annotate@meta.data[,c("seurat_clusters","cellType")]) # get cluster-celltype table
cluster_celltype = cluster_celltype[match(col.order,cluster_celltype$seurat_clusters),] # put table in same order as matrix
cluster_celltype$cellType = factor(cluster_celltype$cellType, levels = c("Neurons","Astrocytes","Oligodendrocytes","OPCs",
                                                                         "Endothelial","Pericytes","Microglia","Mixed"))
# make top color bar annotation                                                                   
column_ha = HeatmapAnnotation( 
  cellType = cluster_celltype$cellType,
  col = list(cellType = c("Neurons" = "#49796b", "Astrocytes" = "#409cb3", "Oligodendrocytes" = "#c6588f",
                          "OPCs" = "#a92a3f", "Endothelial" = "#bc8400", "Pericytes" = "#a99fbf",
                          "Microglia" = "#e4994b","Mixed" = "#928678")))


# complex heatmap
complement_all_heatmap <-Heatmap(integrated_matrix, name = "Expression",
                                 show_row_dend = F, 
                                 show_column_dend = F, 
                                 cluster_columns = F,
                                 cluster_rows = F, 
                                 col = viridis(100),
                                 column_names_side = "top",
                                 top_annotation = column_ha)
complement_all_heatmap

##################################
########## FIGURE 2 ##############
##################################

####Figure 2A
all_dimplot <- DimPlot(ad.combined.annotate, reduction = "umap", group.by = "cellType", cols=palette8)

####Figure 2B
sn <- subset(ad.combined.annotate, subset = experiment == "sn")

sn_barplot <- ggplot(sn@meta.data, aes(fill=cellType, x=group, width=.5)) + geom_bar(position = "fill")+ scale_fill_manual(values = palette8) +
  coord_flip() + theme_bw() + 
  theme(plot.title = element_text(size = 20), axis.title = element_text(size = 16))+ 
  theme(legend.position="bottom", legend.title =element_text(size=14),legend.text = element_text(size=12))+
  theme(axis.text = element_text(size = 18), axis.text.x = element_text(size = 14, color = 'black'), axis.text.y = element_text(size = 14, color = 'black')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + 
  ggtitle("SN-only CellType Porportion")# change the order of the factor levels
gentreat_levels <- c("WT.Vehicle", "WT.PMX205", "Arc.Vehicle", "Arc.PMX205")
levels(sn_barplot$data$cellType)
sn_barplot$data$group <- factor(x = sn_barplot$data$gen_treat, levels = gentreat_levels)
sn_barplot

####Figure 2C
ncol(x=sn) #49058 nuclei
sn_sub = subset(sn, subset = cellType %in% c("Neurons","Astrocytes","Oligodendrocytes","OPCs","Endothelial", "Microglia"))
ncol(x=sn_sub) #47734 nuclei

DefaultAssay(sn_sub) #intergrated

# average expression across clusters
sn_sub$celltype_gentreat = paste0(sn_sub$cellType,"_",sn_sub$gen_treat)
matrix = AverageExpression(sn_sub, group.by = "celltype_gentreat", return.seurat = F)
integrated_matrix = matrix$integrated

ad_markers_complement = c("C1qa","C1qb","C1qc","C1r","C1s", 
                          "C2","Cfb","C4a","C4b","C3","C5",
                          "C6","C7","C8a","C8b","C8g","C9", 
                          "Cfd","Cfp","Mbl2","Masp1","Masp2",
                          "Fcn1","Fcn2","Fcn3",
                          "C4bpa","C4bpb",
                          "Cr1","Cr2",
                          "Cd55", "Cd46",
                          "Cfh","Cfhr1","Cfhr2","Cfhr3","Cfhr4","Cfhr5","Hc",
                          "Csmd1","Csmd2","Csmd3",
                          "Cd59","Cfi","Serping1","Masp1","Masp2",
                          "Itgam","Itgb2","Itgax","C5ar1","C5ar2","C3ar1", "Cd93",
                          "C5a","Hc")
ad_markers_complement = ad_markers_complement[ad_markers_complement %in% rownames(integrated_matrix)] # subset by genes of interest

# set column order
col.order = as.character(c("Neurons_WT.Vehicle","Neurons_WT.PMX205","Neurons_Arc.Vehicle","Neurons_Arc.PMX205",
                           "Astrocytes_WT.Vehicle","Astrocytes_WT.PMX205","Astrocytes_Arc.Vehicle","Astrocytes_Arc.PMX205",
                           "Oligodendrocytes_WT.Vehicle","Oligodendrocytes_WT.PMX205","Oligodendrocytes_Arc.Vehicle","Oligodendrocytes_Arc.PMX205",
                           "OPCs_WT.Vehicle","OPCs_WT.PMX205","OPCs_Arc.Vehicle","OPCs_Arc.PMX205",
                           "Endothelial_WT.Vehicle","Endothelial_WT.PMX205","Endothelial_Arc.Vehicle","Endothelial_Arc.PMX205",
                           "Microglia_WT.Vehicle","Microglia_WT.PMX205","Microglia_Arc.Vehicle","Microglia_Arc.PMX205"))
integrated_matrix = integrated_matrix[ad_markers_complement, col.order] # make sure rows and columns are in order using gene list and col.order
integrated_matrix = t(apply(integrated_matrix, 1, function(x) rescale(x, to = c(-2, 2)))) # normalize matrix between -2 and 2
cellType_gentreat = unique(sn_sub@meta.data[,c("cellType","gen_treat","celltype_gentreat")]) # get cellType-gen_treat table
cellType_gentreat= cellType_gentreat[match(colnames(integrated_matrix), cellType_gentreat$celltype_gentreat),] # put table in same order as matrix

cellType_gentreat= cellType_gentreat[match(col.order, cellType_gentreat$celltype_gentreat),] # put table in same order as matrix

# make top color bar annotation                                                                   
column_cellType = HeatmapAnnotation( 
  cellType = cellType_gentreat$cellType,
  gen_treat = cellType_gentreat$gen_treat,
  col = list(cellType = c("Neurons" = "#49796b", "Astrocytes" = "#409cb3", "Oligodendrocytes" = "#c6588f","OPCs" = "#a92a3f", "Endothelial" = "#bc8400","Microglia" = "#e4994b"),
             gen_treat = c("WT.Vehicle"="#C77CFF", "WT.PMX205"="#00BFC4", "Arc.Vehicle"="#7CAE00", "Arc.PMX205"="#F8766D")))

ad_SNonly_integrated_h <-Heatmap(integrated_matrix, name = "Expression", show_row_dend = F, show_column_dend = F, cluster_columns = F, cluster_rows = F, col = viridis(100), show_column_names = F, top_annotation = c(column_cellType))
ad_SNonly_integrated_h

##################################
########## FIGURE 3 ##############
##################################
mic <- subset(ad.combined.annotate, subset = cellType == "Microglia")
mic <- PercentageFeatureSet(mic, pattern = "^mt-", col.name = "percent.mt", assay = 'SCT')
mic <- SCTransform(mic, verbose = FALSE, assay="RNA", vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA"), return.only.var.genes = F,min_cells = 3)
mic <- RunPCA(mic, npcs = 30, verbose = FALSE, assay = "SCT")
mic <- RunUMAP(mic, reduction = "pca", dims = 1:30, assay = "SCT")
mic <- FindNeighbors(mic, reduction = "pca", dims = 1:30, assay = "SCT")
mic <- FindClusters(mic, resolution = 0.4, algorithm = 4, assay = "SCT") #, algorithm = 4
saveRDS(mic, file ="Splitseq_mic_050522.rds")

mic.markers <- FindAllMarkers(mic, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mic.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(mic,features = top10$gene,draw.lines = T,lines.width = 1 
          ,group.by = 'ident') + #NoLegend() + group.colors = new_palette_micro
  scale_fill_viridis( na.value = "white") +
  theme(text = element_text(size = 20, color = 'black'))

F1 <- FeaturePlot(mic, features = c("Cx3cr1"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
F2 <- FeaturePlot(mic, features = c("Csf1r"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
F3 <- FeaturePlot(mic, features = c("Itgam"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
F4 <- FeaturePlot(mic, features = c("P2ry12"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
F5 <- FeaturePlot(mic, features = c("Apoe"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
F6 <- FeaturePlot(mic, features = c("Trem2"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
F7 <- FeaturePlot(mic, features = c("Clec7a"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
F8 <- FeaturePlot(mic, features = c("Itgax"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
F9 <- FeaturePlot(mic, features = c("Cst7"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
F10 <- FeaturePlot(mic, features = c("Spp1"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
F11 <- FeaturePlot(mic, features = c("Rbfox1"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
F12 <- FeaturePlot(mic, features = c("Meg3"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
F13 <- FeaturePlot(mic, features = c("Cdk8"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
F14 <- FeaturePlot(mic, features = c("Lars2"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)

#highlight cluster 6, 10, 11, 12 back in the big UMAP
mic_6 <- WhichCells(object = mic, ident = 6)
mic_10 <- WhichCells(object = mic, ident = 10)
mic_11 <- WhichCells(object = mic, ident = 11)
mic_12 <- WhichCells(object = mic, ident = 12)
cluster6_highlight <- DimPlot(object = ad.combined.annotate, cells.highlight = mic_6, order = TRUE, group.by="cellType", label = T, label.color = "black", repel =T ) + scale_color_manual(labels = c("unselected", "mic_Cluster6"), values = c("#adacac","#572878")) 
cluster10_highlight <- DimPlot(object = ad.combined.annotate, cells.highlight = mic_10, order = TRUE, group.by="cellType", label = T, label.color = "black", repel =T ) + scale_color_manual(labels = c("unselected", "mic_Cluster10"), values = c("#adacac","#B58CD9")) 

#remove cluster 6 and 10
mic_sub <-subset(mic, idents = c("1","2","3","4","5","7","8","9","11", "12")) 
mic_sub <- RunPCA(mic_sub, npcs = 30, verbose = FALSE, assay = "SCT")
mic_sub <- RunUMAP(mic_sub, reduction = "pca", dims = 1:10, assay = "SCT")
mic_sub <- FindNeighbors(mic_sub, reduction = "pca", dims = 1:10, assay = "SCT")
mic_sub <- FindClusters(mic_sub, resolution = 0.5, algorithm = 4, assay = "SCT") #, algorithm = 4
ncol(x=mic_sub) #4537
DimPlot(mic_sub, reduction = "umap", group.by = "seurat_clusters", label = TRUE,repel = TRUE, label.size = 4)
saveRDS(mic_sub, file ="Splitseq_micsub_062622.rds")

# Microglia clusters Single-cell only complement genes
mic_sub <- readRDS("Splitseq_micsub_062622.rds")
DefaultAssay(mic_sub) = "SCT" 
mic_sub = NormalizeData(mic_sub)
mic_sub = FindVariableFeatures(mic_sub)
mic_sub = ScaleData(mic_sub)

DefaultAssay(mic_sub) = "integrated" 
mic_sub = NormalizeData(mic_sub)
mic_sub = FindVariableFeatures(mic_sub)
mic_sub = ScaleData(mic_sub)

DefaultAssay(mic_sub) = "RNA" 
mic_sub = NormalizeData(mic_sub)
mic_sub = FindVariableFeatures(mic_sub)
mic_sub = ScaleData(mic_sub)
ncol(x=mic_sub)
saveRDS(mic_sub, file ="Splitseq_micsub_042423.rds")

micsub.markers <- FindAllMarkers(mic_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(micsub.markers, 'Markers_integration_microglia_052323.txt',sep="\t",quote = FALSE, row.names = T)
markers = read.delim("Markers_integration_microglia_052323.txt")

####Figure 3A
mic_sub <- readRDS(file = "Splitseq_micsub_042423.rds")
as.matrix(table(mic_sub@meta.data$exp_gen_treat_sex))
ncol(x=mic_sub) #4537
ms_clusters_color <- colorRampPalette(c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99"))(n=11)
ms_dimplot <- DimPlot(mic_sub, reduction = "umap", group.by = "seurat_clusters",label = TRUE,repel = TRUE, label.size = 6, cols=ms_clusters_color)

####Figure 3B
micsub_barplot_cluster <- ggplot(mic_sub@meta.data, aes(fill=seurat_clusters, x=gen_treat, width=.5)) + geom_bar(position = "fill")+
  scale_fill_manual(values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99")) +
  theme(text = element_text(size = 12, color = 'black')) +  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title = element_text(size = 16))+ 
  theme(legend.position="right", legend.title =element_text(size=16),legend.text = element_text(size=14)) +
  theme(axis.text = element_text(size = 18), axis.text.x = element_text(size = 14, color = 'black'), axis.text.y = element_text(size = 14, color = 'black')) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + labs(x="Genotype.Treatment", y="Proportion") +
  ggtitle("Microglia cluster porportion in each genotype.treatment")

mic_sub_clusters_total <- as.data.frame(table(mic_sub@meta.data$gen_treat))

####Figure 3C
seuratcluster_levels <- c(2,6,1,9,5,7,3,8,4,10)
micsub_celltype <- mic_sub
micsub_celltype@meta.data$seurat_clusters <- factor(micsub_celltype@meta.data$seurat_clusters, # Change ordering manually
                                                    levels = c("2","6","1","9","5","7","3","8","4","10"))

micsub_barplot_exp <- ggplot(micsub_celltype@meta.data, aes(fill=experiment, x=seurat_clusters, width=.5)) + geom_bar(position = "fill")+ 
  scale_fill_manual(values = c("#000080","#F19CBB"),labels = c("SC","SN")) + 
  theme_bw() + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 16))+ 
  theme(legend.position="bottom", legend.title =element_text(size=16),legend.text = element_text(size=14))+
  theme(axis.text = element_text(size = 18), axis.text.x = element_text(size = 14, color = 'black'), axis.text.y = element_text(size = 14, color = 'black')) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +
  ggtitle("Microglia - SN+SC Experiment Porportion")

####Figure 3D
micsub_barplot_gentreat <- ggplot(micsub_celltype@meta.data, aes(fill=gen_treat, x=seurat_clusters, width=.5)) + geom_bar(position = "fill")+ 
  theme_bw() + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 16))+ 
  theme(legend.position="bottom", legend.title =element_text(size=16),legend.text = element_text(size=14))+
  theme(axis.text = element_text(size = 18), axis.text.x = element_text(size = 14, color = 'black'), axis.text.y = element_text(size = 14, color = 'black')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),panel.background = element_blank()) +
  ggtitle("Microglia - SN+SC Treatment Porportion")

####Figure 3E
microglia_markers <-c("Spi1",'Mef2c','Tmem119','Selplg', 'Csf1r', 'Tgfbr1','Tgfbr2', 'Hexb','Cx3cr1','P2ry12'
                      , 'Clec7a','Tyrobp', 'Apoe', 'Trem2', 'Fth1','B2m', 'C1qa','C1qb','C1qc','Lyz2' #DAM 1
                      , 'Lpl', 'Cst7','Fgf13','Spp1', 'Csf1', 'Axl', 'Itgax', 'Cd9', #DAM 2
                      'C5ar1',
                      "Cd44","Cdk8","Gphn",#Cluster6
                      "Inpp5d","Mertk", #Cluster9
                      "C5a","Hc")

DefaultAssay(mic_sub) = "integrated"
matrix = AverageExpression(mic_sub, group.by = "seurat_clusters", return.seurat = F)
integrated_matrix = matrix$integrated # get integrated average expression, other options are raw RNA and normalized SCT
microglia_markers = microglia_markers[microglia_markers %in% rownames(integrated_matrix)] # subset by genes of interest
col.order = as.character(c("2","6","1","9","5","7","3","8","4","10"))
integrated_matrix = integrated_matrix[microglia_markers,col.order] # make sure rows and columns are in order using gene list and col.order
integrated_matrix = t(apply(integrated_matrix, 1, function(x) rescale(x, to = c(-2, 2)))) # normalize matrix between -2 and 2
cluster_celltype = unique(mic_sub@meta.data[,c("seurat_clusters","cellType")]) # get cluster-celltype table
cluster_celltype = cluster_celltype[match(col.order,cluster_celltype$seurat_clusters),] # put table in same order as matrix
cluster_celltype$cellType = factor(cluster_celltype$cellType, levels = "Microglia")

# make top color bar annotation                                                                   
column_ha = HeatmapAnnotation( 
  seurat_clusters = cluster_celltype$seurat_clusters,
  col = list(seurat_clusters = c("1"="#a6cee3","2"="#1f78b4","3"="#b2df8a","4"="#33a02c","5"="#fb9a99","6"="#e31a1c","7"="#fdbf6f","8"="#ff7f00","9"="#cab2d6","10"="#6a3d9a","11"="#ffff99")))
# complex heatmap
micsub_clusters_all_heatmap <-Heatmap(integrated_matrix, name = "Expression", show_row_dend = F, show_column_dend = T, cluster_columns = F, cluster_rows = F, col = viridis(100), column_names_side = "top", top_annotation = column_ha)
micsub_clusters_all_heatmap

###Figure 3G
# average expression across clusters
mic_sub$Gen_treat = paste0(mic_sub$genotype,".",mic_sub$treatment)
mic_sub<- subset(mic_sub, subset = experiment == "mic")
ncol(x=mic_sub) #3269

# average expression across clusters
DefaultAssay(mic_sub) = "integrated"
matrix = AverageExpression(mic_sub, group.by = "Gen_treat", return.seurat = F)
integrated_matrix = matrix$integrated # get integrated average expression, other options are raw RNA and normalized SCT
microglia_markers = microglia_markers[microglia_markers %in% rownames(integrated_matrix)] # subset by genes of interest
col.order = as.character(c("WT.Vehicle","WT.PMX205","Arc.Vehicle","Arc.PMX205"))
integrated_matrix = integrated_matrix[microglia_markers,col.order] # make sure rows and columns are in order using gene list and col.order
integrated_matrix = t(apply(integrated_matrix, 1, function(x) rescale(x, to = c(-2, 2)))) # normalize matrix between -2 and 2
genotype_treatment = unique(mic_sub@meta.data[,c("genotype","treatment")]) # get cluster-celltype table
genotype_treatment = genotype_treatment[match(col.order,genotype_treatment$gen_treat),] # put table in same order as matrix

# make top color bar annotation                                                                   
column_ha = HeatmapAnnotation( 
  genotype = genotype_treatment$genotype,
  treatment = genotype_treatment$treatment,
  col = list(genotype = c("WT"="#F68282","Arc"="#1f78b4"),
             treatment = c("Vehicle"="#b2df8a","PMX205"="#33a02c")))

# complex heatmap
micsub_gentreat_heatmap <-Heatmap(integrated_matrix, name = "Expression",
                                  show_row_dend = F, 
                                  show_column_dend = F, 
                                  cluster_columns = F,
                                  cluster_rows = F, 
                                  col = viridis(100),
                                  column_names_side = "top",
                                  top_annotation = column_ha)
micsub_gentreat_heatmap


##################################
########## FIGURE 4 ##############
##################################
#Recluster astrocytes
ad.combined.annotate <- readRDS("sc_sn_postSCT_20kanchors_RPCA_res1_annotated_050522.rds")
ast <- subset(ad.combined.annotate, subset = cellType == "Astrocytes")
ast <- PercentageFeatureSet(ast, pattern = "^mt-", col.name = "percent.mt", assay = 'SCT')
ast <- SCTransform(ast, verbose = FALSE, assay="RNA", vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA"), return.only.var.genes = F,min_cells = 3)
ast <- RunPCA(ast, npcs = 30, verbose = FALSE, assay = "SCT")
ast <- RunUMAP(ast, reduction = "pca", dims = 1:30, assay = "SCT")
ast <- FindNeighbors(ast, reduction = "pca", dims = 1:30, assay = "SCT")
ast <- FindClusters(ast, resolution = 0.3, algorithm = 4, assay = "SCT") #, algorithm = 4
head(Idents(ast),5)
DefaultAssay(ast) <- 'SCT'
clusters <- FetchData(object =ast, vars = c('seurat_clusters'))
head(x = FetchData(object = ast, vars = c('group', 'ident')))
unique(ast@meta.data$experiment)
saveRDS(ast, file ="Splitseq_ast_062622.rds")

ast.markers <- FindAllMarkers(ast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ast.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
write.table(top10, "astrocytes_clusters_top10gene_070122.txt",sep="\t",row.names = F)

aF1 <- FeaturePlot(ast, features = c("Clu"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
aF2 <- FeaturePlot(ast, features = c("Gfap"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
aF4 <- FeaturePlot(ast, features = c("Slc1a3"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)
aF5 <- FeaturePlot(ast, features = c("Apoe"),cols = c("#adacac", "#b30000"), pt.size =0.5, order = T)

ncol(x=ast) #4507

ast <- readRDS("Splitseq_ast_062622.rds")
DefaultAssay(ast) = "SCT" 
ast = NormalizeData(ast)
ast = FindVariableFeatures(ast)
ast = ScaleData(ast)

DefaultAssay(ast) = "integrated" 
ast = NormalizeData(ast)
ast = FindVariableFeatures(ast)
ast = ScaleData(ast)

DefaultAssay(ast) = "RNA" 
ast = NormalizeData(ast)
ast = FindVariableFeatures(ast)
ast = ScaleData(ast)
ncol(x=ast) #4507
saveRDS(ast, file ="Splitseq_ast_042423.rds")

gentreat_levels <- c("WT.Vehicle", "WT.PMX205", "Arc.Vehicle", "Arc.PMX205")

####Figure4A 
ast <- readRDS("Splitseq_ast_042423.rds")
ast_clusters_color <- colorRampPalette(c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928"))(n=12)
ast_dimplot <- DimPlot(ast, reduction = "umap", group.by = "seurat_clusters",label = TRUE,repel = TRUE, label.size = 8, cols =ast_clusters_color)
ast_dimplot

####Figure4B
ast_barplot_clusters <- ggplot(ast@meta.data, aes(x=gen_treat, fill = seurat_clusters)) + geom_bar(position = "fill")+ scale_fill_manual(values = ast_clusters_color) + 
  theme(text = element_text(size = 12, color = 'black')) +   theme_bw() + theme(plot.title = element_text(size = 18), axis.title = element_text(size = 16))+ 
  theme(legend.position="right", legend.title =element_text(size=16),legend.text = element_text(size=14))+
  theme(axis.text = element_text(size = 18), axis.text.x = element_text(size = 14, color = 'black'), axis.text.y = element_text(size = 14, color = 'black')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +
  labs(x="Genotype.Treatment", y="Proportion")+
  ggtitle("Cluster porportion in each genotype.treatment")

ast_barplot_clusters$data$gen_treat <- factor(x = ast_barplot_clusters$data$gen_treat, levels = gentreat_levels) # change the order of the factor levels
ast_barplot_clusters

ast_clusters_total <- as.data.frame(table(ast@meta.data$gen_treat))
#Arc.PMX205 1050
#Arc.Vehicle 990
#WT.PMX205 1207
#WT.Vehicle 1260

####Figure 4C
ast_barplot_genetreat <- ggplot(ast_celltype@meta.data, aes(fill=gen_treat, x=seurat_clusters, width=.5)) + geom_bar(position = "fill")+ 
  theme_bw() + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 16))+ 
  theme(legend.position="right", legend.title =element_text(size=16),legend.text = element_text(size=14))+
  theme(axis.text = element_text(size = 18), axis.text.x = element_text(size = 14, color = 'black'), axis.text.y = element_text(size = 14, color = 'black')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),panel.background = element_blank()) +
  ggtitle("Astrocyte - SN+SC Treatment Porportion")

ast_barplot_genetreat

###Figure 4D
# genes of interest
Ast_genes <-c('Slc1a2','Slc1a3','Aqp4','Apoe','Clu','Aldh1l1','Vim','Glu1','S100a6', # homeostatic
              'Lcn2', 'Cxcl10','Cp','Steap4', 'Timp1','Cd44',  'Aspg', 'Gfap', 'Serpina3n',
              'C3', 'H2-T23', 'H2-D1','Gbp2',  'Srgn',
              'S100a10',  'Cd109', 'Slc10a6', 'Tm4sf1',
              "App") 

# average expression across clusters
DefaultAssay(ast) = "SCT"
matrix = AverageExpression(ast, group.by = "seurat_clusters", return.seurat = F)
SCT_matrix = matrix$SCT
Ast_genes = Ast_genes[Ast_genes %in% rownames(SCT_matrix)]
col.order = as.character(c("9","12","3","11","10","4","5","6","8","2","7","1"))
SCT_matrix = SCT_matrix[Ast_genes,col.order] # make sure rows and columns are in order using gene list and col.order
SCT_matrix = t(apply(SCT_matrix, 1, function(x) rescale(x, to = c(-2, 2)))) # normalize matrix between -2 and 2

cluster_celltype = unique(ast@meta.data[,c("seurat_clusters","cellType")]) # get cluster-celltype table
cluster_celltype = cluster_celltype[match(col.order,cluster_celltype$seurat_clusters),] # put table in same order as matrix
cluster_celltype$cellType = factor(cluster_celltype$cellType, levels = "Astrocytes")

# make top color bar annotation                                                                   
column_ha = HeatmapAnnotation( 
  seurat_clusters = cluster_celltype$seurat_clusters,
  col = list(seurat_clusters = c("1"="#a6cee3","2"="#1f78b4","3"="#b2df8a","4"="#33a02c","5"="#fb9a99","6"="#e31a1c","7"="#fdbf6f","8"="#ff7f00","9"="#cab2d6","10"="#6a3d9a","11"="#ffff99","12"="#b15928")))

# complex heatmap
ast_clusters_all_heatmap <-Heatmap(SCT_matrix, name = "Expression", show_row_dend = F, show_column_dend = T, cluster_columns = F, cluster_rows = F, col = viridis(100), column_names_side = "top", top_annotation = column_ha)
ast_clusters_all_heatmap

####Figure 4E 
# average expression across clusters
DefaultAssay(ast) = "SCT"
ast$Gen_treat = paste0(ast$genotype,".",ast$treatment)
matrix = AverageExpression(ast, group.by = "Gen_treat", return.seurat = F)
SCT_matrix = matrix$SCT # get integrated average expression, other options are raw RNA and normalized SCT
Ast_genes = Ast_genes[Ast_genes %in% rownames(SCT_matrix)] # subset by genes of interest
col.order = as.character(c("WT.Vehicle","WT.PMX205","Arc.Vehicle","Arc.PMX205"))
SCT_matrix = SCT_matrix[Ast_genes,col.order] # make sure rows and columns are in order using gene list and col.order
SCT_matrix = SCT_matrix[ast_c7_genes,col.order] # make sure rows and columns are in order using gene list and col.order
SCT_matrix = t(apply(SCT_matrix, 1, function(x) rescale(x, to = c(-2, 2)))) # normalize matrix between -2 and 2
genotype_treatment = unique(ast@meta.data[,c("genotype","treatment")]) # get cluster-celltype table
genotype_treatment = genotype_treatment[match(col.order,genotype_treatment$gen_treat),] # put table in same order as matrix

# make top color bar annotation                                                                   
column_ha = HeatmapAnnotation( 
  genotype = genotype_treatment$genotype,
  treatment = genotype_treatment$treatment,
  col = list(genotype = c("WT"="#F68282","Arc"="#1f78b4"),
             treatment = c("Vehicle"="#b2df8a","PMX205"="#33a02c")))

# complex heatmap
ast_gentreat_heatmap <-Heatmap(SCT_matrix, name = "Expression",
                               show_row_dend = F, 
                               show_column_dend = F, 
                               cluster_columns = F,
                               cluster_rows = F, 
                               col = viridis(100),
                               column_names_side = "top",
                               top_annotation = column_ha)
ast_gentreat_heatmap


##################################
########## FIGURE 5 ##############
##################################
ad.combined.1 <- readRDS("sc_sn_postSCT_20kanchors_RPCA_res1_annotated_050522.rds")
ad.combined.1 <- subset(ad.combined.1, subset = cellType == "Mixed", invert = T)

##Set up Arc.PMX.cellchat
Arc.PMX <- subset(ad.combined.1, subset = gen_treat == "Arc.PMX205" )
ncol(x=Arc.PMX) #12158
Arc.PMX[["new.ident"]] <- Idents(object = Arc.PMX)
Arc.PMX <- SetIdent(Arc.PMX, value = Arc.PMX@meta.data$new.ident)
Idents(Arc.PMX)
Arc.PMX.data.input <- GetAssayData(Arc.PMX, assay = "SCT", slot = "data") # normalized data matrix
Arc.PMX.labels <- Arc.PMX@active.ident
unique(Arc.PMX@active.ident)
Arc.PMX.meta <- Arc.PMX@meta.data
Arc.PMX.matrix <- as.matrix(table(Arc.PMX@active.ident))

Arc.PMX.cellchat <- createCellChat(object = Arc.PMX.data.input, meta = Arc.PMX.meta, group.by = "cellType")
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
#CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-receptor", "Cell-cell Contact")) # use Secreted Signaling
Arc.PMX.cellchat@DB <- CellChatDB.use

#Subset the expression data of signaling genes for saving computation cost
Arc.PMX.cellchat <- subsetData(Arc.PMX.cellchat) # This step is necessary even if using the whole database
##Issue identified!! Please check the official Gene Symbol of the following genes:  H2-BI H2-Ea-ps
future::plan("multiprocess", workers = 4) # do parallel

Arc.PMX.cellchat <- identifyOverExpressedGenes(Arc.PMX.cellchat)
Arc.PMX.cellchat <- identifyOverExpressedInteractions(Arc.PMX.cellchat)

#drop unused levels using 'droplevels' function
unique(Arc.PMX.cellchat@idents)
Arc.PMX.cellchat@idents = droplevels(Arc.PMX.cellchat@idents, exclude = setdiff(levels(Arc.PMX.cellchat@idents),unique(Arc.PMX.cellchat@idents)))
unique(Arc.PMX.cellchat@idents)

#computing average expression per cell group with "truncatedMean"
Arc.PMX.cellchat <- computeCommunProb(Arc.PMX.cellchat, population.size = F, type = "truncatedMean", trim = 0.2) #, type = "truncatedMean", trim = 0.2

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
#the minmum number of cells required in each cell group for cell-cell communication
Arc.PMX.cellchat <- filterCommunication(Arc.PMX.cellchat, min.cells = 10)

Arc.PMX.cellchat <- computeCommunProbPathway(Arc.PMX.cellchat)
Arc.PMX.cellchat <- aggregateNet(Arc.PMX.cellchat)

saveRDS(Arc.PMX.cellchat, file ="ArcPMX_cellchat_trim02_mincells10_nosubsetDB_100422.rds")


##Setup for Arc.V
Arc.V <- subset(ad.combined.1, subset = gen_treat == "Arc.Vehicle" )
ncol(x=Arc.V)#12521
Arc.V[["new.ident"]] <- Idents(object = Arc.V)
Arc.V <- SetIdent(Arc.V, value = Arc.V@meta.data$new.ident)
Idents(Arc.V)
Arc.V.data.input <- GetAssayData(Arc.V, assay = "SCT", slot = "data") # normalized data matrix
Arc.V.labels <- Arc.V@active.ident
Arc.V.labels
unique(Arc.V@active.ident)
Arc.V.meta <- Arc.V@meta.data
Arc.V.matrix <- as.matrix(table(Arc.V@active.ident))

Arc.V.cellchat <- createCellChat(object = Arc.V.data.input, meta = Arc.V.meta, group.by = "cellType")
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
#CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-receptor", "Cell-cell Contact")) # use Secreted Signaling
Arc.V.cellchat@DB <- CellChatDB.use
Arc.V.cellchat <- subsetData(Arc.V.cellchat) # This step is necessary even if using the whole database
#Issue identified!! Please check the official Gene Symbol of the following genes: H2-BI H2-Ea-ps 
future::plan("multiprocess", workers = 4) # do parallel

Arc.V.cellchat <- identifyOverExpressedGenes(Arc.V.cellchat)
Arc.V.cellchat <- identifyOverExpressedInteractions(Arc.V.cellchat)

#drop unused levels using 'droplevels' function
unique(Arc.V.cellchat@idents)
Arc.V.cellchat@idents = droplevels(Arc.V.cellchat@idents, exclude = setdiff(levels(Arc.V.cellchat@idents),unique(Arc.V.cellchat@idents)))
unique(Arc.V.cellchat@idents)
Arc.V.cellchat <- computeCommunProb(Arc.V.cellchat,  population.size = F, type = "truncatedMean", trim = 0.2) #, type = "truncatedMean", trim = 0.2

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
Arc.V.cellchat <- filterCommunication(Arc.V.cellchat, min.cells = 10)

Arc.V.cellchat <- computeCommunProbPathway(Arc.V.cellchat)
Arc.V.cellchat <- aggregateNet(Arc.V.cellchat)

saveRDS(Arc.V.cellchat, file ="ArcV_cellchat_trim02_mincells10_nosubsetDB_100422.rds")


##Setup for WT.PMX205
WT.PMX <- subset(ad.combined.1, subset = gen_treat == "WT.PMX205" )
ncol(x=WT.PMX) #16321
WT.PMX[["new.ident"]] <- Idents(object = WT.PMX)
WT.PMX <- SetIdent(WT.PMX, value = WT.PMX@meta.data$new.ident)
Idents(WT.PMX)
WT.PMX.data.input <- GetAssayData(WT.PMX, assay = "SCT", slot = "data") # normalized data matrix
WT.PMX.labels <- WT.PMX@active.ident
unique(WT.PMX@active.ident)
WT.PMX.meta <- WT.PMX@meta.data
WT.PMX.matrix <- as.matrix(table(WT.PMX@active.ident))

WT.PMX.cellchat <- createCellChat(object = WT.PMX.data.input, meta = WT.PMX.meta, group.by = "cellType")
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
#CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-receptor", "Cell-cell Contact")) # use Secreted Signaling
WT.PMX.cellchat@DB <- CellChatDB.use
WT.PMX.cellchat <- subsetData(WT.PMX.cellchat) # This step is necessary even if using the whole database
#Issue identified!! Please check the official Gene Symbol of the following genes: H2-BI H2-Ea-ps 
future::plan("multiprocess", workers = 4) # do parallel

WT.PMX.cellchat <- identifyOverExpressedGenes(WT.PMX.cellchat)
WT.PMX.cellchat <- identifyOverExpressedInteractions(WT.PMX.cellchat)

#drop unused levels using 'droplevels' function
unique(WT.PMX.cellchat@idents)
WT.PMX.cellchat@idents = droplevels(WT.PMX.cellchat@idents, exclude = setdiff(levels(WT.PMX.cellchat@idents),unique(WT.PMX.cellchat@idents)))
unique(WT.PMX.cellchat@idents)
WT.PMX.cellchat <- computeCommunProb(WT.PMX.cellchat,  population.size = F, type = "truncatedMean", trim = 0.2) #, type = "truncatedMean", trim = 0.2

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
WT.PMX.cellchat <- filterCommunication(WT.PMX.cellchat, min.cells = 10)

WT.PMX.cellchat <- computeCommunProbPathway(WT.PMX.cellchat)
WT.PMX.cellchat <- aggregateNet(WT.PMX.cellchat)

saveRDS(WT.PMX.cellchat, file ="WTPMX_cellchat_trim02_mincells10_nosubsetDB_100322.rds")


##Setup for WT.V
WT.V <- subset(ad.combined.1, subset = gen_treat == "WT.Vehicle" ) 
ncol(x=WT.V)#12140
WT.V[["new.ident"]] <- Idents(object = WT.V)
WT.V <- SetIdent(WT.V, value = WT.V@meta.data$new.ident)
Idents(WT.V)
WT.V.data.input <- GetAssayData(WT.V, assay = "SCT", slot = "data") # normalized data matrix
WT.V.labels <- WT.V@active.ident
WT.V.labels
unique(WT.V@active.ident)
WT.V.meta <- WT.V@meta.data
WT.V.matrix <- as.matrix(table(WT.V@active.ident))

WT.V.cellchat <- createCellChat(object = WT.V.data.input, meta = WT.V.meta, group.by = "cellType")
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
#CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-receptor", "Cell-cell Contact")) # use Secreted Signaling
WT.V.cellchat@DB <- CellChatDB.use
WT.V.cellchat <- subsetData(WT.V.cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel

WT.V.cellchat <- identifyOverExpressedGenes(WT.V.cellchat)
WT.V.cellchat <- identifyOverExpressedInteractions(WT.V.cellchat)

#drop unused levels using 'droplevels' function
unique(WT.V.cellchat@idents)
WT.V.cellchat@idents = droplevels(WT.V.cellchat@idents, exclude = setdiff(levels(WT.V.cellchat@idents),unique(WT.V.cellchat@idents)))
unique(WT.V.cellchat@idents)
WT.V.cellchat <- computeCommunProb(WT.V.cellchat,  population.size = F, type = "truncatedMean", trim = 0.2) #, type = "truncatedMean", trim = 0.2

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
WT.V.cellchat <- filterCommunication(WT.V.cellchat, min.cells = 10)

WT.V.cellchat <- computeCommunProbPathway(WT.V.cellchat)
WT.V.cellchat <- aggregateNet(WT.V.cellchat)

saveRDS(WT.V.cellchat, file ="WTV_cellchat_trim02_mincells10_nosubsetDB_100322.rds")

####Figure 5A
Arc.PMX.cellchat <- readRDS("ArcPMX_cellchat_trim02_mincells10_nosubsetDB_100422.rds")
Arc.V.cellchat <- readRDS("ArcV_cellchat_trim02_mincells10_nosubsetDB_100422.rds")
WT.PMX.cellchat <- readRDS("WTPMX_cellchat_trim02_mincells10_nosubsetDB_100322.rds")
WT.V.cellchat <- readRDS("WTV_cellchat_trim02_mincells10_nosubsetDB_100322.rds")

object.list <- list(ArcVeh = Arc.V.cellchat, WTVeh = WT.V.cellchat)
Veh.cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
Veh.cellchat
#An object of class CellChat created from a merged object with multiple datasets 
#nosubsetDB
#927 signaling genes.
# 24661 cells.
Veh_r1 <- rankNet(Veh.cellchat, mode = "comparison", measure = c("weight", "count"),comparison = c(1, 2), stacked = T, do.stat = TRUE, color.use = c(ArcVeh="#f8766d",WTVeh="#619CFF"))
Veh_r1

#new WT VS Arc in specific order that will match Arc vs ArcPMX
Veh_r1_1 <- rankNet(Veh.cellchat, mode = "comparison", measure = c("weight", "count"),comparison = c(1, 2), stacked = T, do.stat = TRUE,
                    signaling=c("TENASCIN","ANGPT","PROS","AGRN",
                                "FGF","MPZ","CDH5","CLDN","CDH",
                                "VTN","COLLAGEN","LAMININ","MAG", #significant in ArcVeh
                                "CSF","PECAM1",
                                "FN1","SEMA4","EGF","PTN",  #significant in ArcVeh
                                "NRG"),#significant in WTVeh
                    title = "ArcVeh and WTVeh information flow of significant signaling pathway",
                    font.size = 10, color.use = c(ArcVeh="#f8766d",WTVeh="#619CFF"))
Veh_r1_1

####Figure 5B
object.list <- list(ArcVeh = Arc.V.cellchat, ArcPMX = Arc.PMX.cellchat)
Arc.cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
Arc.cellchat
#after put the threshold to trim=0.2
#657 signaling genes.
#24679 cells.
Arc_r1 <- rankNet(Arc.cellchat, mode = "comparison", measure = c("weight", "count"),comparison = c(1, 2), stacked = T, do.stat = TRUE)
Arc_r1

#new ArcPMX VS Arc in specific order
Arc_r1_1 <- rankNet(Arc.cellchat, mode = "comparison", measure = c("weight", "count"),comparison = c(1, 2), stacked = T, do.stat = TRUE,
                  signaling=c("CLDN","CDH","SEMA5","FGF","MPZ","CDH5","SEMA7","TGFb","BMP", #significant in ArcPMX
                              "PECAM1","TENASCIN","FN1","EPHA","COLLAGEN","PROS","AGRN","LAMININ","ANGPT","CSF","MAG","VTN"),#significant in ArcVeh
                  title = "ArcPMX and ArcVeh information flow of significant signaling pathway",
                  font.size = 15)
Arc_r1_1

####Figure 5C
#Arc.cellchat overall signaling pattern
#WTveh vs ArcVeh vs ArcPMX significant pathways from Arcpmx vs Arcveh Relative information flow
paletteCells_labeled <- c("Neurons" = "#49796b", "Astrocytes" = "#409cb3", "Oligodendrocytes" = "#c6588f",
                          "OPCs" = "#e4994b", "Endothelial" = "#a92a3f", "Pericytes" = "#bc8400",
                          "Microglia" = "#a99fbf")

Arc.PMX.cellchat_a1 <- netAnalysis_computeCentrality(Arc.PMX.cellchat, slot.name = "netP")
Arc.V.cellchat_a1 <- netAnalysis_computeCentrality(Arc.V.cellchat, slot.name = "netP")
object.list <- list(ArcVeh = Arc.V.cellchat_a1, ArcPMX = Arc.PMX.cellchat_a1)
Arc.cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# combining all the identified signaling pathways from different datasets 

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", color.use=paletteCells_labeled, signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", color.use=paletteCells_labeled, signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)

Arc_ht1_out = netAnalysis_signalingRole_heatmap(object.list[[i]], color.use=paletteCells_labeled, pattern = "all"
                                                ,font.size = 12, signaling = c("NRXN","NCAM","NRG","NEGR","CADM","CNTN","PTN","PTPRM","SEMA6","NGL",
                                                                               "SEMA3","LAMININ","MAG","CDH","VEGF","PDGF","VISTA","SEMA5","PECAM1","EPHA",
                                                                               "CSF","EGF","SEMA4","VTN","COLLAGEN","FN1","CLDN","JAM","TENASCIN","ANGPT",
                                                                               "FGF","MPZ","PROS","CDH5","AGRN",
                                                                               "SEMA7","TGFb","BMP"),
                                                title = names(object.list)[i], width = 8, height = 14, color.heatmap = "Greens")
Arc_ht2_out = netAnalysis_signalingRole_heatmap(object.list[[i+1]], color.use=paletteCells_labeled, pattern = "all"
                                                ,font.size = 12, signaling = c("NRXN","NCAM","NRG","NEGR","CADM","CNTN","PTN","PTPRM","SEMA6","NGL",
                                                                               "SEMA3","LAMININ","MAG","CDH","VEGF","PDGF","VISTA","SEMA5","PECAM1","EPHA",
                                                                               "CSF","EGF","SEMA4","VTN","COLLAGEN","FN1","CLDN","JAM","TENASCIN","ANGPT",
                                                                               "FGF","MPZ","PROS","CDH5","AGRN",
                                                                               "SEMA7","TGFb","BMP"), 
                                                title = names(object.list)[i+1], width = 8, height = 14, color.heatmap = "Greens")

Arc_ht1_out + Arc_ht2_out

#Veh.cellchat overall signaling pattern
Arc.V.cellchat_a1 <- netAnalysis_computeCentrality(Arc.V.cellchat, slot.name = "netP")
WT.V.cellchat_a1 <- netAnalysis_computeCentrality(WT.V.cellchat, slot.name = "netP")
object.list <- list(WTVeh = WT.V.cellchat_a1, ArcVeh = Arc.V.cellchat_a1)
Veh.cellchat <- mergeCellChat(object.list, add.names = names(object.list))

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
Veh_ht1_out = netAnalysis_signalingRole_heatmap(object.list[[i]], color.use=paletteCells_labeled, pattern = "all"
                                                ,font.size = 12, signaling = c("NRXN","NCAM","NRG","NEGR","CADM","CNTN","PTN","PTPRM","SEMA6","NGL",
                                                                               "SEMA3","LAMININ","MAG","CDH","VEGF","PDGF","VISTA","SEMA5","PECAM1","EPHA",
                                                                               "CSF","EGF","SEMA4","VTN","COLLAGEN","FN1","CLDN","JAM","TENASCIN","ANGPT",
                                                                               "FGF","MPZ","PROS","CDH5","AGRN",
                                                                               "SEMA7","TGFb","BMP"),
                                                title = names(object.list)[i], width = 8, height = 14, color.heatmap = "Greens")
Veh_ht2_out = netAnalysis_signalingRole_heatmap(object.list[[i+1]], color.use=paletteCells_labeled, pattern = "all"
                                                ,font.size = 12, signaling = c("NRXN","NCAM","NRG","NEGR","CADM","CNTN","PTN","PTPRM","SEMA6","NGL",
                                                                               "SEMA3","LAMININ","MAG","CDH","VEGF","PDGF","VISTA","SEMA5","PECAM1","EPHA",
                                                                               "CSF","EGF","SEMA4","VTN","COLLAGEN","FN1","CLDN","JAM","TENASCIN","ANGPT",
                                                                               "FGF","MPZ","PROS","CDH5","AGRN",
                                                                               "SEMA7","TGFb","BMP"), 
                                                title = names(object.list)[i+1], width = 8, height = 14, color.heatmap = "Greens")
Veh_ht1_out + Veh_ht2_out

