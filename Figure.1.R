library(data.table)
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)
library(patchwork)
library(pheatmap)
library(R.utils)
library(harmony)
library(clusterProfiler)
library(NMF)
library(ggalluvial)
##############################################################################################################################################
#Fig1
rm(list=ls())
#Loading all 14 HCC single-cell samples from our project, and 61 HCC patient samples from public dataset of  GSE149614, GSE151530, GSE156625, GSE189903, and GSE202642
#Integrated all 75 samples and create the single-cell atlas
data <- list(75 samples)
scobj <- merge(x=data[[1]], y = data[-1])
scobj <- JoinLayers(scobj)
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj <- Matrix(scobj@assays$RNA$counts,sparse = T)
scobj <- CreateSeuratObject(counts = scobj, assay = "RNA", min.cells = 3)
scobj <- subset(x = scobj, subset = percent.mt < 31 & nFeature_RNA >299 & nFeature_RNA <8001)

scobj <- NormalizeData(scobj, normalization.method = "LogNormalize", scale.factor = 10000)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 1500)
scobj <- ScaleData(scobj, features = rownames(scobj))
scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj),reduction.name = "pca")
ElbowPlot(scobj,60)
num=30
scobj <- RunHarmony(scobj,reduction = "pca" ,group.by.vars = 'orig.ident',reduction.save = "harmony")
scobj <- RunUMAP(scobj, reduction = "harmony", dims = 1:num, reduction.name = "umap", min.dist = 0.01)
scobj <- FindNeighbors(scobj, reduction = "harmony", dims = 1:num)
scobj <- FindClusters(scobj, resolution = 0.1)

#Create the subset of 14-patient samples of our project
Idents(scobj2) = scobj$Source
scobj2<-scobj[, Idents(scobj) %in% c("OurProject")]
#Save the total single-cell atlas of 75 patient samples and 14-patient samples of our project, respectively
saveRDS(scobj, file = "Total single-cell atlas.rds")
saveRDS(scobj2, file = "Our single-cell atlas.rds")
#####################################################
#UMAPplot
#Source represent the scRNA-seq data from our project and external public datasets
#MajorCellType includes B, endothelial, myeloid, HCC, T/NK cells and fibroblasts
DimPlot(scobj, reduction="umap", pt.size =0.2, repel = TRUE, raster = FALSE,
        group.by=c("MajorCellType"), alpha = 0.05, split.by = 'Source',
        cols = c('#8ecfc9','#ffa510','#f47254','#002c53','#8983bf','#67a583'))
#####################################################
#Dotplot
scobj$MajorCellType = factor(scobj$MajorCellType, levels =  c('B cells','T/NK','Myeloid', 'Fibroblasts','Endothelial cells','HCC'))
Markers <- c('MS4A1','MZB1','CD79A','CD3D','CD3E','CD3G','KLRB1','KLRD1','CD68','C1QA','C1QB','CSF1R',
                        'DCN','COL1A1','LUM','PECAM1','PLVAP','VWF','APOA1','APOB','APOC3','ALB')
DotPlot(scobj,features = Markers, group.by = "MajorCellType")
+RotatedAxis()+scale_color_gradientn(colours = c('#7b95c6','#D2D1D0','#FFCC33'))
#####################################################
#Featureplot
FeaturePlot(scobj,features = c('CD3D','CD3E','KLRB1','KLRD1'), reduction = "umap",cols = c("#d9dee7", "#FFB600"))
FeaturePlot(scobj,features = c('CD79A','MZB1','CD68','C1QA'), reduction = "umap",cols = c("#d9dee7", "#FFB600"))
FeaturePlot(scobj,features = c('ALB','APOA1','PECAM1','TAGLN'), reduction = "umap",cols = c("#d9dee7", "#FFB600"))
#####################################################
#CellChat analysis of  our 14-patietns data at MainCellType level
#MainCellType include B, endothelial, mast, dendritic, HCC, T, NK cells, neutrophils, monocytes, macrophages and fibroblasts
data.input <- scobj2[[scobj2@active.assay]]@layers$data
row.names(data.input) <- row.names(scobj2)
colnames(data.input) <- colnames(scobj2)
meta = scobj2@meta.data
unique(meta$MainCellType)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "MainCellType")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat,thresh.p = 0.05)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean",trim = 0.05,nboot = 1,raw.use = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#CirclePlot of CellChat
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat2@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#NMF analysis of CellChat
selectK(cellchat, pattern = "outgoing")
nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, height = 15, width = 6, font.size = 5)
netAnalysis_river(cellchat2, pattern = "outgoing")
netAnalysis_dot(cellchat2, pattern = "outgoing")
##############################################################################################################################################