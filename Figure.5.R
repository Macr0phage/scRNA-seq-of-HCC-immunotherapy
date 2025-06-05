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
library(rlang)
library(ggunchull)
library(jjAnno)
library(scRNAtoolVis)
library(Rcpp)
library(cellcall)
library(networkD3)
##############################################################################################################################################
#Fig5
rm(list=ls())
#Loading the 14-patient samples of our project
scobj <- readRDS ("Our single-cell atlas.rds")
#CellChat analysis of ligand-receptor interactions between HCC_FA and MacroLA
data.input <- scobj[[scobj@active.assay]]@layers$data
row.names(data.input) <- row.names(scobj)
colnames(data.input) <- colnames(scobj)
meta = scobj@meta.data
unique(meta$FineCellType)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "FineCellType")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat,thresh.p = 0.05)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", raw.use = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netVisual_bubble(cellchat, remove.isolate = FALSE,angle.x = 45,
                 sources.use = c('HCC_FA','MacroLA'),targets.use = c('HCC_FA','MacroLA'))

#Create subset of myeloid cells from 14-patient samples of our project
scobjMye = scobj[,scobj@meta.data$MajorCellType %in% c('Myeloid')]
scobjMye <- NormalizeData(scobjMye, normalization.method = "LogNormalize", scale.factor = 10000)
scobjMye <- FindVariableFeatures(scobjMye, selection.method = "vst", nfeatures = 1500)
scobjMye <- ScaleData(scobjMye, features = rownames(scobjMye))
scobjMye <- RunPCA(scobjMye, features = VariableFeatures(object = scobjMye),reduction.name = "pca")
ElbowPlot(scobjMye,60)
num=25
scobjMye <- RunHarmony(scobjMye,reduction = "pca" ,group.by.vars = 'orig.ident',reduction.save = "harmony")
scobjMye <- RunUMAP(scobjMye, reduction = "harmony", dims = 1:num, reduction.name = "umap", min.dist = 0.5)
scobjMye <- FindNeighbors(scobjMye, reduction = "harmony", dims = 1:num)
scobjMye <- FindClusters(scobjMye, resolution = 0.1)
DimPlot(scobjMye,reduction="umap", pt.size =0.2, repel = TRUE, raster = FALSE,
               group.by=c("seurat_clusters",'Myesubtype'), alpha = 0.1)
#The correlation of gene expression levels in MacroLA
FeaturePlot(scobjMye, features=c("C1QA", "FABP1"),blend = TRUE, repel = TRUE, order = TRUE, alpha = 0.2, 
            cols = c("gray80","red", "green"),pt.size = 0.1, raster = F) +theme(aspect.ratio = 1)
FeaturePlot(scobjMye, features=c("CD68", "APOA1"),blend = TRUE, repel = TRUE, order = TRUE, alpha = 0.2, 
            cols = c("gray80","red", "green"),pt.size = 0.1, raster = F) +theme(aspect.ratio = 1)

#CellCall analysis of transcription factor based on the receptor gene expression of MacroLA
test <- CreateObject_fromSeurat(Seurat.object=scobjMye, slot="counts",cell_type="Myesubtype",data_source="UMI",
                                scale.factor = 10^6,Org = "Homo sapiens") 
mt <- TransCommuProfile(object = test, pValueCor = 0.05, CorValue = 0.1, topTargetCor=1, p.adjust = 0.05,
                        use.type="mean",probs = 0.9,method="weighted",IS_core = TRUE,Org = 'Homo sapiens')
n <- mt@data$expr_l_r_log2_scale
egmt <- mt@data$gsea.list$"MacroLA"
egmt.df <- data.frame(egmt)
flag.index <- which(egmt.df$p.adjust < 0.05)
ridgeplot.DIY(x=egmt, fill="p.adjust", showCategory=flag.index, core_enrichment = T,
              orderBy = "NES", decreasing = FALSE)
p.tf <- names(mt@data$gsea.list$"MacroLA"@geneSets)
getGSEAplot(gsea.list=mt@data$gsea.list, geneSetID=c("PPARD"), 
            myCelltype="MacroLA", fc.list=mt@data$fc.list,  
            selectedGeneID =mt@data$gsea.list$"MacroLA"@geneSets$PPARD,
            mycol = NULL)
##############################################################################################################################################