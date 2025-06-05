library(data.table)
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)
library(patchwork)
library(SingleR)
library(celldex)
library(pheatmap)
library(cellcall)
library(networkD3)
library(R.utils)
library(harmony)
library(hdf5r)
library(SCopeLoomR)
library(SCENIC)
library(SingleCellExperiment)
######################################################################################################################################
#Fig7
rm(list=ls())
#Loading the 14-patient samples of our project
scobj <- readRDS ("Our single-cell atlas.rds")
#Create subset of responders and non-responders from 14-patient samples of our project
scobjR = scobj[,scobj@meta.data$Response %in% c('Responder')]
scobjNR = scobj[,scobj@meta.data$Response %in% c('Non-responder')]
#Analysis of Responders
data.inputR <- scobjR[[scobjR@active.assay]]@layers$data
row.names(data.inputR) <- row.names(scobjR)
colnames(data.inputR) <- colnames(scobjR)
metaR = scobjR@meta.data
unique(metaR$FineCellType)
cellchatR <- createCellChat(object = data.inputR, meta = metaR, group.by = "FineCellType")
CellChatDB <- CellChatDB.human
cellchatR@DB <- CellChatDB.human
cellchatR <- subsetData(cellchatR)
cellchatR <- identifyOverExpressedGenes(cellchatR)
cellchatR <- identifyOverExpressedInteractions(cellchatR)
cellchatR <- computeCommunProb(cellchatR, type = "truncatedMean",trim = 0.05, raw.use = TRUE, population.size = TRUE)
cellchatR <- computeCommunProbPathway(cellchatR)
cellchatR <- aggregateNet(cellchatR)
cellchatR <- netAnalysis_computeCentrality(cellchatR, slot.name = "netP")
df.net <- subsetCommunication(cellchatR)
df.netp <- subsetCommunication(cellchatR, slot.name = "netP")
#CircusPlot of R group
groupSize <- as.numeric(table(cellchatR@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchatR@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",
                 sources.use = c('CD4Th','CD8Tcyto','CD4Treg','MacroLA','MacroInf','Fib_RGS5','HCCT_FA','Monocytes','pDC'), 
                 targets.use = c('CD4Th','CD8Tcyto','CD4Treg','MacroLA','MacroInf','Fib_RGS5','HCCT_FA','Monocytes','pDC'))
#Analysis of Non-responders
data.inputNR <- scobjNR[[scobjNR@active.assay]]@layers$data
row.names(data.inputNR) <- row.names(scobjNR)
colnames(data.inputNR) <- colnames(scobjNR)
metaNR = scobjNR@meta.data
unique(metaNR$FineCellType)
cellchatNR <- createCellChat(object = data.inputNR, meta = metaNR, group.by = "FineCellType")
CellChatDB <- CellChatDB.human
cellchatNR@DB <- CellChatDB.human
cellchatNR <- subsetData(cellchatNR)
cellchatNR <- identifyOverExpressedGenes(cellchatNR)
cellchatNR <- identifyOverExpressedInteractions(cellchatNR)
cellchatNR <- computeCommunProb(cellchatNR, type = "truncatedMean",trim = 0.05, raw.use = TRUE, population.size = TRUE)
cellchatNR <- computeCommunProbPathway(cellchatNR)
cellchatNR <- aggregateNet(cellchatNR)
cellchatNR <- netAnalysis_computeCentrality(cellchatNR, slot.name = "netP")
df.net <- subsetCommunication(cellchatNR)
df.netp <- subsetCommunication(cellchatNR, slot.name = "netP")
#CircusPlot of NR group
groupSize <- as.numeric(table(cellchatNR@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchatNR@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",
                 sources.use = c('CD4Th','CD8Tcyto','CD4Treg','MacroLA','MacroInf','Fib_RGS5','HCCT_FA','Monocytes','pDC'), 
                 targets.use = c('CD4Th','CD8Tcyto','CD4Treg','MacroLA','MacroInf','Fib_RGS5','HCCT_FA','Monocytes','pDC'))
#Cross-comparison between R and NR group
object.list <- list(R = cellchatR, NR = cellchatNR)
cellchatMerge <- mergeCellChat(object.list, add.names = names(object.list))
netVisual_diffInteraction(cellchatMerge, weight.scale = T, measure = "weight",
                          sources.use = c('CD4Th','CD8Tcyto','CD4Treg','MacroLA','MacroInf','Fib_RGS5','HCCT_FA','Monocytes','pDC'), 
                          targets.use = c('CD4Th','CD8Tcyto','CD4Treg','MacroLA','MacroInf','Fib_RGS5','HCCT_FA','Monocytes','pDC'))

#Scatter Plot
netAnalysis_signalingRole_scatter(cellchatMerge,x.measure = "outdeg", y.measure = "indeg", 
                                  xlabel = "Outgoing interaction strength", ylabel = "Incoming interaction strength")

#Heatmap Plot
i=1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
#Incoming signals
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(object.list)[i], width = 10, height = 43, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 10, height = 43, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#Outgoing signals
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(object.list)[i], width = 10, height = 43, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(object.list)[i+1], width = 10, height = 43, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#Hierarchy Plot
data.input <- scobj[[scobj@active.assay]]@layers$data
row.names(data.input) <- row.names(scobj)
colnames(data.input) <- colnames(scobj)
meta = scobj@meta.data
unique(meta$FineCellType)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "FineCellType")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean",trim = 0.05, raw.use = TRUE, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
df.net <- subsetCommunication(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
netVisual_aggregate(cellchat, signaling = 'SPP1',layout = "hierarchy",vertex.receiver = vertex.receiver)
netVisual_aggregate(cellchat, signaling = 'APOA',layout = "hierarchy",vertex.receiver = vertex.receiver)
netVisual_aggregate(cellchat, signaling = 'TNF',layout = "hierarchy",vertex.receiver = vertex.receiver)
netVisual_aggregate(cellchat, signaling = 'IFN-II',layout = "hierarchy",vertex.receiver = vertex.receiver)
netVisual_aggregate(cellchat, signaling = 'MIF',layout = "hierarchy",vertex.receiver = vertex.receiver)
netVisual_aggregate(cellchat, signaling = 'FN1',layout = "hierarchy",vertex.receiver = vertex.receiver)
netVisual_aggregate(cellchat, signaling = 'LIGHT',layout = "hierarchy",vertex.receiver = vertex.receiver)
netVisual_aggregate(cellchat, signaling = 'CD86',layout = "hierarchy",vertex.receiver = vertex.receiver)

#Bubble Plot
netVisual_bubble(cellchatMerge,comparison = c(1, 2), angle.x = 45, 
                 sources.use = c('MacroLA','MacroInf','Monocytes','HCCT_MA','HCCT_FA',
                                 'CD4Th','CD8Tcyto','CD4Treg','pDC','CD4Treg','Fib_RGS5'), 
                 targets.use = c('MacroLA','MacroInf','Monocytes','HCCT_MA','HCCT_FA',
                                 'CD4Th','CD8Tcyto','CD4Treg','pDC','CD4Treg','Fib_RGS5'))

#CellCall analysis for heatmap of KEGG enrichment
test <- CreateObject_fromSeurat(Seurat.object=scobj,slot="counts",cell_type="FineCellType",data_source="UMI",
                                scale.factor = 10^6,Org = "Homo sapiens") 
mt <- TransCommuProfile(object = test,pValueCor = 0.05,CorValue = 0.1,topTargetCor=1,p.adjust = 0.05,
                        use.type="mean",probs = 0.9,method="weighted",IS_core = TRUE,Org = 'Homo sapiens')
n2 <- mt@data$expr_l_r_log2_scale
pathway.hyper.list <- lapply(colnames(n2), function(i){print(i)
  tmp <- getHyperPathway(data = n2, object = mt, cella_cellb = i, Org="Homo sapiens")
  return(tmp)})
myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=colnames(n2))
######################################################################################################################################
