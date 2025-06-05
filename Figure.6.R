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
#Fig6
rm(list=ls())
#The external spatial transcriptomics dataset from Liu, et al. (PMID: 36708811) can be found on Mendeley Data (skrx2fz79n)
P9T <- readRDS(file = "P9T_Spatial.rds")
P10T <- readRDS(file = "P10T_Spatial.rds")
P9T$Celltype <- P9T@active.ident
P10T$Celltype <- P10T@active.ident
SpatialDimPlot(P9T, label = TRUE, label.size = 3,group.by = c('Celltype'), alpha = 0.5)
SpatialFeaturePlot(P9T, features = "APOA1")
SpatialFeaturePlot(P9T, features = "FABP1")
SpatialFeaturePlot(P9T, features = "C1QA")
P1|P2|P3|P4
SpatialDimPlot(P10T, label = TRUE, label.size = 3,group.by = c('Celltype'), alpha = 0.5)
SpatialFeaturePlot(P10T, features = "APOA1")
SpatialFeaturePlot(P10T, features = "FABP1")
SpatialFeaturePlot(P10T, features = "C1QA")
P1|P2|P3|P4

#The SCENIC analysis of transcription factor within myeloid subtypes
#R
#Loading the 14-patient samples of our project
scobj <- readRDS ("Our single-cell atlas.rds")
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
MacroLA = scobjMye[,scobjMye@meta.data$Myesubtype %in% c('MacroLA')]
exprMat <- as.matrix(GetAssayData(MacroLA,slot ="counts"))
write.csv(t(exprMat),file="scenic.data.csv")

#The following part shall be delivered in Python
import os,sys
import loompy as lp
import numpy as np
import scanpy as sc
import anndata as ad
x=sc.read_csv("scenic.data.csv")
row_attrs={"Gene":np.array(x.var_names),}
col_attrs={"CellID":np.array(x.obs_names)}
lp.create("scenic.loom",x.X.transpose(),row_attrs,col_attrs)
pyscenic grn --num_workers 8 --output grn.tsv --method grnboost2 scenic.loom hgnc_tfs.txt
pyscenic ctx grn.tsv hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname scenic.loom --mode "dask_multiprocessing" --output ctx.csv --num_workers 8 --mask_dropouts
pyscenic aucell scenic.loom ctx.csv --output aucell.loom --num_workers 8

#Back to R
scenicLoomPath=file.path(inputDir,'aucell.loom')
loom <- open_loom(scenicLoomPath) 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])
embeddings <- get_embeddings(loom)  
close_loom(loom)
rownames(regulonAUC)
names(regulons)
sub_regulonAUC <- regulonAUC[,match(colnames(MacroLA),colnames(regulonAUC))]
dim(sub_regulonAUC)
regulonsToPlot = c('ARNTL(+)','E2F1(+)','ELF1(+)','EP300(+)','NFKB1(+)','NFKB2(+)','RUNX2(+)','ESRRA(+)','PPARA(+)','NR1H4(+)')
sce <- MacroLA
sce@meta.data = cbind(sce@meta.data,t(assay(sub_regulonAUC)))
Idents(sce) <- sce$sub_celltype
table(Idents(sce))
DotPlot(sce,features = regulonsToPlot,group.by = "Response",cols = c('lightgrey','red'))+coord_flip()
######################################################################################################################################
