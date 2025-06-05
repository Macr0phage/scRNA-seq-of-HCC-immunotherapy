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
#Fig3
rm(list=ls())
#Loading the total single-cell atlas of 75 patient samples
scobj <- readRDS ("Total single-cell atlas.rds")
#UMAP of All T/NK cells
scobjTNK = scobj[,scobj@meta.data$MainCellType %in% c('T/NK')]
scobjTNK <- NormalizeData(scobjTNK, normalization.method = "LogNormalize", scale.factor = 10000)
scobjTNK <- FindVariableFeatures(scobjTNK, selection.method = "vst", nfeatures = 1000)
scobjTNK <- ScaleData(scobjTNK, features = rownames(scobjTNK))
scobjTNK <- RunPCA(scobjTNK, features = VariableFeatures(object = scobjTNK),reduction.name = "pca")
ElbowPlot(scobjTNK,50)
num=20
scobjTNK <- RunHarmony(scobjTNK,reduction = "pca" ,group.by.vars = 'orig.ident',reduction.save = "harmony")
scobjTNK <- RunUMAP(scobjTNK, reduction = "harmony", dims = 1:num, reduction.name = "umap")
scobjTNK <- FindNeighbors(scobjTNK, reduction = "harmony", dims = 1:num)
scobjTNK <- FindClusters(scobjTNK, resolution = 1)
DimPlot(scobjTNK,reduction="umap", pt.size =0.2, repel = TRUE, raster = FALSE,
               group.by=c("seurat_clusters"), alpha = 0.1)

#DotPlot of All T/NK cells
scobjTNK$TNKsubtype = factor(scobjTNK$TNKsubtype, levels =  c('CD4Tnaive','CD8Tnaive','CD4Thelper','CD8Tcyto','CD8Trm',
                                                'CD4Tmemory','CD4Treact','CD8Tex','CD4Treg','CD3Tprolif','NK1','NK2'))
TNKMarkers<-c('CD3G','CD4','CD8A','CCR7','LEF1','SELL','MAL','RACK1','LTB','TNF','GPR183','NR4A1','DUSP2','GZMK','NKG7',
               'GNB2L1','TCEB2','SELK','CD200','FABP5','DUSP4','NR3C1','CXCL13','CCL4','LAG3','PDCD1','CTLA4','TIGIT','LYST',
               'BATF','FOXP3', 'STMN1','MKI67','TOP2A','GNLY','FGFBP2','FCER1G')
DotPlot(scobjTNK,features = TNKMarkers, group.by = "TNKsubtype")+
               RotatedAxis()+RotatedAxis()+scale_color_gradientn(colours = c('#7b95c6','#D2D1D0','#FFCC33'))

#Monocle analysis of CD8+T cells of 14-patient samples in our project
scobjTNK2 = scobjTNK[,scobjTNK@meta.data$Source %in% c('OurProject')]
CD8T = scobjTNK2[,scobjTNK2@meta.data$TNKsubtype %in% c('CD8Tnaive', 'CD8Tcyto', 'CD8Trm', 'CD8Tex')]
expr_matrix<-as(as.matrix(CD8T@assays$RNA$counts), 'sparseMatrix')
p_data <- CD8T@meta.data
f_data <- data.frame(gene_short_name = row.names(CD8T), row.names = row.names(CD8T))
pd <-new('AnnotatedDataFrame', data = p_data)
fd <-new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds,min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed>=3))
diff<-differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = '~TNKsubtype',cores = 1)
deg<-subset(diff, qval<0.05)
deg<-deg[order(deg$qval,decreasing = F),]
ordergene<-rownames(deg)
ordergene <- ordergene[1:900]
cds <- setOrderingFilter(cds,ordergene)
cds@featureData@data[['use_for_ordering']]
plot_ordering_genes(cds)
cds<-reduceDimension(cds,max_components = 2,method='DDRTree')
cds<-orderCells(cds)
plot_cell_trajectory(cds2,color_by = 'State',size=1,show_backbone = TRUE)
plot_cell_trajectory(cds2,color_by = 'TNKsubtype',size=1,show_backbone = TRUE)

#Monocle analysis of CD4+T cells of 14-patient samples in our project
CD4T = scobjTNK2[,scobjTNK2@meta.data$TNKsubtype %in% c('CD4Tnaive', 'CD4Thelper', 'CD4Tmemory', 'CD4Treact', 'CD4Treg')]
expr_matrix<-as(as.matrix(CD4T@assays$RNA$counts), 'sparseMatrix')
p_data <- CD4T@meta.data
f_data <- data.frame(gene_short_name = row.names(CD4T), row.names = row.names(CD4T))
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.1,expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds<-detectGenes(cds,min_expr = 0.1)
print(head(fData(cds)))
expressed_genes<-row.names(subset(fData(cds),num_cells_expressed>=3))
diff<-differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = '~TNKsubtype',cores = 1)
deg<-subset(diff, qval<0.05)
deg<-deg[order(deg$qval,decreasing = F),]
ordergene<-rownames(deg)
ordergene <- ordergene[1:1000]
cds <- setOrderingFilter(cds,ordergene)
cds@featureData@data[['use_for_ordering']]
plot_ordering_genes(cds)
cds<-reduceDimension(cds,max_components = 2,method='DDRTree')
cds<-orderCells(cds)
plot_genes_in_pseudotime(cds['IL7R'],color_by = 'State')
plot_genes_in_pseudotime(cds['IFNG'],color_by = 'State')
plot_genes_in_pseudotime(cds['PDCD1'],color_by = 'State')

#Pathway enrichment of DEGs between Responder/Non-responder within CD4+/CD8+T cells
CD8Tmarkers <- FindMarkers(CD8T, ident.1 =c('Responder'),ident.2 =c('Non-responder'), group.by = 'Response')
genelist <- CD8Tmarkers$X
x <- genelist
keytypes(org.Hs.eg.db)
ids <- bitr(x, fromType="SYMBOL", toType=c("ENSEMBL","UNIPROT", "ENSEMBL","ENTREZID"), OrgDb="org.Hs.eg.db")
ego <- enrichGO(gene = ids$ENTREZID,universe= names(x),OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05, readable  = TRUE, qvalueCutoff  = 0.05,)
barplot(ego, drop=TRUE, title = "CD8+T GO enrichment pathway")
CD4Tmarkers <- FindMarkers(CD4T, ident.1 =c('Responder'),ident.2 =c('Non-responder'), group.by = 'Response')
genelist <- CD4Tmarkers$X
x <- genelist
keytypes(org.Hs.eg.db)
ids <- bitr(x, fromType="SYMBOL", toType=c("ENSEMBL","UNIPROT", "ENSEMBL","ENTREZID"), OrgDb="org.Hs.eg.db")
ego <- enrichGO(gene = ids$ENTREZID,universe= names(x),OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05, readable  = TRUE, qvalueCutoff  = 0.05,)
barplot(ego, drop=TRUE, title = "CD4+T GO enrichment pathway")
##############################################################################################################################################