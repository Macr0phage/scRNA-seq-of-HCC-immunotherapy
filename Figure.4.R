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
##############################################################################################################################################
#Fig4
rm(list=ls())
#Loading the total single-cell atlas of 75 patient samples
scobj <- readRDS ("Total single-cell atlas.rds")
#UMAP of All Myeloid cells
scobjMye = scobj[,scobj@meta.data$MajorCellType %in% c('Myeloid')]
scobjMye <- NormalizeData(scobjMye, normalization.method = "LogNormalize", scale.factor = 10000)
scobjMye <- FindVariableFeatures(scobjMye, selection.method = "vst", nfeatures = 1500)
scobjMye <- ScaleData(scobjMye, features = rownames(TotalMye2))
scobjMye <- RunPCA(scobjMye, features = VariableFeatures(object = TotalMye2),reduction.name = "pca")
ElbowPlot(scobjMye,60)
num=25
scobjMye <- RunHarmony(scobjMye,reduction = "pca" ,group.by.vars = 'orig.ident',reduction.save = "harmony")
scobjMye <- RunUMAP(scobjMye, reduction = "harmony", dims = 1:num, reduction.name = "umap", min.dist = 0.5)
scobjMye <- FindNeighbors(scobjMye, reduction = "harmony", dims = 1:num)
scobjMye <- FindClusters(scobjMye, resolution = 0.1)
DimPlot(scobjMye,reduction="umap", pt.size =0.3, repel = TRUE, raster = FALSE,
        group.by=c('Myesubtype'),label=FALSE,alpha = 0.1,split.by = 'Origin')
DimPlot(scobjMye,reduction="umap", pt.size =0.3, repel = TRUE, raster = FALSE,
        group.by=c('Response'),label=FALSE,alpha = 0.1)

#DotPlot of All Myeloid cells
scobjMye$Myesubtype = factor(scobjMye$Myesubtype, levels =  c('MacroAngio','MacroProlif','MacroLA',
                                                           'MacroTR','MacroInflam','MacroIFN','Neutrophils','Mast','Monocytes',
                                                           'cDC1','cDC2','cDC3','pDC'))
MyeMarkers<-c('C1QA','CD68','GPNMB','TREM2','SPP1','TOP2A','MKI67','APOB','APOC2','APOC3',
                              'VCAM1','SLC40A1','CXCL2','CCL3L1','CXCL10','GBP1','IFITM1','IL32','IL7R','S100A8',
                              'TREM1','CPA3','TPSB2','VCAN','EREG','CLEC9A','CADM1','CD1C','CLEC10A','LAMP3','CCR7','GZMB','LILRA4')
DotPlot(scobjMye,features = MyeMarkers, group.by = "Myesubtype")+
               RotatedAxis()+RotatedAxis()+scale_color_gradientn(colours = c('#7b95c6','#D2D1D0','#FFCC33'))

#Monocle analysis of monocytes and macrophages of 14-patient samples in our project
MonoMacro = scobjMye[,scobjMye@meta.data$Myesubtype %in% c('MacroAngio','MacroProlif','MacroLA','MacroTR','MacroInflam','MacroIFN','Monocytes')]
expr_matrix<-as(as.matrix(MonoMacro@assays$RNA$counts), 'sparseMatrix')
p_data <- MonoMacro@meta.data
f_data<-data.frame(gene_short_name = row.names(MonoMacro), row.names = row.names(MonoMacro))
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds,min_expr = 0.1)
expressed_genes<-row.names(subset(fData(cds),num_cells_expressed>=3))
diff<-differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = '~Myesubtype',cores = 1)
deg<-subset(diff, qval<0.05)
deg<-deg[order(deg$qval,decreasing = F),]
ordergene<-rownames(deg)
ordergene<-row.names(deg)[order(deg$qval)][1:500]
cds<-setOrderingFilter(cds,ordergene)
cds@featureData@data[['use_for_ordering']]
plot_ordering_genes(cds)
cds<-reduceDimension(cds,max_components = 2,method='DDRTree')
cds<-orderCells(cds)
plot_cell_trajectory(cds,color_by = 'Myesubtype',size=1,show_backbone = TRUE)
plot_cell_trajectory(cds,color_by = 'State',size=1,show_backbone = TRUE)
plot_genes_in_pseudotime(cds['CXCL8'],color_by = 'State')
plot_genes_in_pseudotime(cds['APOC1'],color_by = 'State')
plot_genes_in_pseudotime(cds['CD86'],color_by = 'State')

#DEG analysis of among three states of monocytes and macrophages from Monocle
State2_1markers <- FindMarkers(MonoMacro, ident.1 =c('State2'),ident.2 =c('State1'), group.by = 'State')
State3_1markers <- FindMarkers(MonoMacro, ident.1 =c('State3'),ident.2 =c('State1'), group.by = 'State')
State3_2markers <- FindMarkers(MonoMacro, ident.1 =c('State3'),ident.2 =c('State2'), group.by = 'State')
jjVolcano(diffData = JJPatientDEGs, log2FC.cutoff = 0.75, size  = 3.5)

#Pathway enrichment of DEGs among three states of monocytes and macrophages from Monocle
genelist <- State2markers$X
x <- genelist
keytypes(org.Hs.eg.db)
ids <- bitr(x, fromType="SYMBOL", toType=c("ENSEMBL","UNIPROT", "ENSEMBL","ENTREZID"), OrgDb="org.Hs.eg.db")
kk<-enrichKEGG(gene=ids$ENTREZID, organism="hsa", pvalueCutoff=0.05, qvalueCutoff=0.05)
barplot(kk, drop=TRUE, title = "State2 KEGG enrichment pathway")
genelist <- State3markers$X
x <- genelist
keytypes(org.Hs.eg.db)
ids <- bitr(x, fromType="SYMBOL", toType=c("ENSEMBL","UNIPROT", "ENSEMBL","ENTREZID"), OrgDb="org.Hs.eg.db")
kk<-enrichKEGG(gene=ids$ENTREZID, organism="hsa", pvalueCutoff=0.05, qvalueCutoff=0.05)
barplot(kk, drop=TRUE, title = "State3 KEGG enrichment pathway")
##############################################################################################################################################