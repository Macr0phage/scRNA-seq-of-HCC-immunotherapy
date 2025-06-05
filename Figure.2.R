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
library(ggalluvial)
library(infercnv)
library(Matrix)
library(tidyverse)
library(NMF)
library(scMetabolism)
library(rsvd)
library(ComplexHeatmap)
##############################################################################################################################################
#Fig2
rm(list=ls())
#Loading the 14-patient samples of our project
scobj <- readRDS ("Our single-cell atlas.rds")
#UMAP of All HCC cells
scobjHCC1 = scobj[,scobj@meta.data$MainCellType %in% c('HCC')]
scobjHCC1 <- NormalizeData(scobjHCC1, normalization.method = "LogNormalize", scale.factor = 10000)
scobjHCC1 <- FindVariableFeatures(scobjHCC1, selection.method = "vst", nfeatures = 2500)
scobjHCC1 <- ScaleData(scobjHCC1, features = rownames(scobjHCC1))
scobjHCC1 <- RunPCA(scobjHCC1, features = VariableFeatures(object = scobjHCC1),reduction.name = "pca")
ElbowPlot(scobjHCC1,60)
num=25
scobjHCC1 <- RunHarmony(scobjHCC1,reduction = "pca" ,group.by.vars = 'orig.ident',reduction.save = "harmony", lambda=0.5)
scobjHCC1 <- RunUMAP(scobjHCC1, reduction = "harmony", dims = 1:num, reduction.name = "umap")
scobjHCC1 <- FindNeighbors(scobjHCC1, reduction = "harmony", dims = 1:num)
scobjHCC1 <- FindClusters(scobjHCC1, resolution = 0.2)
DimPlot(scobjHCC1, reduction="umap", pt.size =0.2, repel = TRUE, raster = FALSE,
        group.by=c("seurat_clusters"), alpha = 0.5,
        cols = c('#F9B3AD','#5CB3DD','#70CDBE','#AC99D2','#F5AA61','#81B21F','#6087AD','#BD8566'))

#InferCNV analysis for identification of malignant HCC cells
species = "human"
cnvgroup = "orig.ident"
count_matrix <- scobjHCC1@assays$RNA$counts
cell_metadata <- data.frame(clusters = paste("", scobjHCC1@meta.data$RNA_snn_res.0.2, sep = ""))
cell_metadata[, 1] <- as.character(cell_metadata[, 1])
rownames(cell_metadata) <- rownames(scobjHCC1[[]])
gene_order_file <- read.table("gene_order_human.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = count_matrix, 
                                    annotations_file = cell_metadata, delim = "\t", gene_order_file = gene_order_file, 
                                    ref_group_names = NULL)
dir.create("infercnv_result")
out_dir = "infercnv_result"
infercnv_obj = infercnv::run(infercnv_obj, cutoff = 0.1,
                             out_dir = "infercnv_result", cluster_by_groups = TRUE, plot_steps = FALSE,
                             denoise = F, HMM = FALSE, no_prelim_plot = TRUE,write_phylo = T,write_expr_matrix = T, output_format = "pdf")
cnv_table <- read.table("infercnv_result/infercnv.observations.txt", header=T)
cnv_score_table <- as.matrix(cnv_table)
cnv_score_mat <- as.matrix(cnv_table)
#Scoring
cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < 0.3] <- "A" #complete loss. 2pts
cnv_score_table[cnv_score_mat >= 0.3 & cnv_score_mat < 0.7] <- "B" #loss of one copy. 1pts
cnv_score_table[cnv_score_mat >= 0.7 & cnv_score_mat < 1.3] <- "C" #Neutral. 0pts
cnv_score_table[cnv_score_mat >= 1.3 & cnv_score_mat <= 1.5] <- "D" #addition of one copy. 1pts
cnv_score_table[cnv_score_mat > 1.5 & cnv_score_mat <= 2] <- "E" #addition of two copies. 2pts
cnv_score_table[cnv_score_mat > 2] <- "F" #addition of two copies. 2pts
cnv_score_table_pts <- cnv_table
cnv_score_table_pts[cnv_score_table == "A"] <- 2
cnv_score_table_pts[cnv_score_table == "B"] <- 1
cnv_score_table_pts[cnv_score_table == "C"] <- 0
cnv_score_table_pts[cnv_score_table == "D"] <- 1
cnv_score_table_pts[cnv_score_table == "E"] <- 2
cnv_score_table_pts[cnv_score_table == "F"] <- 2
cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
colnames(cell_scores_CNV) <- "cnv_score"
head(cell_scores_CNV)
write.csv(x = cell_scores_CNV, file = "cnvscore.csv")

#UMAP of malignant HCC cells
scobjHCC2 = scobjHCC1[,scobjHCC1@meta.data$HCC %in% c('Aneuploid')]
scobjHCC2 <- NormalizeData(scobjHCC2, normalization.method = "LogNormalize", scale.factor = 10000)
scobjHCC2 <- FindVariableFeatures(scobjHCC2, selection.method = "vst", nfeatures = 2000)
scobjHCC2 <- ScaleData(scobjHCC2, features = rownames(scobjHCC2))
scobjHCC2 <- RunPCA(scobjHCC2, features = VariableFeatures(object = scobjHCC2),reduction.name = "pca")
ElbowPlot(scobjHCC2,60)
num=25
scobjHCC2 <- RunHarmony(scobjHCC2,reduction = "pca" ,group.by.vars = 'orig.ident',reduction.save = "harmony", lambda=0.5)
scobjHCC2 <- RunUMAP(scobjHCC2, reduction = "harmony", dims = 1:num, reduction.name = "umap")
scobjHCC2 <- FindNeighbors(scobjHCC2, reduction = "harmony", dims = 1:num)
scobjHCC2 <- FindClusters(scobjHCC2, resolution = 0.06)

#NMF analysis of malignant HCC cells in Python
import scanpy as sc
import numpy as np
from cnmf import cNMF
cnmf_obj = cNMF(output_dir="NMF", name="example_cNMF")
cnmf_obj.prepare(counts_fn="subdat.txt", components=np.arange(4,10), n_iter=5, seed=14)
cnmf_obj.factorize(worker_i=0, total_workers=1)
cnmf_obj.combine()
cnmf_obj.k_selection_plot()
cnmf_obj.consensus(k=5, density_threshold=0.01)
usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(K=5, density_threshold=0.01)

#Back to R
NMFcluster <- read.delim("example_cNMF.gene_spectra_score.k_4.dt_0_01.txt",header = T,row.names = 1)
rownames(NMFcluster)<-c("C1","C2","C3","C4")
NMFcluster <- t(NMFcluster)
write.csv(NMFcluster,file = "NMFcluster.csv")
DimPlot(scobjHCC2, reduction="umap", pt.size =0.2, repel = TRUE, raster = FALSE,
        group.by=c("seurat_clusters",'NMFcluster'),label=FALSE,alpha = 0.5,
        cols = c('#A4C3DD','#F9B3AD','#F5AA61','#70CDBE','#AC99D2','#81B21F'))
VlnPlot(scobjHCC2, group.by = 'NMFcluster',features = c('GRN','PLIN2','CTSA'), pt.size = 0, ncol = 3)
VlnPlot(scobjHCC2, group.by = 'NMFcluster',features = c('TOP2A','MKI67','STMN1'), pt.size = 0, ncol = 3)
VlnPlot(scobjHCC2, group.by = 'NMFcluster',features = c('ABCA8','ABCB11','GPAM'), pt.size = 0, ncol = 3)
VlnPlot(scobjHCC2, group.by = 'NMFcluster',features = c('SPARC','IGFBP7','VIM'), pt.size = 0, ncol = 3)

#scMetabolism of malignant HCC cells
new_function1 <- function (obj, method = "VISION", imputation = F, ncores = 2, 
                           metabolism.type = "KEGG") {
  countexp <- obj@assays$RNA$counts
  countexp <- data.frame(as.matrix(countexp))
  signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", 
                                       package = "scMetabolism")
  signatures_REACTOME_metab <- system.file("data", "REACTOME_metabolism.gmt", 
                                           package = "scMetabolism")
  if (metabolism.type == "KEGG") {
    gmtFile <- signatures_KEGG_metab
    cat("Your choice is: KEGG\n")
  }
  if (metabolism.type == "REACTOME") {
    gmtFile <- signatures_REACTOME_metab
    cat("Your choice is: REACTOME\n")
  }
  if (imputation == F) {
    countexp2 <- countexp
  }
  if (imputation == T) {
    cat("Start imputation...\n")
    cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")
    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]
    row.names(countexp2) <- row.names(countexp)}
  cat("Start quantify the metabolism activity...\n")
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2)/n.umi) * median(n.umi)
    vis <- Vision(scaled_counts, signatures = gmtFile)
    options(mc.cores = ncores)
    vis <- analyze(vis)
    signature_exp <- data.frame(t(vis@SigScores))}
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), 
                                           nCores = ncores, plotStats = F)
    geneSets <- getGmt(gmtFile)
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
    signature_exp <- data.frame(getAUC(cells_AUC))}
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile)
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("ssgsea"), 
                    kcdf = c("Poisson"), parallel.sz = ncores)
    signature_exp <- data.frame(gsva_es)}
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile)
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("gsva"), 
                    kcdf = c("Poisson"), parallel.sz = ncores)
    signature_exp <- data.frame(gsva_es)}
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  obj@assays$METABOLISM$score <- signature_exp
  obj}
scobjHCC2 <- new_function1(scobjHCC2)
new_function2 <- function (obj, pathway, dimention.reduction.type = "umap", dimention.reduction.run = T, 
                           size = 1) 
{cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  if (dimention.reduction.type == "umap") {
    if (dimention.reduction.run == T) 
      obj <- Seurat::RunUMAP(obj, reduction = "harmony", dims = 1:40)
    umap.loc <- obj@reductions$umap@cell.embeddings
    row.names(umap.loc) <- colnames(obj)
    signature_exp <- obj@assays$METABOLISM$score
    input.pathway <- pathway
    signature_ggplot <- data.frame(umap.loc, t(signature_exp[input.pathway, 
    ]))
    library(wesanderson)
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    library(ggplot2)
    plot <- ggplot(data = signature_ggplot, aes(x = umap_1, 
                                                y = umap_2, color = signature_ggplot[, 3])) + geom_point(size = size) + 
      scale_fill_gradientn(colours = pal) + scale_color_gradientn(colours = pal) + 
      labs(color = input.pathway) + xlab("UMAP 1") + ylab("UMAP 2") + 
      theme(aspect.ratio = 1) + theme(panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                      axis.line = element_line(colour = "black"))}
  if (dimention.reduction.type == "tsne") {
    if (dimention.reduction.run == T) 
      obj <- Seurat::RunTSNE(obj, reduction = "pca", dims = 1:40)
    tsne.loc <- obj@reductions$tsne@cell.embeddings
    row.names(tsne.loc) <- colnames(obj)
    signature_exp <- obj@assays$METABOLISM$score
    input.pathway <- pathway
    signature_ggplot <- data.frame(tsne.loc, t(signature_exp[input.pathway, 
    ]))
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    library(ggplot2)
    plot <- ggplot(data = signature_ggplot, aes(x = tSNE_1, 
                                                y = tSNE_2, color = signature_ggplot[, 3])) + geom_point(size = size) + 
      scale_fill_gradientn(colours = pal) + scale_color_gradientn(colours = pal) + 
      labs(color = input.pathway) + xlab("tSNE 1") + ylab("tSNE 2") + 
      theme(aspect.ratio = 1) + theme(panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                      axis.line = element_line(colour = "black"))}
  plot}
new_function2(obj = scobjHCC2,                    
              pathway = "Fatty acid degradation",                    
              dimention.reduction.type = "umap",                    
              dimention.reduction.run = F, size = 1.5)
new_function2(obj = scobjHCC2,                    
              pathway = "Glycolysis / Gluconeogenesis",                    
              dimention.reduction.type = "umap",                    
              dimention.reduction.run = F, size = 1.5)

#Pathway enrichment of DEGs among four NMF subtypes of malignant HCC cells
NMFC1markers <- FindMarkers(scobjHCC2, ident.1 =c('C1'), group.by = 'NMF')
genelist <- NMFC1markers$X
x <- genelist
keytypes(org.Hs.eg.db)
ids <- bitr(x, fromType="SYMBOL", toType=c("ENSEMBL","UNIPROT", "ENSEMBL","ENTREZID"), OrgDb="org.Hs.eg.db")
ego <- enrichGO(gene = ids$ENTREZID,universe= names(x),OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",
                pvalueCutoff  = 0.05, readable  = TRUE, qvalueCutoff  = 0.05,)
barplot(ego, drop=TRUE, title = "NMFC1 GO enrichment pathway")
NMFC2markers <- FindMarkers(scobjHCC2, ident.1 =c('C2'), group.by = 'NMF')
genelist <- NMFC1markers$X
x <- genelist
keytypes(org.Hs.eg.db)
ids <- bitr(x, fromType="SYMBOL", toType=c("ENSEMBL","UNIPROT", "ENSEMBL","ENTREZID"), OrgDb="org.Hs.eg.db")
ego <- enrichGO(gene = ids$ENTREZID,universe= names(x),OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",
                pvalueCutoff  = 0.05, readable  = TRUE, qvalueCutoff  = 0.05,)
barplot(ego, drop=TRUE, title = "NMFC2 GO enrichment pathway")
NMFC3markers <- FindMarkers(scobjHCC2, ident.1 =c('C3'), group.by = 'NMF')
genelist <- NMFC1markers$X
x <- genelist
keytypes(org.Hs.eg.db)
ids <- bitr(x, fromType="SYMBOL", toType=c("ENSEMBL","UNIPROT", "ENSEMBL","ENTREZID"), OrgDb="org.Hs.eg.db")
ego <- enrichGO(gene = ids$ENTREZID,universe= names(x),OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",
                pvalueCutoff  = 0.05, readable  = TRUE, qvalueCutoff  = 0.05,)
barplot(ego, drop=TRUE, title = "NMFC3 GO enrichment pathway")
NMFC4markers <- FindMarkers(scobjHCC2, ident.1 =c('C4'), group.by = 'NMF')
genelist <- NMFC1markers$X
x <- genelist
keytypes(org.Hs.eg.db)
ids <- bitr(x, fromType="SYMBOL", toType=c("ENSEMBL","UNIPROT", "ENSEMBL","ENTREZID"), OrgDb="org.Hs.eg.db")
ego <- enrichGO(gene = ids$ENTREZID,universe= names(x),OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",
                pvalueCutoff  = 0.05, readable  = TRUE, qvalueCutoff  = 0.05,)
barplot(ego, drop=TRUE, title = "NMFC4 GO enrichment pathway")
##############################################################################################################################################