rm(list=ls())
library(Seurat)


ct=fread("scRNA-seq/TAC_raw_umi_matrix.csv",data.table = F)
rownames(ct) <- ct$V1
ct <- ct[,-1]
head(ct)[1:5,1:2]

sce=CreateSeuratObject(counts =  ct ,
                       project =  sc ,
                       min.cells = 5,
                       min.features = 300,)
dim(sce)

sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")

##### 标准化数据 #####
sce <- NormalizeData(sce, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000, verbose = FALSE)
##### 识别高度可变的特征（特征选择）#####
sce <- FindVariableFeatures(sce, selection.method = "vst", 
                            nfeatures = 2000, verbose = FALSE)



# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sce), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sce) + 
  theme(legend.position="top")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + 
  theme(legend.position="none")

##### 缩放数据 #####
all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes, verbose = FALSE)

##### 线性降维 #####
sce <- RunPCA(sce, features = VariableFeatures(object = sce), verbose = FALSE)
VizDimLoadings(sce, dims = 1:2, reduction = "pca")
DimPlot(sce, reduction = "pca")

##### 细胞聚类 ##### 
sce <- FindNeighbors(sce, dims = 1:20, verbose = FALSE)
sce <- FindClusters(sce, resolution = 0.5, verbose = FALSE)


##### 运行非线性降维 UMAP/tSNE #####

sce <- RunUMAP(sce, dims = 1:20, verbose = FALSE)
DimPlot(sce, reduction = "umap")

sce <- RunTSNE(sce, dims = 1:20, verbose = FALSE)
DimPlot(sce, reduction = "tsne")

DimPlot(sce, reduction = "umap", label = TRUE)


# save data
saveRDS(sce, file = "mmNHT_processed.rds")


