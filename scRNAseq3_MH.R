# 加载必要的包
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(celldex)
library(SingleR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

# 1. 读取数据
heart <- Read10X(data.dir = "./scRNAseq3_MouseHeart/scdata_MH")
heart <- CreateSeuratObject(counts = heart, project = "ccr2tdt-heart", min.cells = 3, min.features = 200)
table(heart$orig.ident) 

heart[["percent.mt"]] <- PercentageFeatureSet(heart, pattern = "^mt-")
head(heart)
sample <- sapply(colnames(GetAssayData(object = heart, slot = "counts")),FUN=function(x){substr(x,18,19)})

sample <- as.numeric(as.character(sample))
names(sample) <- colnames(GetAssayData(object = heart, slot = "counts"))
heart <- AddMetaData(heart, sample, "sample")
table(heart$orig.ident) 

VlnPlot(object = heart, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(heart, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(heart, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

heart_filtered <- subset(heart, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
VlnPlot(object = heart_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

heart_filtered <- NormalizeData(heart_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
heart_filtered <- FindVariableFeatures(heart_filtered, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(heart_filtered), 10)
all.genes <- rownames(heart_filtered)
heart_filtered <- ScaleData(heart_filtered, features = all.genes)
heart_filtered <- RunPCA(heart_filtered, features = VariableFeatures(object = heart_filtered))
DimHeatmap(heart_filtered, dims = 5:15, cells = 500, balanced = TRUE)

heart_filtered <- FindNeighbors(heart_filtered, dims = 1:13)
heart_filtered <- FindClusters(heart_filtered, resolution = 0.4)
heart_filtered <- RunUMAP(heart_filtered, dims = 1:13)

library(RColorBrewer)
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Set1", n = 8), brewer.pal(name="Accent", n = 6))
DimPlot(heart_filtered, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8)
DimPlot(heart_filtered, reduction = "umap", split.by =  "sample", group.by = "sample")

heartok <- heart_filtered 
new.grouping <- c("group")
heartok[[new.grouping]] <- new.grouping
colnames(heartok@meta.data)

heartok$group[heartok$sample == "2"] <- "stroke"
heartok$group[heartok$sample == "1"] <- "control"
heartok$group[heartok$sample == "3"] <- "control"
heartok$group[heartok$sample == "4"] <- "stroke"
head(heartok@meta.data, 10)

markers <- FindAllMarkers(object = heartok, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
write.table((markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)), './scRNAseq3_MouseHeart/markers_by_cluster.tsv', sep='\t')


# 细胞类型注释
ref <- MouseRNAseqData()
annotations <- SingleR(test = GetAssayData(heartok, assay = "RNA"), 
                       ref = ref, labels = ref$label.main)
heartok$celltype <- annotations$labels

# 数据可视化
features <- c("Ccl25","Il22ra1","Il24","Tslp")
RidgePlot(heartok, features = features, group.by = "celltype", ncol = 3)
VlnPlot(heartok, features = features, group.by = "celltype", ncol = 3)
FeaturePlot(heartok, features = features, cols = c("grey", "blue"), combine = TRUE,, ncol = 3)



# 功能富集分析
top_markers <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 200, wt = avg_log2FC)

# 按cluster拆分基因列表
gene_list <- split(top_markers$gene, top_markers$cluster)

# 1. GO富集分析 ------------------------------------------------------------
go_results <- list()

for(cluster_id in names(gene_list)){
  # 基因符号转Entrez ID
  entrez_ids <- bitr(gene_list[[cluster_id]], 
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Mm.eg.db)
  
  # 去除重复的Entrez ID
  unique_entrez <- unique(entrez_ids$ENTREZID)
  
  # GO分析）
  ego <- enrichGO(gene          = unique_entrez,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2)
  
  go_results[[cluster_id]] <- ego
}

# 可视化GO结果（示例展示第一个cluster）
clusterProfiler::dotplot(ego, showCategory=10) + 
  ggtitle("GO Enrichment Analysis")


# 2. KEGG通路分析 ----------------------------------------------------------
kegg_results <- list()

for(cluster_id in names(gene_list)){
  # 基因符号转Entrez ID
  entrez_ids <- bitr(gene_list[[cluster_id]], 
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Mm.eg.db)
  
  # KEGG分析
  kk <- enrichKEGG(gene         = unique(entrez_ids$ENTREZID),
                   organism     = 'mmu',
                   pvalueCutoff = 0.05)
  
  kegg_results[[cluster_id]] <- kk
}

# 可视化KEGG结果（示例展示第一个cluster）
dotplot(kk, showCategory=10) + 
  ggtitle("KEGG Pathways")



# 提取标准化表达数据
exp_data <- FetchData(heartok, vars = c(features, "group"))

# 长格式转换
plot_data <- exp_data %>%
  tidyr::pivot_longer(cols = -group, names_to = "gene", values_to = "expression")

ggplot(plot_data, aes(x = gene, y = expression, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.8),outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             size = 0.1, alpha = 0.3) +
  ggpubr::stat_compare_means(aes(group = group), 
                             method = "wilcox.test",
                             label = "p.format",
                             label.y = max(plot_data$expression) * 1.1) +
  labs(title = "Gene Expression by Group", y = "Log-Normalized Expression") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



