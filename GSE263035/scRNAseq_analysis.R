# 加载必要的包
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(monocle3)
library(CellChat)
library(ggraph)

library(celldex)
library(SingleR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(SCINA)

## load data
heartok <- readRDS('./heartok.rds')
DimPlot(heartok, group.by = "celltype",split.by = "group",label = TRUE) 

# UMAP
DimPlot(heartok, reduction = "umap", label=TRUE)

# DimPlot(heartok, 
#         reduction = "umap",
#         group.by = "orig.ident",  # 按样本来源分组
#         pt.size = 0.5) +
#   ggtitle("Control vs Stroke")

# tSNE
DimPlot(sce, reduction = "tsne", label=TRUE)

sce <- heartok

# 为每个细胞类型寻找 marker genes
sce <- JoinLayers(sce) # 合并数据层
Idents(sce) <- "celltype"  # 设置细胞类型为当前标识
all.markers <- FindAllMarkers(sce,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

# 为每个细胞类型选择 top 5 marker genes
top_markers <- all.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

# 创建 dot plot
celltype_order <- names(sort(table(sce$celltype), decreasing = TRUE))
sce$celltype <- factor(sce$celltype, levels = celltype_order)
DotPlot(sce,features = unique(top_markers$gene),
            group.by = "celltype",
            cols = c("lightblue", "darkblue"),
            dot.scale = 6,
            cluster.idents = TRUE) +  # 启用聚类排序
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 11)) +
  scale_size_continuous(name = "Percent\nExpressed",
                        range = c(1, 8),
                        breaks = c(25, 50, 75)) +
  scale_color_gradient(name = "Average\nExpression",
                       low = "lightblue",
                       high = "darkblue") 

# 计算各细胞类型在两组中的比例
celltype_counts <- as.data.frame(table(sce$group, sce$celltype))
colnames(celltype_counts) <- c("Condition", "CellType", "Count")

celltype_proportions <- celltype_counts %>%
  group_by(Condition) %>%
  mutate(Proportion = Count / sum(Count) * 100)

# 创建堆叠条形图
ggplot(celltype_proportions, aes(x = Condition, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(celltype_proportions$CellType)))) +
  theme_classic() +
  theme(panel.background = element_blank(),
        legend.position = "right",
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 12)) +
  labs(x = "Experimental Group", 
       y = "Percentage of Cells", 
       fill = "Cell Type") +
  guides(fill = guide_legend(ncol = 1)) +
  geom_text(aes(label = paste0(round(Proportion, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 3, color = "black", alpha = 0.9)


##########################
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggpubr)  # 用于添加统计学标注

# Feature Plot 
# ---------------------------------------------
# 设置配色方案
gradient_colors <- c("lightgrey", "purple4")

# 为每个基因和每个组别创建feature plot
feature_plots <- lapply(c("Ccl25","Il22ra1","Tslp"), function(gene) {
  plot_control <- FeaturePlot(sce, 
                              features = gene, 
                              split.by = NULL, 
                              reduction = "umap", 
                              pt.size = 0.2,
                              cols = gradient_colors) +
    ggtitle(paste("Control:", gene)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  
  plot_stroke <- FeaturePlot(sce, 
                             features = gene, 
                             split.by = NULL, 
                             reduction = "umap", 
                             pt.size = 0.2,
                             cols = gradient_colors) +
    ggtitle(paste("Stroke:", gene)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  
  plot_control + plot_stroke
})

# 组合所有feature plot
wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Gene Expression in UMAP Space",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

# 保存结果
# ggsave("feature_plots.png", feature_grid, width = 16, height = 8, dpi = 300)

# 小提琴图
library(ggpubr)
vln_data <- FetchData(sce, vars = c("celltype", "group", "Ccl25", "Il22ra1", "Tslp")) %>%
  tidyr::pivot_longer(cols = c("Ccl25", "Il22ra1", "Tslp"),
                      names_to = "Gene", 
                      values_to = "Expression")

# 数据预处理
# ---------------------------
# 确保celltype是因子且水平正确
vln_data$celltype <- factor(vln_data$celltype, 
                                 levels = sort(unique(vln_data$celltype)))

# 过滤掉数据点不足的组合
vln_data_filtered <- vln_data %>%
  group_by(celltype, Gene, group) %>%
  filter(n() >= 3) %>%  # 每组至少3个数据点
  ungroup()

# 创建分面小提琴图
# ---------------------------
ggplot(vln_data_filtered, aes(x = celltype, y = Expression, fill = group)) +
  # 小提琴图主体
  geom_violin(scale = "width", 
              position = position_dodge(width = 0.9), 
              trim = TRUE, 
              width = 0.8,
              alpha = 0.7) +
  
  # 箱形图叠加
  geom_boxplot(width = 0.15, 
               position = position_dodge(width = 0.9),
               outlier.size = 0.5, 
               show.legend = FALSE,
               alpha = 0.5) +
  
  # 分面设置
  facet_wrap(~ Gene, nrow = 1, scales = "free_y") +
  
  # 统计检验（仅在有足够数据的组间进行）
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = F,
    label.y = max(vln_data_filtered$Expression, na.rm = TRUE) * 0.95,
    size = 4,
    vjust = 0.5,
    show.legend = F) +
  
  # 颜色和主题设置
  scale_fill_manual(values = c("control" = "grey70", "stroke" = "#E69F00")) +
  labs(x = "Cell Type", 
       y = "Expression Level", 
       fill = "Group",
       title = "Gene Expression Across Cell Types") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.spacing = unit(1, "lines")
  )

# ggsave("violin_plots_optimized.pdf", vln_plot, width = 16, height = 8, dpi = 300)


# Dot Plot 
# ---------------------------------------------
# 创建分组变量（细胞类型+条件）
sce$celltype_condition <- paste0(sce$celltype, "_", sce$group)

# 创建dot plot
DotPlot(sce, features = c("Ccl25","Il22ra1","Tslp"), dot.scale = 6) +
  scale_color_gradientn(colors = c("lightgrey", "navy")) +
  scale_size_continuous(range = c(1, 8)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right") +
  labs(x = "Genes", y = "Cell Types and Conditions",
       color = "Average Expression", size = "Percent Expressed") 

# # 添加条件标签
# cell_conditions <- data.frame(
#   celltype_condition = unique(sce$celltype_condition),
#   Condition = ifelse(grepl("Control", unique(sce$celltype_condition)), "Control", "Stroke")
# )
# 
# dot_plot$data <- left_join(dot_plot$data, cell_conditions, by = "celltype_condition")
# dot_plot
# 保存结果
# ggsave("dot_plot.png", dot_plot, width = 12, height = 8, dpi = 300)



# 拟时序分析 (monocle3) ---------------------------
library(dplyr)
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(SeuratObject)
library(SeuratWrappers)
# 导入注释好的Seurat对象
# heartok <- readRDS('./heartok.rds')
scRNA_tpm <- heartok

# 清理内存
gc(full = TRUE, reset = TRUE)  # 强制释放内存
rm(list = ls()[!ls() %in% "scRNA_tpm"])  # 只保留必需对象

library(monocle3)
library(SeuratWrappers)
library(ggplot2)

# 创建CellDataSet
cds <- as.cell_data_set(scRNA_tpm)
unique(colData(cds)$celltype)  # 确认列名和内容

# 提取成纤维细胞
fb_cell <- colnames(cds)[colData(cds)$celltype == "Fibroblasts"]
cds_fd <- cds[, fb_cell]

# 数据降维
cds_fd <- reduce_dimension(cds_fd, reduction_method = "UMAP")
cds_fd <- cluster_cells(cds_fd)
cds_fd <- learn_graph(cds_fd)

# 重新指定根节点（例如选择某个高表达基因的细胞）
if ("Tslp" %in% rownames(cds_fd)) {
  tslp_expr <- exprs(cds_fd)["Tslp", ]
  progenitor_cells <- names(tslp_expr[tslp_expr > quantile(tslp_expr, 0.95)])
  cds_fd <- order_cells(cds_fd, root_cells = progenitor_cells)
} else {
  # 备选方案：手动选择集群
  plot_cells(cds_fd, color_cells_by = "cluster")
  cds_fd <- order_cells(cds_fd) # 交互式选择
}

# 可视化成纤维细胞的轨迹
plot_cells(cds_fd,
           color_cells_by = "group",
           label_cell_groups = F,
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           cell_size=1,
           graph_label_size = 2,
           trajectory_graph_color='black',
           trajectory_graph_segment_size = 1)
           

################### 绘制基因表达随伪时间变化的趋势图 ######################
# 1. 确保基因存在并创建子集对象
genes_to_plot <- c("Tslp", "Il22ra1", "Ccl25")

# 检查基因是否存在
available_genes <- genes_to_plot[genes_to_plot %in% rownames(cds_fd)]

if (length(available_genes) == 0) {
  # 尝试大小写不敏感的匹配
  all_genes <- rownames(cds_fd)
  available_genes <- unique(unlist(sapply(genes_to_plot, function(g) {
    grep(paste0("^", g, "$"), all_genes, ignore.case = TRUE, value = TRUE)
  })))
}

if (length(available_genes) == 0) {
  stop("指定的基因在数据集中不存在。请检查基因名称。")
}

# 创建包含这些基因的子集对象
cds_subset <- cds_fd[available_genes, ]

# 2. 添加基因符号作为gene_short_name（如果不存在）
if (!"gene_short_name" %in% colnames(rowData(cds_subset))) {
  rowData(cds_subset)$gene_short_name <- rownames(cds_subset)
}

# 3. 绘制伪时间表达趋势图
plot_genes_in_pseudotime(
  cds_subset,
  min_expr = 0.1,  # 降低阈值以显示更多数据点
  ncol = min(length(available_genes), 2),
  cell_size = 1.5,  # 点的大小
  color_cells_by = "group",  # 按实验组着色
  trend_formula = "~ splines::ns(pseudotime, df=3)",  # 添加平滑曲线
  label_by_short_name = TRUE
) +
  labs(title = "基因表达随伪时间变化",
       x = "伪时间",
       y = "标准化表达量") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12))

# 4. 分面显示每个基因的表达
library(ggplot2)
library(dplyr)

# 提取伪时间和表达数据
pseudotime_values <- pseudotime(cds_subset)
expression_data <- exprs(cds_subset)[available_genes, ]

# 创建绘图数据框
plot_data <- data.frame(
  Pseudotime = rep(pseudotime_values, length(available_genes)),
  Expression = as.vector(expression_data),
  Gene = rep(available_genes, each = ncol(cds_subset)),
  Group = rep(colData(cds_subset)$group, length(available_genes))
)

# 绘制分面图
ggplot(plot_data, aes(x = Pseudotime, y = Expression, color = Group)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("control" = "blue", "stroke" = "red")) +
  labs(title = "",
       x = "Presudotime",
       y = "Relative Expression",
       color = "Group") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"))

# 5. 在轨迹图上可视化基因表达
plot_cells(cds_fd,
           genes = available_genes,
           show_trajectory_graph = TRUE,
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           cell_size = 1,
           trajectory_graph_color = "grey50",
           trajectory_graph_segment_size = 0.5) +
  facet_wrap(~ feature_label, ncol = min(length(available_genes), 3)) +
  scale_color_viridis_c(option = "C") +
  labs(title = " ") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

################### 绘制gene表达趋势图 ######################
# 1. 检查现有基因元数据
colnames(rowData(cds_fd))  # 查看可用的列名

# 2. 添加gene_short_name列（如果不存在）
if(!"gene_short_name" %in% colnames(rowData(cds_fd))) {
  # 假设目前的rownames就是基因ID
  rowData(cds_fd)$gene_short_name <- rownames(cds_fd)
  
  # 或者手动添加标准基因符号（示例）
  # rowData(cds_fd)$gene_short_name <- ifelse(
  #   rownames(cds_fd) == "some_id", "Tslp", rownames(cds_fd)
  # )
}

# 3. 查找基因的标准符号（关键步骤）
# 将您的基因ID转换为标准基因符号
all_genes <- rownames(cds_fd)

# 查找类似基因（不区分大小写）
tslp_id <- grep("^tslp$", all_genes, ignore.case = TRUE, value = TRUE)
ccl25_id <- grep("^ccl25$", all_genes, ignore.case = TRUE, value = TRUE)
il22ra1_id <- grep("^il22ra1$", all_genes, ignore.case = TRUE, value = TRUE)

# 如果找到，使用实际ID
gene_ids <- c(tslp_id, ccl25_id, il22ra1_id)
gene_ids <- unique(gene_ids[!is.na(gene_ids) & nchar(gene_ids) > 0])

# 4. 绘制基因表达趋势图
if(length(gene_ids) > 0) {
  # 方法1：使用实际找到的基因ID
  plot_genes_in_pseudotime(cds_fd[gene_ids,], color_cells_by = "group",
                           min_expr = 1e-05,  # 降低阈值以提高灵敏度
                           ncol = 3,
                           label_by_short_name = FALSE)  # 使用原始ID作为标签
  
  # 方法2：添加友好的标签
  for_plot <- cds_fd[gene_ids,]
  rowData(for_plot)$custom_label <- c("Tslp", "Ccl25", "Il22ra1")[1:length(gene_ids)]
  
  plot_genes_in_pseudotime(for_plot, color_cells_by = "group",
                           min_expr = 1e-05,
                           ncol = min(length(gene_ids), 2),
                           label_by_short_name = FALSE) +
    facet_wrap(~custom_label)
} 
