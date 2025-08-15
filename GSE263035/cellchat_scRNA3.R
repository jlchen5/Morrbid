library(CellChat)
library(Seurat)
library(patchwork)  # 用于组合图形
library(dplyr)

# 0. 准备分组数据 ------------------------------------------------------------
# 添加分组信息到元数据
scRNA_tpm@meta.data$condition_group <- paste(scRNA_tpm@meta.data$celltype, 
                                             scRNA_tpm@meta.data$group, 
                                             sep = "_")

# 创建condition_group因子变量
scRNA_tpm@meta.data$condition_group <- factor(scRNA_tpm@meta.data$condition_group)

# 1. 创建函数运行CellChat分析 -------------------------------------------------
run_cellchat_analysis <- function(sce_subset, group_name, db) {
  message("\n运行CellChat分析：", group_name)
  
  # 创建CellChat对象
  cellchat <- createCellChat(sce_subset, 
                             group.by = "celltype", 
                             meta = sce_subset@meta.data)
  
  # 设置数据库
  cellchat@DB <- db
  
  # 预处理
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # 计算通信概率
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # 计算路径级别的通信
  cellchat <- computeCommunProbPathway(cellchat)
  
  # 聚合网络
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}

# 2. 分别对control和stroke组进行分析 -----------------------------------------
# 使用Secreted Signaling数据库
CellChatDB.ss <- subsetDB(CellChatDB.mouse, search = "Secreted Signaling", key = 'annotation')

# 创建分组对象
control_cells <- subset(scRNA_tpm, group == "control")
stroke_cells <- subset(scRNA_tpm, group == "stroke")

# 运行分析
cellchat_control <- run_cellchat_analysis(control_cells, "Control", CellChatDB.ss)
cellchat_stroke <- run_cellchat_analysis(stroke_cells, "Stroke", CellChatDB.ss)

# 3. 比较两组通讯网络 --------------------------------------------------------
# 可视化比较
par(mfrow = c(2, 2))

# 控制组
groupSize_control <- as.numeric(table(cellchat_control@idents))
netVisual_circle(cellchat_control@net$weight, 
                 vertex.weight = groupSize_control,
                 weight.scale = T,
                 vertex.label.cex = 0.8,
                 title.name = "Control Group Interaction Strength")
# 整合的互作对数量展示
netVisual_circle(cellchat_control@net$count, vertex.weight = groupSize_control, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Control Group Number of interactions")
# 中风组
groupSize_stroke <- as.numeric(table(cellchat_stroke@idents))
netVisual_circle(cellchat_stroke@net$weight, 
                 vertex.weight = ,
                 weight.scale = T,
                 vertex.label.cex = 0.8,
                 title.name = "Stroke Group Interaction Strength")

netVisual_circle(cellchat_stroke@net$count, vertex.weight = groupSize_stroke, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Stroke Group Number of interactions")


# 4. 差异分析 ----------------------------------------------------------------
# 合并结果进行比较
cellchat.list <- list(Control = cellchat_control, Stroke = cellchat_stroke)
cellchat.merged <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))

# 比较相互作用数量
gg1 <- compareInteractions(cellchat.merged, show.legend = F, group = c(1, 2))
gg2 <- compareInteractions(cellchat.merged, show.legend = F, group = c(1, 2), measure = "weight")
gg1 + gg2  # 并列显示两个图

# 差异热图
par(mfrow = c(1, 1))
netVisual_heatmap(cellchat.merged)

# 识别差异最大的信号通路
rankNet(cellchat.merged, mode = "comparison")

# 5. 细胞类型特异性分析 ------------------------------------------------------
# 定义可视化函数
plot_celltype_specific <- function(cellchat_obj, group_name) {
  weight_mat <- cellchat_obj@net$weight
  groupSize <- as.numeric(table(cellchat_obj@idents))
  
  par(mfrow = c(3, 3), mar = c(1, 1, 2, 1))
  
  for (cel in unique(cellchat_obj@idents)) {
    cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
    cir_mat[cel, ] <- weight_mat[cel, ]
    
    netVisual_circle(
      cir_mat,
      vertex.weight = groupSize,
      weight.scale = T,
      edge.weight.max = max(weight_mat),
      vertex.label.cex = 0.6,
      title.name = paste(group_name, "-", cel)
    )
  }
}

# 绘制各组细胞类型特异性通讯
plot_celltype_specific(cellchat_control, "Control")
plot_celltype_specific(cellchat_stroke, "Stroke")

