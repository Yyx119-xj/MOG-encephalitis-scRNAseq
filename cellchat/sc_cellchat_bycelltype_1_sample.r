#!/usr/bin/env Rscript

### Set VERSION
VERSION <- "1.0.0"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1 && (args[1] == "-v" | args[1] == "--version")) {
  message("Version: \n\t", VERSION, "\n")
  quit(status = 0)
}

### 修改记录
####################################################################################
# Date           Version       Author            ChangeLog

#####################################################################################
#####################################################################################
##### 参考说明
# https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
# https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
#####################################################################################




#####################################################################################
##### 参数获取
#####################################################################################
# windows系统检查当前程序路径
# script.dir <- dirname(sys.frame(1)$ofile)

library(tidyverse)
get_this_file <- function() {
  commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col = value, into = c("key", "value"), sep = "=", fill = "right") %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
}
script.dir <- get_this_file()

setwd("./")

suppressPackageStartupMessages(library(optparse)) ## Options
suppressPackageStartupMessages(library(futile.logger)) ## logger
suppressPackageStartupMessages(library("stats"))
suppressPackageStartupMessages(library(dplyr))
`%ni%` <- Negate(`%in%`)

library("ggplot2")
options(stringsAsFactors = FALSE)
library(grDevices)
library(gplots)
library(pheatmap)
library(Seurat)
library(patchwork)
library(SeuratWrappers)
library(cowplot)
library("future.apply")
suppressPackageStartupMessages(library(doParallel))
library(dplyr)
library(CellChat)
library(NMF)
library(ggalluvial)
library(ComplexHeatmap)

## 参数读取
option_list <- list(
  make_option(
    c("-n", "--ncpus"),
    action = "store",
    type = "integer",
    default = 5,
    help = "parallel run number of cores"
  ),
  make_option(
    c("--MaxMemMega"),
    action = "store",
    type = "integer",
    default = 80000,
    help = "Parallel Max Memory size megabytes"
  ),
  make_option(
    c("-i", "--infile"),
    action = "store",
    type = "character",
    help = "The cell analysis rds file"
  ),
  make_option(
    c("-g", "--groupname"),
    action = "store",
    type = "character",
    help = "The groupname name"
  ),
  make_option(
    c("-s", "--species"),
    action = "store",
    type = "character",
    default = "human",
    help = "human or mouse"
  ),
  make_option(
    c("-k", "--keep"),
    action = "store",
    type = "character",
    default = "",
    help = "ex: 3,5,7 clusters need to keep, default is all. "
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

registerDoParallel(cores = opt$ncpus)
plan("multicore", workers = opt$ncpus)
options(future.globals.maxSize = opt$MaxMemMega * 1024^2)
gname = opt$groupname

#####################################################################################

# 直接通过seurat生成的rds文件的方式进行转换为CellDataSet对象
sc.combined <- readRDS(file = opt$infile)

#剔除掉任意细胞数量小于10的类群
sc.combined$celltype <- paste0(sc.combined$group,"_",sc.combined$customclassif)
celltype_counts <- table(sc.combined$celltype)
celltypes_to_remove <- names(celltype_counts[celltype_counts < 10])
sc.combined <- subset(sc.combined, cells = colnames(sc.combined)[!sc.combined$celltype %in% celltypes_to_remove])

#
sc.combined$cellcluster <- sc.combined$customclassif

CellTalk_dir <- paste("./", "", sep = "")
if (!file.exists(CellTalk_dir)) {
  dir.create(CellTalk_dir)
}
setwd(CellTalk_dir)

# 从seurat对象读取表达量
data <- sc.combined@assays$RNA@data
head(data)

# 替换
original = as.matrix(sc.combined$cellcluster)
# 构建cellchat对象需要知道每个细胞对应的细胞类型和分组
meta <- data.frame(cell_type = as.character(original), cell_group = sc.combined@meta.data$group)
rownames(meta) <- rownames(sc.combined@meta.data)
head(meta)
# 提取某一组的细胞用于后续的分析
cell.use <- rownames(meta)[meta$cell_group == gname] 
data.input <- data[, cell.use]
meta.input <- meta[cell.use, ]
unique(meta.input$cell_group) 
unique(meta.input$cell_type) 
# 创建cellchat对象
cellchat.gname <- createCellChat(object = data.input, meta = meta.input, group.by = "cell_type")
table(cellchat.gname@idents) 

unique(cellchat.gname@idents) # 确认细胞类型的顺序，如果需要做某些指定target的图，需要在这一步进行确认
# 创建文件夹存放结果
setwd(CellTalk_dir)
tmpdir <- paste("", "sample_",gname, sep = "")
if (!file.exists(tmpdir)) {
  dir.create(tmpdir)
}
setwd(tmpdir)

#### 开始作图

#确定使用的物种，目前cellchatv2只支持人和小鼠
CellChatDB <- CellChatDB.human 
if(opt$species == "mouse"){
  CellChatDB <- CellChatDB.mouse
}

CellChatDB.use <- subsetDB(CellChatDB) # 使用除了Non-protein Signaling以外的信号数据
cellchat.gname@DB <- CellChatDB.use


cellchat.gname <- subsetData(cellchat.gname) # 提取表达量
future::plan("multicore", workers = opt$ncpus) # 默认为5，应该可以适当调高

cellchat.gname <- identifyOverExpressedGenes(cellchat.gname)
cellchat.gname <- identifyOverExpressedInteractions(cellchat.gname) #两个官方默认参数运行即可，有自定义选项

#if(opt$species == "human"){
#  cellchat.gname <- projectData(cellchat.gname, PPI.human)
#}
#if(opt$species == "mouse"){
#  cellchat.gname <- projectData(cellchat.gname, PPI.mouse)
#}
# 这个不需要投影吧，而且原脚本后面也没有加对应的参数调整啊

cellchat.gname <- computeCommunProb(cellchat.gname) # 如果是确定数据内零值较多，则需要将上面的注释去掉，并在这里加个参数raw.use = FALSE
cellchat.gname <- filterCommunication(cellchat.gname, min.cells = 10) # 用默认值就行

df.net <- subsetCommunication(cellchat.gname) # returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
write.table(df.net, file = "net_LR.txt", quote = F, sep = '\t', col.names = T,row.names = F)

cellchat.gname <- computeCommunProbPathway(cellchat.gname) # 计算通路信号强度，用的默认trimean方法

df.netp <- subsetCommunication(cellchat.gname, slot.name = "netP")
write.table(df.netp, file = "net_pathway.txt", quote = F, sep = '\t', col.names = T,row.names = F)


cellchat.gname <- aggregateNet(cellchat.gname)# 计算通讯网络，可以自定义source和target
cellchat.gname <- netAnalysis_computeCentrality(cellchat.gname, slot.name = "netP") # 计算四种角色的值
groupSize <- as.numeric(table(cellchat.gname@idents)) # 数字向量，用于后续作图

# df.net <- subsetCommunication(cellchat.gname) # 提取表格
setwd(CellTalk_dir)
# write.table(df.net,file = paste0("cellchat.",gname,".xls"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE) # 保存整体信息表格
saveRDS(cellchat.gname, file = paste0("cellchat.",gname,".rds"))# 保存数据用于后续比较分析

# heatmap
gg1 <- netVisual_heatmap(cellchat.gname)
gg2 <- netVisual_heatmap(cellchat.gname, measure = "weight")
pdf(file = paste("netVisual_heatmap", "pdf", sep = "."), width = 10, height = 5)
print(gg1 + gg2)
dev.off()
png(file = paste("netVisual_heatmap", ".png", sep = ""), width = 250, height = 125, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
print(gg1 + gg2)
dev.off()

# 圈图，整体展示细胞间通讯情况，这个应该是单纯的配体-受体对信号强度的加和
pdf(file = paste("cell_cell_commu_network", "pdf", sep = "."), width = 8, height = 6)
netVisual_circle(cellchat.gname@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat.gname@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
dev.off()

png(file = paste("cell_cell_commu_network-1", ".png", sep = ""), width = 200, height = 150, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
netVisual_circle(cellchat.gname@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
dev.off()
# 互作数量
png(file = paste("cell_cell_commu_network-2", ".png", sep = ""), width = 200, height = 150, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
netVisual_circle(cellchat.gname@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
dev.off()
# 互作强度

# 每个细胞类型的信号强度图
mat <- cellchat.gname@net$weight
pdf(file = paste("network_strength_from_each_celltype", "pdf", sep = "."), width = 8, height = 8)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
png(file = paste("network_strength_from_each_celltype", ".png", sep = ""), width = 200, height = 200, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# 每个细胞类型的信号数量图感觉也可以试一下
mat <- cellchat.gname@net$count
pdf(file = paste("network_count_from_each_celltype", "pdf", sep = "."), width = 8, height = 8)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
png(file = paste("network_count_from_each_celltype", ".png", sep = ""), width = 200, height = 200, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

## 每个信号通路的可视化图
tmpdir2 <- paste("", "all_signaling_pathways", sep = "")
if (!file.exists(tmpdir2)) {
  dir.create(tmpdir2)
}
setwd(tmpdir2)

pathways.show.all <- cellchat.gname@netP$pathways # 通路中的所有信号通路
levels(cellchat.gname@idents)
# vertex.receiver <- seq(1, 4) 
# 感觉这种不可控变量干脆去掉算了，有需求自己读rds进来手动画
for (i in 1:length(pathways.show.all)) {
#  netVisual(cellchat.gname, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy", out.format = c("pdf"))
#  netVisual(cellchat.gname, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, out.format = c("pdf"), layout = "circle")
  m <- subset.data.frame(df.net,df.net$pathway_name == pathways.show.all[i])
  num1 <- length(unique(m$receptor)) +length(unique(m$ligand))
  num2 <- length(unique(m$interaction_name))

  pdf(file = paste0(pathways.show.all[i], "_heatmap.pdf"), width = 4, height = 4)
  ht <- list()
  par(mfrow=c(1,1))
  ht[[1]] <- netVisual_heatmap(cellchat.gname, signaling = pathways.show.all[i], color.heatmap = "Reds")
  ComplexHeatmap::draw(ht[[1]])
  dev.off()

  gg <- netAnalysis_contribution(cellchat.gname, signaling = pathways.show.all[i])
  ggsave(filename = paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot = gg, width = 400, height = 30*num2, units = "mm", dpi = 300,limitsize = FALSE)

  gg <- plotGeneExpression(cellchat.gname, signaling = pathways.show.all[i])
  ggsave(filename = paste0(pathways.show.all[i], "_GeneExpression.pdf"), plot = gg, width = 400, height = 30*num1, units = "mm", dpi = 300,limitsize = FALSE)
}
setwd("../")
# 这个循环可以考虑再继续加图




## 气泡图 这里主要是长宽设置不好自适应，看后续怎么调整吧，指定source和target倒是可以，或者干脆不展示
gg <- netVisual_bubble(cellchat.gname, remove.isolate = TRUE)
ggsave(filename = paste0("", "netVisual_bubble.pdf"), plot = gg, width = 500, height = 800, units = "mm", dpi = 300, limitsize = FALSE)
# 热图，展示四种角色的关键程度
pdf(file = paste0("", "signalingRole.pdf"), width = 10, height = 6)
for (i in 1:length(pathways.show.all)) {
  netAnalysis_signalingRole_network(cellchat.gname, signaling = pathways.show.all[i], width = 10, height = 2.5, font.size = 10)
}
dev.off()
# 气泡图，展示各个细胞的输入强度和输出强度，大小代表数目
pdf(file = paste0("", "signalingRole_scatter.pdf"), width = 8, height = 8)
netAnalysis_signalingRole_scatter(cellchat.gname)
dev.off()
# 热图，展示输入输出的细胞类型通路强度
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.gname, pattern = "outgoing", height = 22)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.gname, pattern = "incoming", height = 22)
pdf(file = paste0("", "netAnalysis_signalingRole_heatmap.pdf"), width = 10, height = 12)
print(ht1 + ht2)
dev.off()

# 出图判断可能需要做的冲积图nPattern是多少
pdf(file = paste0("", "outgoing_selectK.pdf"), width = 6, height = 4)
selectK(cellchat.gname, pattern = "outgoing")
dev.off()
#下面这部分图感觉可以放在别的脚本里当做自定义来做
#nPatterns <- 3 # 这个数值最好是根据上面的selectK图来判断，所以新版可以去掉，放在另一个脚本里面
#pdf(file = paste0("", "outgoing_heatmap.pdf"), width = 10, height = 20)
#cellchat.gname <- identifyCommunicationPatterns(cellchat.gname, pattern = "outgoing", k = nPatterns, height = 22)
#dev.off()

#pdf(file = paste0("", "outgoing_river.pdf"), width = 10, height = 12)
#netAnalysis_river(cellchat.gname, pattern = "outgoing")
#dev.off()

#pdf(file = paste0("", "outgoing_dot.pdf"), width = 12, height = 12)
#netAnalysis_dot(cellchat.gname, pattern = "outgoing")
#dev.off()
# 这三张图展示的结果是类似的冲积图感觉最容易看懂

#pdf(file = paste0("", "incoming_selectK.pdf"), width = 6, height = 4)
#selectK(cellchat.gname, pattern = "incoming")
#dev.off()

#nPatterns <- 6
#pdf(file = paste0("", "incoming_heatmap.pdf"), width = 10, height = 20)
#cellchat.gname <- identifyCommunicationPatterns(cellchat.gname, pattern = "incoming", k = nPatterns, height = 22)
#dev.off()

#pdf(file = paste0("", "incoming_river.pdf"), width = 10, height = 12)
#netAnalysis_river(cellchat.gname, pattern = "incoming")
#dev.off()

#pdf(file = paste0("", "incoming_dot.pdf"), width = 12, height = 12)
#netAnalysis_dot(cellchat.gname, pattern = "incoming")
#dev.off()
