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
# https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
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

# The Bioconductor package SingleCellExperiment provides the SingleCellExperiment classfor usage.While the package is implicitly installed and loaded when using any package that depends on the SingleCellExperiment class, it can be explicitly installed(and loaded) as follows: #install.packages("BiocManager")# BiocManager::install(c('SingleCellExperiment', 'Rtsne', 'scater', 'scran', 'uwot'))
# library(SingleCellExperiment)
# library(scater)

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
    default = 4,
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
    c("-a", "--postname"),
    action = "store",
    type = "character",
    help = "post sample name, post vs pre"
  ),
  make_option(
    c("-b", "--prename"),
    action = "store",
    type = "character",
    help = "pre sample name, post vs pre"
  ),
  make_option(
    c("-s", "--source"),
    action = "store",
    type = "character",
    help = "selected cluster: 1,2,3,4"
  ),
  make_option(
    c("-t", "--target"),
    action = "store",
    type = "character",
    help = "target cluster: 5,6,7,8"
  ),
  make_option(
    c("-o", "--outdir"),
    action = "store",
    type = "character",
    help = "outdir",
    default = "./"
  ),
  make_option("--print",action = "store_true",default = FALSE,
    help = "输出cluster"),
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

postname <- opt$postname
prename <- opt$prename

CellTalk_dir <- opt$outdir
if (!file.exists(CellTalk_dir)) {
  dir.create(CellTalk_dir)
}
setwd(CellTalk_dir)
cellchat.post <- readRDS(file = paste("sample_",postname,"/cellchat.",postname, ".rds", sep = ""))
cellchat.pre <- readRDS(file = paste("sample_",prename,"/cellchat.",prename, ".rds", sep = ""))

##########################
# 两组样本比较 opt$postname vs opt$prename
##########################

setwd(CellTalk_dir)
pos.dataset <- postname
pre.dataset <- prename
cmpdir <- paste0("compare-", pos.dataset, "_vs_", pre.dataset)
if (!file.exists(cmpdir)) {
  dir.create(cmpdir)
}
setwd(cmpdir)

#分两组取并集
a=levels(cellchat.pre@idents)
b=levels(cellchat.post@idents)
group.new = sort(unique(c(a,b)))

cellchat.pre = liftCellChat(cellchat.pre, group.new)
cellchat.post = liftCellChat(cellchat.post, group.new)

#将list合并成cellchat对象
object.list <- list(pre = cellchat.pre, post = cellchat.post)
names(object.list) <- c(pre.dataset,pos.dataset)
cellchat <- mergeCellChat(object.list, add.names = c(pre.dataset,pos.dataset))

##这里只保留在两组样本中同时存在的细胞类型
group.new = intersect(a,b)
cellchat <- subsetCellChat(cellchat, idents.use = group.new)

#
cat("共同的细胞类型：",group.new,"\n")

# start
selected_cluster = c(1:length(group.new))
chat_cluster = c(1:length(group.new))

# 用热图展示两组共有的信号通路 输入 输出 输入输出加和 的信号强度
i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 22)
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i + 1], width = 10, height = 22)
pdf(file = paste("netAnalysis_signalingRole_heatmap_outgoing", "pdf", sep = "."), width = 10, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
png(file = paste("netAnalysis_signalingRole_heatmap_outgoing", ".png", sep = ""), width = 250, height = 250, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 22, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i + 1], width = 10, height = 22, color.heatmap = "GnBu")
pdf(file = paste("netAnalysis_signalingRole_heatmap_incoming", "pdf", sep = "."), width = 10, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
png(file = paste("netAnalysis_signalingRole_heatmap_incoming", ".png", sep = ""), width = 250, height = 250, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 25, color.heatmap = "OrRd")
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i + 1], width = 10, height = 25, color.heatmap = "OrRd")
pdf(file = paste("netAnalysis_signalingRole_heatmap_all", "pdf", sep = "."), width = 10, height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
png(file = paste("netAnalysis_signalingRole_heatmap_all", ".png", sep = ""), width = 250, height = 375, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

# 用柱状图展示两组的互作强度和互作数量
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), measure = "weight")
pdf(file = paste("compareInteractions", "pdf", sep = "."), width = 6, height = 3)
print(gg1 + gg2)
dev.off()
png(file = paste("compareInteractions", ".png", sep = ""), width = 150, height = 75, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
print(gg1 + gg2)
dev.off()

# 圈图展示各细胞输入/输出的强度/数量关系
pdf(file = paste("netVisual_diffInteraction", "pdf", sep = "."), width = 8, height = 6)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()
png(file = paste("netVisual_diffInteraction-1", ".png", sep = ""), width = 200, height = 150, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
netVisual_diffInteraction(cellchat, weight.scale = T)
dev.off()
png(file = paste("netVisual_diffInteraction-2", ".png", sep = ""), width = 200, height = 150, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()
# 可以直接展示特定细胞之间的互作具体数量和强度
# group.cellType <- c(rep("Treg", 2), rep("CD8+TE", 2), rep("CD4+TN", 2)) # 有几组就重复几次
# group.cellType <- factor(group.cellType, levels = c("Treg", "CD8+TE", "CD4+TN")) # 因子转化
# object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)}) 
# cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
# par(mfrow = c(1,2), xpd=TRUE)
# for (i in 1:length(object.list)) {
#  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
#}
# 比较适合自定的时候做

# 热图展示不同细胞之间的通讯强度（定量），左边是数量，右边是强度
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
pdf(file = paste("netVisual_heatmap", "pdf", sep = "."), width = 10, height = 5)
print(gg1 + gg2)
dev.off()
png(file = paste("netVisual_heatmap", ".png", sep = ""), width = 250, height = 125, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
print(gg1 + gg2)
dev.off()

# 这里的结果是想要做上方备注是的那部分内容，但是并没有定量，可能是因为版本更新的原因，也可能是因为没有确定特定细胞类型的原因

weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))

pdf(file = paste("netVisual_circle", "pdf", sep = "."), width = 10, height = 5)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge = F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()
png(file = paste("netVisual_circle", ".png", sep = ""), width = 250, height = 125, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge = F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

# 常规点图，展示各细胞类型输入输出强度，大小和数量成正比

num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
pdf(file = paste("netAnalysis_signalingRole_scatter", "pdf", sep = "."), width = 10, height = 5)
patchwork::wrap_plots(plots = gg)
dev.off()
png(file = paste("netAnalysis_signalingRole_scatter", ".png", sep = ""), width = 250, height = 125, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
patchwork::wrap_plots(plots = gg)
dev.off()

# 柱状图展示两组之间不同通路的权重分配，红色颜色字体的通路代表在A组富集，蓝色颜色字体的通路代表在B组富集

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
pdf(file = paste("rankNet", "pdf", sep = "."), width = 10, height = 10)
print(gg1 + gg2)
dev.off()
png(file = paste("rankNet", ".png", sep = ""), width = 250, height = 250, units = "mm", res = 300, pointsize = 10)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
print(gg1 + gg2)
dev.off()


# 气泡图展示不同细胞间通路的配体-受体的信号强度 两个组
pdf(file = paste("netVisual_bubble", "pdf", sep = "."), width = 8, height = 10)
netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45)
dev.off()

# pdf(file = paste("netVisual_bubble_in_selected_cell", "pdf", sep = "."), width = 8, height = 10)
# netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45)
# dev.off()

# 气泡图展示差异 往往因为图片太大不具备可展示性，倾向于确定了目标细胞类型互作之后再自定义展示，
gg1 <- netVisual_bubble(cellchat, sources.use = selected_cluster, targets.use = chat_cluster, comparison = c(1, 2), max.dataset = 2, title.name = paste0("Increased signaling in ", names(object.list)[2]), angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = selected_cluster, targets.use = chat_cluster, comparison = c(1, 2), max.dataset = 1, title.name = paste0("Decreased signaling in ", names(object.list)[2]), angle.x = 45, remove.isolate = T)

# pdf(file = paste("netVisual_bubble_updown_out_selected_cell", "pdf", sep = "."), width = 30, height = 15)
# print(gg1 + gg2)
# dev.off()

# 输出表格的信息是从图形数据中抓取出来的
cat(c(''),file=paste(postname,'_vs_',prename,'.up.prob','.xls',sep=''))
write.table(gg1$data,file=paste(postname,'_vs_',prename,'.up.prob','.xls',sep=''),append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

cat(c(''),file=paste(postname,'_vs_',prename,'.down.prob','.xls',sep=''))
write.table(gg2$data,file=paste(postname,'_vs_',prename,'.down.prob','.xls',sep=''),append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

# 气泡图展示差异 往往因为图片太大不具备可展示性，倾向于确定了目标细胞类型互作之后再自定义展示
gg1 <- netVisual_bubble(cellchat, sources.use = chat_cluster, targets.use = selected_cluster, comparison = c(1, 2), max.dataset = 2, title.name = paste0("Increased signaling in ", names(object.list)[2]), angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = chat_cluster, targets.use = selected_cluster, comparison = c(1, 2), max.dataset = 1, title.name = paste0("Decreased signaling in ", names(object.list)[2]), angle.x = 45, remove.isolate = T)

#pdf(file = paste("netVisual_bubble_updown_in_selected_cell", "pdf", sep = "."), width = 30, height = 15)
#print(gg1 + gg2)
#dev.off()
# 重复表格不用输出
# cat(c(''),file=paste(postname,'_vs_',prename,'.up.probability.sub.in','.xls',sep=''))
# write.table(gg1$data,file=paste(postname,'_vs_',prename,'.up.probability.sub.in','.xls',sep=''),append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

# cat(c(''),file=paste(postname,'_vs_',prename,'.down.probability.sub.in','.xls',sep=''))
# write.table(gg2$data,file=paste(postname,'_vs_',prename,'.down.probability.sub.in','.xls',sep=''),append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

# 计算logfc和pvalue，提取一个总表
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "cell_group", pos.dataset = pos.dataset, features.name = pos.dataset, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = pos.dataset)
# cat(c(''),file=paste(postname,'_vs_',prename,'_net_deg','.xls',sep=''))
write.table(net,file=paste(postname,'_vs_',prename,'_net_raw','.xls',sep=''),append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)
# 设定阈值，挑选上调通路 
net.up <- subsetCommunication(cellchat, net = net, datasets = pos.dataset,ligand.logFC = 0.1, receptor.logFC = 0.1)
cat(c(''),file=paste(postname,'_vs_',prename,'_net.up','.xls',sep=''))
write.table(net.up,file=paste(postname,'_vs_',prename,'_net.up','.xls',sep=''),append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)
# 设定阈值，挑选下调通路
net.down <- subsetCommunication(cellchat, net = net, datasets = pre.dataset,ligand.logFC = -0.1, receptor.logFC = -0.1)
cat(c(''),file=paste(postname,'_vs_',prename,'_net.down','.xls',sep=''))
write.table(net.down,file=paste(postname,'_vs_',prename,'_net.down','.xls',sep=''),append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)
# 这里筛选的阈值官方文档中并没有给出建议，理论上来说只要>0就是上调，<0就是下调，且受体FC和配体FC是分开的，可以进行分开挑选
# 提取上调和下调的配体 受体symbol
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# 和上面的表是一致的，没有必要
# net.up.sub <- subsetCommunication(cellchat, net = net, datasets = pos.dataset,ligand.logFC = 0.1, receptor.logFC = 0.1, sources.use = selected_cluster, targets.use = chat_cluster)
# cat(c(''),file=paste(postname,'_vs_',prename,'_net.up.sub.out','.xls',sep=''))
# write.table(net.up.sub,file=paste(postname,'_vs_',prename,'_net.up.sub.out','.xls',sep=''),append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

pairLR.use.up <- net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, comparison = c(1, 2), angle.x = 90, remove.isolate = T, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

# net.down.sub <- subsetCommunication(cellchat, net = net, datasets = pre.dataset,ligand.logFC = -0.1, receptor.logFC = -0.1, sources.use = selected_cluster, targets.use = chat_cluster)
# cat(c(''),file=paste(postname,'_vs_',prename,'_net.down.sub.out','.xls',sep=''))
# write.table(net.down.sub,file=paste(postname,'_vs_',prename,'_net.down.sub.out','.xls',sep=''),append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

pairLR.use.down <- net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, comparison = c(1, 2), angle.x = 90, remove.isolate = T, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
# 通过气泡图展示差异配体-受体,数量太多，不展示
#pdf(file = paste("netVisual_bubble_deg", "pdf", sep = "."), width = 20, height = 15)
#print(gg1 + gg2)
#dev.off()
#png(file = paste("netVisual_bubble_deg", ".png", sep = ""), width = 500, height = 375, units = "mm", res = 300, pointsize = 10)
#par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
#print(gg1 + gg2)
#dev.off()

# net.up.sub <- subsetCommunication(cellchat, net = net, datasets = pos.dataset,ligand.logFC = 0.1, receptor.logFC = 0.1, sources.use = chat_cluster, targets.use = selected_cluster)
# cat(c(''),file=paste(postname,'_vs_',prename,'_net.up.sub.in','.xls',sep=''))
# write.table(net.up.sub,file=paste(postname,'_vs_',prename,'_net.up.sub.in','.xls',sep=''),append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

# pairLR.use.up <- net.up.sub[, "interaction_name", drop = F]
# gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = chat_cluster, targets.use = selected_cluster, comparison = c(1, 2), angle.x = 90, remove.isolate = T, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
# net.down.sub <- subsetCommunication(cellchat, net = net, datasets = pre.dataset,ligand.logFC = -0.1, receptor.logFC = -0.1, sources.use = chat_cluster, targets.use = selected_cluster)
# cat(c(''),file=paste(postname,'_vs_',prename,'_net.down.sub.in','.xls',sep=''))
# write.table(net.down.sub,file=paste(postname,'_vs_',prename,'_net.down.sub.in','.xls',sep=''),append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

#pairLR.use.down <- net.down.sub[, "interaction_name", drop = F]
#gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = chat_cluster, targets.use = selected_cluster, comparison = c(1, 2), angle.x = 90, remove.isolate = T, title.name = paste0("Down-regulated signaling in ", names(object.list)[2])) # nolint
#pdf(file = paste("netVisual_bubble_deg_in_selected_cell", "pdf", sep = "."), width = 20, height = 15)
#print(gg1 + gg2)
#dev.off()
#png(file = paste("netVisual_bubble_deg_in_selected_cell", ".png", sep = ""), width = 500, height = 375, units = "mm", res = 300, pointsize = 10)
#par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
#print(gg1 + gg2)
#dev.off()
# 重复代码，没有意义

# 使用弦图可视化上调和下调信号配体-受体对,经常会因为矩阵太大报错退出，先取消了
net.up <- subsetCommunication(cellchat, net = net, datasets = pos.dataset,ligand.logFC = 0.1, receptor.logFC = 0.1, sources.use = selected_cluster, targets.use = chat_cluster)
top50_ligand_up <- net.up %>% arrange(desc(ligand.logFC)) %>% slice_head(n = 50)
top50_receptor_up <- net.up %>% arrange(desc(receptor.logFC)) %>% slice_head(n = 50)
net.up.select <- bind_rows(top50_ligand_up, top50_receptor_up) %>% distinct()

net.down <- subsetCommunication(cellchat, net = net, datasets = pre.dataset,ligand.logFC = -0.1, receptor.logFC = -0.1, sources.use = selected_cluster, targets.use = chat_cluster)
top50_ligand_down <- net.down %>% arrange(ligand.logFC) %>% slice_head(n = 50)
top50_receptor_down <- net.down %>% arrange(receptor.logFC) %>% slice_head(n = 50)
net.down.select <- bind_rows(top50_ligand_down, top50_receptor_down) %>% distinct()

pdf(file = paste("netVisual_chord_up_deg", "pdf", sep = "."), width = 10, height = 10)
netVisual_chord_gene(object.list[[2]], slot.name = "net", net = net.up.select, lab.cex = 0.8, small.gap = 1, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()
pdf(file = paste("netVisual_chord_down_deg", "pdf", sep = "."), width = 10, height = 10)
netVisual_chord_gene(object.list[[1]], slot.name = "net", net = net.down.select, lab.cex = 0.8, small.gap = 1, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()
#net.up <- subsetCommunication(cellchat, net = net, datasets = pos.dataset,ligand.logFC = 0.1, receptor.logFC = 0.1, sources.use = chat_cluster, targets.use = selected_cluster)
#net.down <- subsetCommunication(cellchat, net = net, datasets = pre.dataset,ligand.logFC = -0.1, receptor.logFC = -0.1, sources.use = chat_cluster, targets.use = selected_cluster)

#pdf(file = paste("netVisual_chord_gene_deg", "pdf", sep = "."), width = 10, height = 10)
#netVisual_chord_gene(object.list[[2]], sources.use = chat_cluster, targets.use = selected_cluster, slot.name = "net", net = net.up.sub, lab.cex = 0.8, small.gap = 1, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#netVisual_chord_gene(object.list[[1]], sources.use = chat_cluster, targets.use = selected_cluster, slot.name = "net", net = net.down.sub, lab.cex = 0.8, small.gap = 1, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#dev.off()



# 和弦图展示两组的通路互作关系 一般不展示
#pdf(file = paste("netVisual_chord_netP", "pdf", sep = "."), width = 20, height = 20)
# par(mfrow = c(1, 2), xpd=TRUE)
#for (i in 1:length(object.list)) {
#  netVisual_chord_gene(object.list[[i]],slot.name = "netP", lab.cex = 0.8, small.gap = 1, title.name = paste0("signaling received by selected cell - ", names(object.list)[i]))
#} 
#dev.off()

#png(file = paste("netVisual_chord_netP", "png", sep = "."), width = 1000, height = 500, units = "mm", res = 300, pointsize = 10)
#par(mfrow = c(1, 2), xpd=TRUE)
#for (i in 1:length(object.list)) {
#  netVisual_chord_gene(object.list[[i]], sources.use = chat_cluster, targets.use = selected_cluster,slot.name = "netP", lab.cex = 0.8, small.gap = 1, title.name = paste0("signaling received by selected cell - ", names(object.list)[i]))
#} 
#dev.off()

# png(file = paste("netVisual_chord_gene_from_selected_cell_netP", "png", sep = "."), width = 1000, height = 500, units = "mm", res = 300, pointsize = 10)
# par(mfrow = c(1, 2), xpd=TRUE)
# for (i in 1:length(object.list)) {
#  netVisual_chord_gene(object.list[[i]], sources.use = selected_cluster, targets.use = chat_cluster,slot.name = "netP", lab.cex = 0.8, small.gap = 1, title.name = paste0("signaling sending from selected cell - ", names(object.list)[i]))
#} 
#dev.off()
#pdf(file = paste("netVisual_chord_gene_from_selected_cell_netP", "pdf", sep = "."), width = 20, height = 20)
# par(mfrow = c(1, 2), xpd=TRUE)
#for (i in 1:length(object.list)) {
#  netVisual_chord_gene(object.list[[i]], sources.use = selected_cluster, targets.use = chat_cluster,slot.name = "netP", lab.cex = 0.8, small.gap = 1, title.name = paste0("signaling sending from selected cell - ", names(object.list)[i]))
#} 
#dev.off()

## 每个信号通路的可视化图
tmpdir2 <- paste("", "all_signaling_pathways", sep = "")
if (!file.exists(tmpdir2)) {
  dir.create(tmpdir2)
}
setwd(tmpdir2)


pathways.show.all <- object.list[[1]]@netP$pathways
pathways.show.all2 <- object.list[[2]]@netP$pathways
# 两组共同的通路展示
# levels(object.list[[1]]@idents)
# vertex.receiver <- seq(1, 4)
for (i in 1:length(pathways.show.all)) {
  if (pathways.show.all[i] %ni% pathways.show.all2) {
    next
  }
  pdf(file = paste0(pathways.show.all[i], "_circle.pdf"), width = 8, height = 8)
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show.all[i]) # control the edge weights across different datasets
  par(mfrow = c(1, 2), xpd = TRUE)
  for (j in 1:length(object.list)) {
    netVisual_aggregate(object.list[[j]], signaling = pathways.show.all[i], layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show.all[i], names(object.list)[j]))
  }
  dev.off()

  pdf(file = paste0(pathways.show.all[i], "_heatmap.pdf"), width = 8, height = 8)
  par(mfrow = c(1, 2), xpd = TRUE)
  ht <- list()
  for (j in 1:length(object.list)) {
    ht[[j]] <- netVisual_heatmap(object.list[[j]], signaling = pathways.show.all[i], color.heatmap = "Reds", title.name = paste(pathways.show.all[i], "signaling ", names(object.list)[j]))
  }
  ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
  dev.off()

  cellchat@meta$cell_group <- factor(cellchat@meta$cell_group, levels = c(pre.dataset, pos.dataset)) # set factor level
  gg <- plotGeneExpression(cellchat, signaling = pathways.show.all[i], split.by = "cell_group", colors.ggplot = T)
  ggsave(filename = paste0(pathways.show.all[i], "_GeneExpression.pdf"), plot = gg, width = 500, height = 500, units = "mm", dpi = 300)
}


setwd("../")

# end


# Save the merged CellChat object
saveRDS(cellchat, file = paste0("cellchat_", pos.dataset, "_vs_", pre.dataset, ".rds"))

