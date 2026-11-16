rm(list = ls())
#传参
library(optparse)
option_list <- list(
  make_option(c("-r", "--rdsfile"),action = "store",type = "character",
              help = "The Input rdsfile"),
  make_option(c("-d", "--defaultassay"),action = "store",type = "character",
              default = "RNA",help = "The defaultassay"),
  make_option(c("-s", "--splitby"),action = "store",type = "character",
              default = "group",help = "The splitby"),
  make_option(c("-i", "--ident"),action = "store",type = "character",
              default = "seurat_clusters",help = "The ident"),
  make_option(c("-g", "--genelist"),action = "store",type = "character",
              default = "",help = "The genelist"),
  make_option(c("-m", "--markerlist"),action = "store",type = "character",
              default = "",help = "The markerlist"),
  make_option(c("-o", "--outdir"),action = "store",type = "character",
	      default = "./",help = "The outdir,the default is './'")
)

opt <- parse_args(OptionParser(option_list = option_list)) ####解释器

#加载包
library(tidyverse)
library(Seurat)
library(patchwork)
library(future.apply)
library(RColorBrewer)
library(circlize)
library(stringr)
library(gtools)

# test
#setwd("D:/Workspace/R_learning/sc_analysis/common_plot/")
#opt$rdsfile <- "cells_analysis.rds"
#opt$markerlist <- "top3_combined_markers_export.xls"
#opt$ident <- "customclassif"

#设置并行核心数
plan("multicore", workers = 4)
options(future.globals.maxSize = 60 * 1024^3)

#配色
mycol <- colorRampPalette(brewer.pal(9, "Set1"))(30) #最后的数字可以任意修改
mycol2 <- colorRampPalette(c('blue','red'))(10)

#创建文件夹
if(!file.exists(opt$outdir)){
 dir.create(opt$outdir)
}
setwd(opt$outdir)

#读入
sc <- readRDS(opt$rdsfile)
if(opt$genelist != ""){
  genes <- read.table(opt$genelist,header = F,sep = "\t",
                      check.names = F,com = "",quote = "")
  genes_to_show <- genes[,1] %>% unique()
}else if(opt$markerlist != ""){
  genes <- read.table(opt$markerlist,header = T,sep = "\t",
                      check.names = F,com = "",quote = "")
  genes_to_show <- genes[,2] %>% unique()
}

#基础设置
DefaultAssay(sc) <- opt$defaultassay
Idents(sc) <- opt$ident
#levels(Idents(sc)) <- sort(levels(Idents(sc)))
Idents(sc) <- factor(Idents(sc),levels=mixedsort(levels(Idents(sc))))
unit <- 4 #作图基本长度单位
identnum <- length(levels(Idents(sc)))
splitnum <- length(unique(sc@meta.data[,opt$splitby]))

#设定图形宽度
mywidth <- identnum
if(mywidth <= 2){
  mywidth <- 2 * mywidth
}else if(mywidth <= 5){
  mywidth <- round(1.5 * mywidth)
}else if(mywidth <= 10){
  mywidth <- 1 * mywidth
}else if(mywidth <= 15){
  mywidth <- round(mywidth / 1.5)
}else{
  mywidth <- round(mywidth / 2)
}

#作图函数
#使用tryCatch捕获异常，即使出现error也会继续执行
#VlnPlot
vlnplot <- function(gene,splitby){
  tryCatch(
    {
      #split by 
      p1 <- VlnPlot(sc, features = gene, split.by = splitby, pt.size = 0)
      ggsave(paste0("VlnPlot_split_",gene,".png"),p1,width = mywidth,height = unit)
      ggsave(paste0("VlnPlot_split_",gene,".pdf"),p1,width = mywidth,height = unit)
      #not split by
      p1 <- VlnPlot(sc, features = gene, pt.size = 0)+NoLegend()
      ggsave(paste0("VlnPlot_",gene,".png"),p1,width = mywidth,height = unit)
      ggsave(paste0("VlnPlot_",gene,".pdf"),p1,width = mywidth,height = unit)      
    },
    # warning = function(w){message("warning")},
    error = function(e){message("基因不存在")}
  )
}

#FeaturePlot umap
featureplot_umap <- function(gene,splitby,splitnum){
  tryCatch(
    {
      #分面
      p1 <- FeaturePlot(sc, reduction = "umap", features = gene, split.by = splitby, cols = c("grey", "red"))
      ggsave(paste0("FeaturePlot_umap_split_",gene,".png"),p1,width = (unit+1)*splitnum,height = unit)
      ggsave(paste0("FeaturePlot_umap_split_",gene,".pdf"),p1,width = (unit+1)*splitnum,height = unit)
      #不分面
      p1 <- FeaturePlot(sc, reduction = "umap", features = gene, cols = c("grey", "red"))
      ggsave(paste0("FeaturePlot_umap_",gene,".png"),p1,width = unit+1,height = unit)
      ggsave(paste0("FeaturePlot_umap_",gene,".pdf"),p1,width = unit+1,height = unit)
    },
    # warning = function(w){message("warning")},
    error = function(e){message("基因不存在")}
  )
}

#DotPlot
dotplot <- function(gene){
  #使用seurat_clusters
  p = DotPlot(sc, features = gene,group.by = "seurat_clusters") + RotatedAxis()
  w = length(gene)/2
  h = length(unique(sc$seurat_clusters))/2
  ggsave('DotPlot_bycluster.png',p,width = w,height = h)
  ggsave('DotPlot_bycluster.pdf',p,width = w,height = h)
  
  #使用customclassif
  p = DotPlot(sc, features = gene,group.by = "customclassif") + RotatedAxis()
  w = length(gene)/2
  h = length(unique(sc$customclassif))/2
  ggsave('DotPlot_bycelltype.png',p,width = w,height = h)
  ggsave('DotPlot_bycelltype.pdf',p,width = w,height = h)
  
  # #使用customclassif和group
  # sc$celltype_group <- paste0(sc$customclassif,"_",sc$group)
  # p = DotPlot(sc, features = gene,group.by = "celltype_group") + RotatedAxis()
  # w = length(gene)/2
  # h = length(unique(sc$customclassif))/2
  # ggsave('DotPlot_bycelltype_bygroup.png',p,width = w * 1.5,height = h * 2.5)
  # ggsave('DotPlot_bycelltype_bygroup.pdf',p,width = w * 1.5,height = h * 2.5)
  
} 

#DoHeatmap
doheatmap <- function(gene){

  scRNA <- sc
  scRNA <- ScaleData(scRNA, features = rownames(scRNA)) #防止没有scale 

  #使用seurat_clusters
  p = DoHeatmap(scRNA,features = gene,group.by = "seurat_clusters")+NoLegend()
  w = length(gene)/2
  h = length(unique(scRNA$seurat_clusters))/2
  ggsave("DoHeatmap_bycluster.png",p,width = w,height = h)
  ggsave("DoHeatmap_bycluster.pdf",p,width = w,height = h)
  
  #使用customclassif
  p = DoHeatmap(scRNA,features = gene,group.by = "customclassif")+NoLegend()
  w = length(gene)/2
  h = length(unique(scRNA$customclassif))/2
  ggsave("DoHeatmap_bycelltype.png",p,width = w,height = h)
  ggsave("DoHeatmap_bycelltype.pdf",p,width = w,height = h)
  
} 

#根据需要运行相应的作图函数
#并行作图,mapply是多参数的sapply
future_mapply(vlnplot,genes_to_show,opt$splitby,future.seed=TRUE)
future_mapply(featureplot_umap,genes_to_show,opt$splitby,splitnum,future.seed=TRUE)
#dotplot(genes_to_show)
#doheatmap(genes_to_show)

