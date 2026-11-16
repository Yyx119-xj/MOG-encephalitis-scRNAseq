#!/usr/bin/env Rscript

###Set VERSION
VERSION = "3.0.0"

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1 && (args[1] == "-v" | args[1] == "--version")){
  message("Version: \n\t", VERSION, "\n")
  quit(status = 0)
}

####################################################################################
### Copyright (C) 2019-2022 by ABLIFE
####################################################################################
# 名称：single cell pipeline
# 描述：单细胞分析流程
# 作者：CHAO CHENG
# 创建时间：2020-10-15
# 联系方式：chaocheng@ablife.cc
####################################################################################
### 修改记录
####################################################################################
# Date           Version       Author            ChangeLog
# 2019-12-27      v1.0          chaocheng         修改测试版本
# 2020-10-15      v2.0          chaocheng         整合seurat
# 2021-3-8        v3.0          chaocheng         流程整合
#####################################################################################
#####################################################################################
#####参考说明
# https://osca.bioconductor.org/overview.html#quick-start
# https://satijalab.org/seurat/v3.2/mca.html
# https://satijalab.org/seurat/sc.combined3k_tutorial.html
# https://github.com/CostaLab/scrna_seurat_pipeline
#####################################################################################


#####################################################################################
#####参数获取
#####################################################################################
# windows系统检查当前程序路径
# script.dir <- dirname(sys.frame(1)$ofile)

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(optparse))      ## Options
suppressPackageStartupMessages(library(futile.logger)) ## logger
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(doParallel))
library(SingleCellExperiment)
library(scater)
library(tidyverse)
library(grDevices)
library(gplots)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(metap)
library(clustree)
library(presto)
library(gtools)
library(future.apply)
`%ni%` <- Negate(`%in%`)

## 参数读取
option_list <- list(
    make_option(
        c("-n","--ncpus"),
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
        help = "sc processing rds file"
    ),
    make_option(
        c("-s","--samplelist"),
        action = "store",
        type = "character",
        default = "",
        help = "samplename1,samplename2,samplename3"
    ),
    make_option(
        c("-l","--listname"),
        action = "store",
        type = "character",
        default = "",
        help = "listname1,listname2,listname3"
    ),
    make_option(
        c("-d", "--dims"),
        action = "store",
        type = "integer",
        default = 30,
        help = "PCA dims number"
    ),
    make_option(
        c("--resolution"),
        action = "store",
        type = "double",
        default = 0.4,
        help = "cluster resolution"
    ),    
    make_option(
        c("--sctypedb"),
        action = "store",
        type = "character",
        default = "/public/pipeline/ablifepipe/latest/src/SingleCell/sctype/ScTypeDB_full.xlsx",
        help = "sctypedb"
    ),    
    make_option(
        c("-t","--tissue"),
        action = "store",
        type = "character",
        default = "Immune system",
        help = "Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Embryo,Gastrointestinal tract,Muscle,Skin,Heart,Ovary,Testis,White adipose tissue,Teeth"
    ),
    make_option(
        c("--rdsfile"),
        action = "store",
        type = "character",
        default = "",
        help = "The data object rds file containing clusters"
    ),
    make_option(
        c("--rmcluster"),
        action = "store",
        type = "character",
        default = "",
        help = "remove clusters from raw seurat data. eg. 1,3,5"
    ),
    # make_option(
    #     c("-k","--keeps"),
    #     action = "store",
    #     type = "character",
    #     default = "",
    #     help = "c1,c2,c3"
    # ),
    # make_option(
    #     c("--defaultassay"),
    #     action = "store",
    #     type = "character",
    #     default = "RNA",
    #     help = "default assay"
    # ),
    make_option(
        c("--skiptsne"),
        action = "store",
        type = "character",
        default = "false",
        help = "skip tsne"
    ),
    make_option(
        c("--skipanno"),
        action = "store",
        type = "character",
        default = "false",
        help = "skip cell anno, true or false"
    ),
    # make_option(
    #     c("--export"),
    #     action = "store",
    #     type = "character",
    #     default = "false",
    #     help = "export count and exp file, skip when it is false"
    # ),
    make_option(
        c("--cellanno"),
        action = "store",
        type = "character",
        default = "",
        help = "The cluster to cell type Annotation file"
    ),
    make_option(
        c("--nmethod"),
        action = "store",
        type = "character",
        default = "log",
        help = "normalize methods, log or sct, default is log"
    ),
    make_option(
        c("--imethod"),
        action = "store",
        type = "character",
        default = "cca",
        help = "integrate methods, cca or harmony, default is cca"
    ),
    make_option(
        c("--mmethod"),
        action = "store",
        type = "character",
        default = "seurat",
        help = "find marker methods, seurat or presto, default is seurat"
    ) 
)
opt <- parse_args(OptionParser(option_list = option_list))

registerDoParallel(cores=opt$ncpus)
plan("multicore", workers = opt$ncpus)
options(future.globals.maxSize = opt$MaxMemMega * 1024^2)

# ## logger记录
# if(!file.exists("logs")){
#   dir.create("logs")
# }
# LOGNAME = "log-"
# cur_date <- as.character(Sys.Date())
# errorLog  =  file.path("logs", sprintf("%sERROR-%s.log", LOGNAME, cur_date))
# warnLog   =  file.path("logs", sprintf("%sWARN-%s.log", LOGNAME, cur_date))
# infoLog   =  file.path("logs",  sprintf("%sINFO-%s.log", LOGNAME, cur_date))
# traceLog   =  file.path("logs",  sprintf("%sTRAC-%s.log", LOGNAME, cur_date))

# invisible(flog.logger("error", ERROR, appender.file(errorLog)))
# invisible(flog.logger("warn", WARN, appender.file(warnLog)))
# invisible(flog.logger("info", INFO, appender.file(infoLog)))
# invisible(flog.logger("trace", TRACE, appender.file(traceLog)))
# invisible(flog.appender(appender.console(), name = "ROOT"))

# logger.info <- function(msg, ...) {
#   flog.info(msg, ..., name = "ROOT")
#   flog.info(msg, ..., name = "info")
# }

# logger.warn <- function(msg, ...) {
#   flog.warn(msg, ..., name = "ROOT")
#   flog.warn(msg, ..., name = "info")
#   flog.warn(msg, ..., name = "warn")
# }

# logger.error <- function(msg, ...) {
#   flog.error(msg, ..., name = "ROOT")
#   flog.error(msg, ..., name = "info")
#   flog.error(msg, ..., name = "warn")
#   flog.error(msg, ..., name = "error")
#   flog.error(msg, ..., name = "trace")
# }

## 主题设置
theme_paper <- theme(
    # panel.border = element_rect(fill = NA,colour = "black"),
    panel.grid.major = element_line(size = 0.1, linetype = "dashed"),
    # panel.grid.minor = element_line(
    #     colour = "grey90",
    #     size = 0.2,
    #     linetype = "dashed"
    # ),
    # axis.text.x= element_text(vjust = 1,hjust = 1, angle = 45, size = 5),
    # legend.position = "top",
    # legend.direction = "horizontal",
    legend.key = element_rect(fill = 'white', color = 'white'),
    legend.key.size = unit(0.5, "cm"),
    legend.position = c(0.2, 0.88),
    # legend.position = "top",
    # legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.background = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0.5, 1, 0.5, 1, unit = "cm"),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, colour = "black")
)

## 测试参数
# We start by reading in the data. The Read10X function reads in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).
# opt$samplelist <- "HC_1,HC_2,LTBI_1,LTBI_2"
# opt$listname <- "HC,HC,LTBI,LTBI"
# opt$infile <- "/users/wuyf/public/Singlecell/pipeline/v2/1_preprocessing_log_cca/sc_preprocessing.rds"
# opt$sctypedb <- "/users/wuyf/public/Singlecell/pipeline/sct_v1_harmony/marker/ScTypeDB_PBMC.xlsx"
# opt$tissue <- "PBMC"
# opt$nmethod <- "log"
# opt$imethod <- "cca"

samplenames <- strsplit(opt$samplelist,",")[[1]]
listnames <- strsplit(opt$listname,",")[[1]]
groupnum = length(unique(listnames))
samplenum = length(unique(samplenames))
res = opt$resolution

#####################################################################################
## preprocessing数据读取，数据对象建立
#####################################################################################

# Load the sc.combined dataset
if(opt$rdsfile == ""){
  sc.combined <- readRDS(file = opt$infile)
}else{
  sc.combined <- readRDS(file = opt$rdsfile)
  Idents(sc.combined) <- "seurat_clusters"
  if(opt$rmcluster != ""){
    rmcluster <- strsplit(opt$rmcluster,",")[[1]]
    keeps <- setdiff(unique(factor(Idents(sc.combined))), rmcluster)
    print("keep clusters:")
    print(keeps)
    print("remove clusters:")
    print(rmcluster)

    sc.combined <- subset(sc.combined, idents= keeps)
  }
}

workdir <- getwd()

#####################################################################################
## 降维聚类
#####################################################################################
sc_cluster_dir = paste(workdir,"/","sc_cluster", sep="")
if(!file.exists(sc_cluster_dir)){
  dir.create(sc_cluster_dir)
}

#降维聚类的函数
clustering <- function(d,res,ass="RNA",red="pca"){
  #参数d 代表 dims
  #参数res 代表 resolution
  #参数ass 代表 assay
  #参数red 代表 reduction

  DefaultAssay(sc.combined) <- ass
  sc.combined <- FindNeighbors(sc.combined, reduction = red, dims = 1:d)
  sc.combined <- FindClusters(sc.combined, resolution = seq(0.3,1.6,by=0.1),seed=TRUE)
  p1 <- clustree(sc.combined)
  pdf(file=paste("resolution_clustree",".pdf",sep=""),width = 20,height = 20)
  print(p1)
  dev.off()

  png(file = paste("resolution_clustree",".png",sep=""),width = 500,height = 500,units = "mm",res = 300,pointsize = 1.5)
  par(mar = c(5, 5, 2, 2), cex.axis = 1.5, cex.lab = 2)
  print(p1)
  dev.off()

  Idents(object = sc.combined) <- "RNA_snn_res.0.4"
  print("resolution:0.4")
  length(unique(factor(Idents(sc.combined))))
  Idents(object = sc.combined) <- "RNA_snn_res.0.6"
  print("resolution:0.6")
  length(unique(factor(Idents(sc.combined))))
  Idents(object = sc.combined) <- "RNA_snn_res.0.8"
  print("resolution:0.8")
  length(unique(factor(Idents(sc.combined))))
  Idents(object = sc.combined) <- "RNA_snn_res.1"
  print("resolution:1")
  length(unique(factor(Idents(sc.combined))))
  Idents(object = sc.combined) <- "RNA_snn_res.1.2"
  print("resolution:1.2")
  length(unique(factor(Idents(sc.combined))))
  Idents(object = sc.combined) <- "RNA_snn_res.1.4"
  print("resolution:1.4")
  length(unique(factor(Idents(sc.combined))))
  Idents(object = sc.combined) <- "RNA_snn_res.1.6"
  print("resolution:1.6")
  length(unique(factor(Idents(sc.combined))))

  sc.combined <- FindClusters(sc.combined, resolution = res)

  #
  # sc.combined <- RunTSNE(sc.combined, dims = 1:d, tsne.method = "FIt-SNE", nthreads = 16, max_iter = 2000, reduction = red)
  sc.combined <- RunTSNE(sc.combined, dims = 1:d, nthreads = 16, max_iter = 2000, reduction = red)
  sc.combined <- RunUMAP(sc.combined, dims = 1:d, reduction = red)
  return(sc.combined)
}

setwd(sc_cluster_dir)
## 降维聚类

if(length(unique(sc.combined$replicate))==1 & opt$nmethod=="log" & opt$rdsfile==""){
  sc.combined <- clustering(opt$dims,opt$resolution,"RNA","pca")
}else if(length(unique(sc.combined$replicate))==1 & opt$nmethod=="sct" & opt$rdsfile==""){
  sc.combined <- clustering(opt$dims,opt$resolution,"SCT","pca")
}else{
  
  if(opt$nmethod=="log" & opt$imethod=="cca" & opt$rdsfile==""){
    sc.combined <- clustering(opt$dims,opt$resolution,"integrated","pca")
  }
  if(opt$nmethod=="log" & opt$imethod=="harmony" & opt$rdsfile==""){
    sc.combined <- clustering(opt$dims,opt$resolution,"RNA","harmony")
  }
  if(opt$nmethod=="sct" & opt$imethod=="harmony" & opt$rdsfile==""){
    sc.combined <- clustering(opt$dims,opt$resolution,"SCT","harmony")
  }
  if(opt$nmethod=="log" & opt$imethod=="false" & opt$rdsfile==""){
    sc.combined <- clustering(opt$dims,opt$resolution,"RNA","pca")
  }  
  if(opt$nmethod=="sct" & opt$imethod=="false" & opt$rdsfile==""){
    sc.combined <- clustering(opt$dims,opt$resolution,"SCT","pca")
  }    
   
}

#作图
plot_cluster <- function(){
  if(opt$skiptsne=="false"){
    ## tsne
    p1 <- DimPlot(sc.combined, reduction = "tsne", group.by = "replicate", pt.size = 0.1, raster = F )
    p2 <- DimPlot(sc.combined, reduction = "tsne", group.by = "group", pt.size = 0.1, raster = F )
    p3 <- DimPlot(sc.combined, reduction = "tsne", label = TRUE, pt.size = 0.1, raster = F )
    # p1 <- AugmentPlot(plot = p1)
    # p2 <- AugmentPlot(plot = p2)

    pdf(file=paste("tsne_cluster","pdf",sep="."),width = 18,height = 5)
    # p1 <- DimPlot(sc.combined, reduction = "tsne", pt.size = 0.1, label = TRUE) + ggtitle(label = "FIt-SNE")
    # p1 <- AugmentPlot(plot = p1)
    print((p1+p2+p3) + theme_paper)
    dev.off()

    png(file = paste("tsne_cluster",".png",sep=""),width = 450,height = 125,units = "mm",res = 300,pointsize = 10)
    par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
    # p1 <- DimPlot(sc.combined, reduction = "tsne", pt.size = 0.1, label = TRUE) + ggtitle(label = "FIt-SNE")
    # p1 <- AugmentPlot(plot = p1)
    print((p1+p2+p3) + theme_paper)
    dev.off()
  }

  ## umap
  p1 <- DimPlot(sc.combined, reduction = "umap", group.by = "replicate", pt.size = 0.1, raster = F )
  p2 <- DimPlot(sc.combined, reduction = "umap", group.by = "group", pt.size = 0.1, raster = F )
  p3 <- DimPlot(sc.combined, reduction = "umap", label = TRUE, pt.size = 0.1, raster = F )
  # p1 <- AugmentPlot(plot = p1)
  # p2 <- AugmentPlot(plot = p2)

  pdf(file=paste("umap_cluster","pdf",sep="."),width = 18,height = 5)
  print((p1+p2+p3) + theme_paper)
  # p1 <- DimPlot(sc.combined, reduction = "umap", pt.size = 0.1, label = TRUE) + ggtitle(label = "UMAP")
  # p1 <- AugmentPlot(plot = p1)
  # (p1)
  dev.off()

  png(file = paste("umap_cluster",".png",sep=""),width = 450,height = 125,units = "mm",res = 300,pointsize = 10)
  par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
  print(p1+p2+p3)
  # p1 <- DimPlot(sc.combined, reduction = "umap", pt.size = 0.1, label = TRUE) + ggtitle(label = "UMAP")
  # p1 <- AugmentPlot(plot = p1)
  # (p1)
  dev.off()

  if(opt$skiptsne=="false"){
    ## tsne和umap一起展示
    p1 <- DimPlot(sc.combined, reduction = "tsne", pt.size = 0.1, label = TRUE, raster = F ) + ggtitle(label = "FIt-SNE")
    p2 <- DimPlot(sc.combined, reduction = "umap", pt.size = 0.1, label = TRUE, raster = F ) + ggtitle(label = "UMAP")
    # p1 <- AugmentPlot(plot = p1)
    # p2 <- AugmentPlot(plot = p2)

    pdf(file=paste("tsne_umap","pdf",sep="."),width = 12,height = 5)
    print((p1 + p2) + theme_paper)
    dev.off()

    png(file = paste("tsne_umap",".png",sep=""),width = 300,height = 125,units = "mm",res = 300,pointsize = 10)
    par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
    print((p1 + p2) + theme_paper)
    dev.off()
  }
  #############
  ## 按照group展示cluster
  p1 <- DimPlot(sc.combined, reduction = "umap", group.by = "replicate", split.by = "group", pt.size = 0.1, label = TRUE, raster = F )
  # p1 <- AugmentPlot(plot = p1)

  pdf(file=paste("umap_byreplicate_bygroup","pdf",sep="."),width = 4 * groupnum, height = 4)
  print(p1)
  dev.off()

  png(file = paste("umap_byreplicate_bygroup",".png",sep=""),width = 100 * groupnum, height = 100,units = "mm",res = 300,pointsize = 10)
  par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
  print(p1)
  dev.off()

  if(opt$skiptsne=="false"){
    p1 <- DimPlot(sc.combined, reduction = "tsne", group.by = "replicate", split.by = "group", pt.size = 0.1, label = TRUE, raster = F )
    # p1 <- AugmentPlot(plot = p1)

    pdf(file=paste("tsne_byreplicate_bygroup","pdf",sep="."),width = 4 * groupnum, height = 4)
    print(p1)
    dev.off()

    png(file = paste("tsne_byreplicate_bygroup",".png",sep=""),width = 100 * groupnum, height = 100,units = "mm",res = 300,pointsize = 10)
    par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
    print(p1)
    dev.off()
  }
  ##############
  p1 <- DimPlot(sc.combined, reduction = "umap", split.by = "replicate", pt.size = 0.1, label = TRUE, raster = F )
  # p1 <- AugmentPlot(plot = p1)

  pdf(file=paste("umap_bycluster_byreplicate","pdf",sep="."),width = 4 * samplenum, height = 4)
  print(p1)
  dev.off()

  png(file = paste("umap_bycluster_byreplicate",".png",sep=""),width = 100 * samplenum, height = 100,units = "mm",res = 300,pointsize = 10)
  par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
  print(p1)
  dev.off()

  if(opt$skiptsne=="false"){
    p1 <- DimPlot(sc.combined, reduction = "tsne", split.by = "replicate", pt.size = 0.1, label = TRUE, raster = F )
    # p1 <- AugmentPlot(plot = p1)

    pdf(file=paste("tsne_bycluster_byreplicate","pdf",sep="."),width = 4 * groupnum, height = 4)
    print(p1)
    dev.off()

    png(file = paste("tsne_bycluster_byreplicate",".png",sep=""),width = 100 * groupnum, height = 100,units = "mm",res = 300,pointsize = 10)
    par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
    print(p1)
    dev.off()
  }
  ###############
  p1 <- DimPlot(sc.combined, reduction = "umap", split.by = "group", pt.size = 0.1, label = TRUE, raster = F )
  # p1 <- AugmentPlot(plot = p1)

  pdf(file=paste("umap_bycluster_bygroup","pdf",sep="."),width = 4 * groupnum, height = 4)
  print(p1)
  dev.off()

  png(file = paste("umap_bycluster_bygroup",".png",sep=""),width = 100 * groupnum, height = 100,units = "mm",res = 300,pointsize = 10)
  par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
  print(p1)
  dev.off()

  if(opt$skiptsne=="false"){
    p1 <- DimPlot(sc.combined, reduction = "tsne", split.by = "group", pt.size = 0.1, label = TRUE, raster = F )
    # p1 <- AugmentPlot(plot = p1)

    pdf(file=paste("tsne_bycluster_bygroup","pdf",sep="."),width = 4 * groupnum, height = 4)
    print(p1)
    dev.off()

    png(file = paste("tsne_bycluster_bygroup",".png",sep=""),width = 100 * groupnum, height = 100,units = "mm",res = 300,pointsize = 10)
    par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
    print(p1)
    dev.off()
  }
  #Phase umap and tsne plot
  tryCatch(
    {
      p1 <- DimPlot(sc.combined, reduction = "tsne", group.by = "Phase", pt.size = 0.1, raster = F )
      p2 <- DimPlot(sc.combined, reduction = "umap", group.by = "Phase", pt.size = 0.1, raster = F )      
      ggsave("tsne_umap_byphase.png",p1|p2,width = 9,height = 4,dpi = 300)
      ggsave("tsne_umap_byphase.pdf",p1|p2,width = 9,height = 4)
    },
    # warning = function(w){message("warning")},
    error = function(e){message("error")}
  )
}

plot_cluster()

## 回到工作根目录
setwd(workdir)

#####################################################################################
## 细胞注释
#####################################################################################
sctype_dir = paste(workdir,"/","sctype", sep="")
if(!file.exists(sctype_dir)){
    dir.create(sctype_dir)
}
setwd(sctype_dir)

scRNAseqData <- sc.combined
DefaultAssay(scRNAseqData) <- "RNA"

if(opt$cellanno == "" & opt$skipanno == "false"){
  
  # load libraries and functions
  lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
  source("/public/pipeline/ablifepipe/latest/src/SingleCell/sctype/gene_sets_prepare.R"); 
  source("/public/pipeline/ablifepipe/latest/src/SingleCell/sctype/sctype_score_.R")

  # get cell-type-specific gene sets from our in-built database (DB)
  # gs_list = gene_sets_prepare("/users/chengc/work/pipeline/ablifepipe/src/SingleCell/sctype/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
  # gs_list = gene_sets_prepare("/users/chengc/work/pipeline/ablifepipe/src/SingleCell/sctype/ScTypeDB_full.xlsx", "Immune system") # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Embryo,Gastrointestinal tract,Muscle,Skin

  # DB file
  db_ = opt$sctypedb
  tissue = opt$tissue # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Embryo,Gastrointestinal tract,Muscle,Skin,Heart,Ovary,Testis,White adipose tissue,Teeth
  # tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Embryo,Gastrointestinal tract,Muscle,Skin,Heart,Ovary,Testis,White adipose tissue,Teeth

  # prepare gene sets
  gs_list = gene_sets_prepare(db_, tissue)

  # scRNAseqData <- FindVariableFeatures(scRNAseqData, selection.method = "vst", nfeatures = 2000)

  # # scale and run PCA
  # scRNAseqData <- ScaleData(scRNAseqData, features = rownames(scRNAseqData))

  # # scale and run PCA
  scRNAseqData <- ScaleData(scRNAseqData)

  es.max = sctype_score(scRNAseqData = scRNAseqData[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

  # es.max = sctype_score(scRNAseqData = scRNAseqData[[opt$defaultassay]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

  # View results, cell-type by cell matrix. See the complete example below
  # View(es.max)

  # saveRDS(es.max, file = "es.max.rds")

  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(scRNAseqData@meta.data$seurat_clusters), function(cl){
      es.max.cl = sort(rowSums(es.max[ ,rownames(scRNAseqData@meta.data[scRNAseqData@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNAseqData@meta.data$seurat_clusters==cl)), 10)
  }))

  cL_resutls <- cL_resutls %>%
              dplyr::arrange(cluster, -scores)

  cat('',file=paste('cL_resutls','.xls',sep=''))
  write.table(cL_resutls,file=paste('cL_resutls','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[,1:4])
  cat('',file=paste('sctype_scores','.xls',sep=''))
  write.table(sctype_scores,file=paste('sctype_scores','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)
  # We can also overlay the identified cell types on UMAP plot:

  scRNAseqData@meta.data$customclassif = ""
  sc.combined@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    scRNAseqData@meta.data$customclassif[scRNAseqData@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
    sc.combined@meta.data$customclassif[sc.combined@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
}else if(opt$cellanno == "" & opt$skipanno == "true"){
  sc.combined@meta.data$customclassif = "None"
  setwd(workdir)
  unlink(sctype_dir,recursive = T)
}else{
  sctype_scores=read.table(opt$cellanno,header = T,com = '',sep = "\t",quote = NULL,check.names = F)

  scRNAseqData@meta.data$customclassif = ""
  sc.combined@meta.data$customclassif = ""

  for(j in unique(sctype_scores$cluster)){
      cl_type = sctype_scores[sctype_scores$cluster==j,]
      print(j)
      scRNAseqData@meta.data$customclassif[scRNAseqData@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
      sc.combined@meta.data$customclassif[sc.combined@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
}

#作图
plot_cluster_sctype <- function(){
  ## tsne
  p1 <- DimPlot(scRNAseqData, reduction = "tsne", group.by = "group", pt.size = 0.1, raster = F )
  p2 <- DimPlot(scRNAseqData, reduction = "tsne", label = TRUE, pt.size = 0.1, raster = F )
  p3 <- DimPlot(scRNAseqData, reduction = "tsne", repel = TRUE, group.by = 'customclassif', label = TRUE, pt.size = 0.1, raster = F )
  # p1 <- AugmentPlot(plot = p1)
  # p2 <- AugmentPlot(plot = p2)

  pdf(file=paste("tsne_cluster_sctype","pdf",sep="."),width = 22,height = 5)
  # p1 <- DimPlot(scRNAseqData, reduction = "tsne", pt.size = 0.1, label = TRUE) + ggtitle(label = "FIt-SNE")
  # p1 <- AugmentPlot(plot = p1)
  print((p1+p2+p3) + theme_paper)
  dev.off()

  png(file = paste("tsne_cluster_sctype",".png",sep=""),width = 550,height = 125,units = "mm",res = 300,pointsize = 10)
  par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
  # p1 <- DimPlot(scRNAseqData, reduction = "tsne", pt.size = 0.1, label = TRUE) + ggtitle(label = "FIt-SNE")
  # p1 <- AugmentPlot(plot = p1)
  print((p1+p2+p3) + theme_paper)
  dev.off()

  ## umap
  p1 <- DimPlot(scRNAseqData, reduction = "umap", group.by = "group", pt.size = 0.1, raster = F )
  p2 <- DimPlot(scRNAseqData, reduction = "umap", label = TRUE, pt.size = 0.1, raster = F )
  p3 <- DimPlot(scRNAseqData, reduction = "umap", repel = TRUE, group.by = 'customclassif', label = TRUE, pt.size = 0.1, raster = F )
  # p1 <- AugmentPlot(plot = p1)
  # p2 <- AugmentPlot(plot = p2)

  pdf(file=paste("umap_cluster_sctype","pdf",sep="."),width = 22,height = 5)
  print((p1+p2+p3) + theme_paper)
  # p1 <- DimPlot(scRNAseqData, reduction = "umap", pt.size = 0.1, label = TRUE) + ggtitle(label = "UMAP")
  # p1 <- AugmentPlot(plot = p1)
  # (p1)
  dev.off()

  png(file = paste("umap_cluster_sctype",".png",sep=""),width = 550,height = 125,units = "mm",res = 300,pointsize = 10)
  par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
  print(p1+p2+p3)
  # p1 <- DimPlot(scRNAseqData, reduction = "umap", pt.size = 0.1, label = TRUE) + ggtitle(label = "UMAP")
  # p1 <- AugmentPlot(plot = p1)
  # (p1)
  dev.off()

  ## tsne和umap一起展示
  p1 <- DimPlot(scRNAseqData, reduction = "tsne", repel = TRUE, group.by = 'customclassif', pt.size = 0.1, label = TRUE, raster = F ) + ggtitle(label = "FIt-SNE")
  p2 <- DimPlot(scRNAseqData, reduction = "umap", repel = TRUE, group.by = 'customclassif', pt.size = 0.1, label = TRUE, raster = F ) + ggtitle(label = "UMAP")
  # p1 <- AugmentPlot(plot = p1)
  # p2 <- AugmentPlot(plot = p2)

  pdf(file=paste("tsne_umap_sctype","pdf",sep="."),width = 22,height = 5)
  print((p1 + p2) + theme_paper)
  dev.off()

  png(file = paste("tsne_umap_sctype",".png",sep=""),width = 550,height = 125,units = "mm",res = 300,pointsize = 10)
  par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
  print((p1 + p2) + theme_paper)
  dev.off()

}

if(opt$skipanno == "false"){
  plot_cluster_sctype()
}

rm(scRNAseqData)
# scRNAseqData<-ls()

#####################################################################################
## 统计和输出
#####################################################################################
## 统计
setwd(sc_cluster_dir)

clu = cbind(Idents(sc.combined),sc.combined[["replicate"]],sc.combined[["group"]],sc.combined[["customclassif"]])

cat(c('CellIndex\tClusterID\tsample\tgroup\tcellclass\n'),file=paste('Clusters','.xls',sep=''))
write.table(clu,file=paste('Clusters','.xls',sep=''),
            append=T,quote=F,sep='\t',row.names = TRUE,col.names = FALSE)
system("perl /public/pipeline/ablifepipe/latest/src/SingleCell/tools/sc_stat_cluster.pl")

## 回到工作根目录
setwd(workdir)

sc.combined$cellcluster <- paste("C", Idents(sc.combined), sep = "")
sc.combined$cellgroup <- paste("C", Idents(sc.combined), sc.combined$group, sep = "_")

## 保存
saveRDS(sc.combined, file = "cells_analysis.rds")

#
sc <- sc.combined
# 将 sc@meta.data 中所有因子型的列转换为字符型
sc@meta.data <- sc@meta.data %>% mutate(across(where(is.factor), as.character))
sc <- DietSeurat(sc,dimreducs=c("pca","tsne","umap"))
DefaultAssay(sc) <- "RNA"
SaveH5Seurat(sc, filename = "sc.h5Seurat",overwrite = TRUE)
Convert("sc.h5Seurat", dest = "h5ad",overwrite = TRUE)
system("rm -rf sc.h5Seurat")
rm(sc)

if(opt$nmethod == "sct"){
  #准备cell2location的文件
  scRNA <- sc.combined
  DefaultAssay(scRNA) <- "RNA"
  scRNA <- DietSeurat(scRNA,assays = "RNA")
  # scRNA[['SCT']] <- NULL
  scRNA@assays$RNA@data <- scRNA@assays$RNA@counts
  SaveH5Seurat(scRNA, filename = "sc_cell2location.h5Seurat",overwrite = TRUE)
  Convert("sc_cell2location.h5Seurat", dest = "h5ad",overwrite = TRUE)
  system("rm -rf sc_cell2location.h5Seurat")
  rm(scRNA)
}


#####################################################################################
## marker基因鉴定
#####################################################################################

#marker基因鉴定就用RNA槽的数据
DefaultAssay(sc.combined) <- "RNA"

if(opt$mmethod == "seurat"){

  ####1. FindAllMarkers 找到各组中所有的markers

  # https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_marker_identification.html

  # find markers for every cluster compared to all remaining cells, report only the positive ones

  combined_markers_dir = paste(workdir,"/","combined_markers", sep="")
  if(!file.exists(combined_markers_dir)){
    dir.create(combined_markers_dir)
  }
  setwd(combined_markers_dir)

  DefaultAssay(sc.combined) <- "RNA"
  Idents(sc.combined) <- "seurat_clusters"
  combined_markers <- FindAllMarkers(object = sc.combined,
                            test.use = "wilcox",
                            slot = "data",
                            only.pos = TRUE,
                            min.pct = 0,
                            logfc.threshold = 0)

  # Order the rows by p-adjusted values
  combined_markers <- combined_markers[ , c(6, 7, 2:4, 1, 5)]
  combined_markers <- combined_markers %>%
          dplyr::arrange(cluster, -avg_log2FC, p_val)

  cat(c(''),file=paste('combined_markers_export','.xls',sep=''))
  write.table(combined_markers,file=paste('combined_markers_export','.xls',sep=''),
              append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

  # saveRDS(combined_markers, file = "../combined_markers.rds")

  # 提取每个类群排名靠前的3个或10个保守标记物
  top3 <- combined_markers %>% 
    group_by(cluster) %>% 
    top_n(n = 3, 
          wt = avg_log2FC)

  cat(c(''),file=paste('top3_combined_markers_export','.xls',sep=''))
  write.table(top3,file=paste('top3_combined_markers_export','.xls',sep=''),
              append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

  top10 <- combined_markers %>% 
    group_by(cluster) %>% 
    top_n(n = 10, 
          wt = avg_log2FC)
  cat(c(''),file=paste('top10_combined_markers_export','.xls',sep=''))
  write.table(top10,file=paste('top10_combined_markers_export','.xls',sep=''),
              append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

  top50 <- combined_markers %>% 
    group_by(cluster) %>% 
    top_n(n = 50, 
          wt = avg_log2FC)
  cat(c(''),file=paste('top50_combined_markers_export','.xls',sep=''))
  write.table(top50,file=paste('top50_combined_markers_export','.xls',sep=''),
              append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

  top100 <- combined_markers %>% 
    group_by(cluster) %>% 
    top_n(n = 100, 
          wt = avg_log2FC)
  cat(c(''),file=paste('top100_combined_markers_export','.xls',sep=''))
  write.table(top100,file=paste('top100_combined_markers_export','.xls',sep=''),
              append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

  deg <- subset(combined_markers,avg_log2FC>=0.5)
  deg <- subset(deg,pct.1>=pct.2)
  deg <- subset(deg,pct.1>=0.2)
  deg <- subset(deg,p_val_adj<=0.05)

  cat(c(''),file=paste('deg_combined_markers_export','.xls',sep=''))
  write.table(deg,file=paste('deg_combined_markers_export','.xls',sep=''),
              append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

  ## plot top3 markers heatmap

  # diffgenes=subset(top3,avg_log2FC>=1|avg_log2FC<=-1)
  # diffgenes=subset(diffgenes,pct.1-pct.2>=0.3|pct.1-pct.2<=-0.3)
  # diffgenes=subset(diffgenes,p_val_adj<=0.05)
  diffgenes  = top3
  showgene = unique(diffgenes$gene)
  showgenenum = length(showgene)
  prefix = paste("top3_combined_markers",sep="")

  nowdir = getwd()
  # show_genes(showgene,showgenenum,groupnum,clusternum,prefix,FALSE)
  # if(opt$nmethod == "log" && opt$imethod == "cca"){
  #   system("p3 /users/chengc/work/pipeline/ablifegraphic/honeycomb/apps/single-cell/scanpy_plot.py -s ../sc.h5ad -m ./top3_combined_markers_export.xls -O top3_scanpy")
  # }else{
  #   system("p3 /users/chengc/work/pipeline/ablifegraphic/honeycomb/apps/single-cell/scanpy_plot2.py -s ../sc.h5ad -m ./top3_combined_markers_export.xls -O top3_scanpy")
  # }

  #作图
  Idents(sc.combined) <- "seurat_clusters"
  tryCatch(
    {
      #按照top3作图
      #heatmap
      sc.combined_sub <- ScaleData(sc.combined, features = top3$gene)
      p1 <- DoHeatmap(sc.combined_sub, features = top3$gene, size = 5) + 
        theme(axis.text.y = element_text(size = 10,color = "black")) + NoLegend()
      ggsave("top3_makrer_Heatmap.png", p1, width = 15, height = 15)
      ggsave("top3_makrer_Heatmap.pdf", p1, width = 15, height = 15)
      
      #VlnPlot
      p1 <- VlnPlot(
        object = sc.combined,
        features = top3$gene,
        pt.size = 0.1,
        stack = T,
        flip = T
      ) + NoLegend()
      ggsave("top3_makrer_VlnPlot.png",p1,width = 15,height = 15,dpi = 300)
      ggsave("top3_makrer_VlnPlot.pdf",p1,width = 15,height = 15)

      #DotPlot
      p1 <- DotPlot(
        object = sc.combined,
        features = unique(top3$gene),
      ) + RotatedAxis()
      ggsave("top3_makrer_DotPlot.png",p1,width = 15,height = 12,dpi = 300)
      ggsave("top3_makrer_DotPlot.pdf",p1,width = 15,height = 12)
            
    },
    error = function(e){message("markers不存在")}
  )  

  setwd(nowdir)
  setwd(workdir)

  ##combined_markers_bycelltype
  if(opt$skipanno == "false"){
    combined_markers_bycelltype_dir = paste(workdir,"/","combined_markers_bycelltype", sep="")
    if(!file.exists(combined_markers_bycelltype_dir)){
      dir.create(combined_markers_bycelltype_dir)
    }
    setwd(combined_markers_bycelltype_dir)

    DefaultAssay(sc.combined) <- "RNA"
    Idents(sc.combined) <- "customclassif"
    combined_markers <- FindAllMarkers(object = sc.combined,
                              test.use = "wilcox",
                              slot = "data",
                              only.pos = TRUE,
                              min.pct = 0,
                              logfc.threshold = 0)
    # Order the rows by p-adjusted values
    combined_markers <- combined_markers[ , c(6, 7, 2:4, 1, 5)]
    combined_markers <- combined_markers %>%
            dplyr::arrange(cluster, -avg_log2FC, p_val)

    cat(c(''),file=paste('combined_markers_export','.xls',sep=''))
    write.table(combined_markers,file=paste('combined_markers_export','.xls',sep=''),
                append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

    # saveRDS(combined_markers, file = "../combined_markers.rds")

    # 提取每个类群排名靠前的3个或10个保守标记物
    top3 <- combined_markers %>% 
      group_by(cluster) %>% 
      top_n(n = 3, 
            wt = avg_log2FC)

    cat(c(''),file=paste('top3_combined_markers_export','.xls',sep=''))
    write.table(top3,file=paste('top3_combined_markers_export','.xls',sep=''),
                append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

    top10 <- combined_markers %>% 
      group_by(cluster) %>% 
      top_n(n = 10, 
            wt = avg_log2FC)
    cat(c(''),file=paste('top10_combined_markers_export','.xls',sep=''))
    write.table(top10,file=paste('top10_combined_markers_export','.xls',sep=''),
                append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

    top50 <- combined_markers %>% 
      group_by(cluster) %>% 
      top_n(n = 50, 
            wt = avg_log2FC)
    cat(c(''),file=paste('top50_combined_markers_export','.xls',sep=''))
    write.table(top50,file=paste('top50_combined_markers_export','.xls',sep=''),
                append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

    top100 <- combined_markers %>% 
      group_by(cluster) %>% 
      top_n(n = 100, 
            wt = avg_log2FC)
    cat(c(''),file=paste('top100_combined_markers_export','.xls',sep=''))
    write.table(top100,file=paste('top100_combined_markers_export','.xls',sep=''),
                append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

    deg <- subset(combined_markers,avg_log2FC>=0.5)
    deg <- subset(deg,pct.1>=pct.2)
    deg <- subset(deg,pct.1>=0.2)
    deg <- subset(deg,p_val_adj<=0.05)

    cat(c(''),file=paste('deg_combined_markers_export','.xls',sep=''))
    write.table(deg,file=paste('deg_combined_markers_export','.xls',sep=''),
                append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

    ## plot top3 markers heatmap

    # diffgenes=subset(top3,avg_log2FC>=1|avg_log2FC<=-1)
    # diffgenes=subset(diffgenes,pct.1-pct.2>=0.3|pct.1-pct.2<=-0.3)
    # diffgenes=subset(diffgenes,p_val_adj<=0.05)
    diffgenes  = top3
    showgene = unique(diffgenes$gene)
    showgenenum = length(showgene)
    prefix = paste("top3_combined_markers",sep="")

    nowdir = getwd()
    # show_genes(showgene,showgenenum,groupnum,clusternum,prefix,FALSE)
    # if(opt$nmethod == "log" && opt$imethod == "cca"){
    #   system("p3 /users/chengc/work/pipeline/ablifegraphic/honeycomb/apps/single-cell/scanpy_plot.py -s ../sc.h5ad -m ./top3_combined_markers_export.xls -O top3_scanpy")
    # }else{
    #   system("p3 /users/chengc/work/pipeline/ablifegraphic/honeycomb/apps/single-cell/scanpy_plot2.py -s ../sc.h5ad -m ./top3_combined_markers_export.xls -O top3_scanpy")
    # }

    sc.combined$customclassif <- factor(sc.combined$customclassif,levels=mixedsort(unique(sc.combined$customclassif))) 
    Idents(sc.combined) <- "customclassif"   
    tryCatch(
      {
        #按照top3作图
        #heatmap
        sc.combined_sub <- ScaleData(sc.combined, features = top3$gene)
        p1 <- DoHeatmap(sc.combined_sub, features = top3$gene, size = 5) + 
          theme(axis.text.y = element_text(size = 10,color = "black")) + NoLegend()
        ggsave("top3_makrer_Heatmap.png", p1, width = 15, height = 15)
        ggsave("top3_makrer_Heatmap.pdf", p1, width = 15, height = 15)
        
        #VlnPlot
        p1 <- VlnPlot(
          object = sc.combined,
          features = top3$gene,
          pt.size = 0.1,
          stack = T,
          flip = T
        ) + NoLegend()
        ggsave("top3_makrer_VlnPlot.png",p1,width = 15,height = 15,dpi = 300)
        ggsave("top3_makrer_VlnPlot.pdf",p1,width = 15,height = 15)

        #DotPlot
        p1 <- DotPlot(
          object = sc.combined,
          features = unique(top3$gene),
        ) + RotatedAxis()
        ggsave("top3_makrer_DotPlot.png",p1,width = 15,height = 12,dpi = 300)
        ggsave("top3_makrer_DotPlot.pdf",p1,width = 15,height = 12)
              
      },
      error = function(e){message("markers不存在")}
    )  

    setwd(nowdir)
    setwd(workdir)
  }

}else if(opt$mmethod == "presto"){

  #一个简单的方法是执行Wilcoxon秩和测试，presto包实现了一个可以在Seurat对象上运行的极快的Wilcoxon测试。
  #使用FindAllMarkers非常慢，即便是开了多线程也非常慢
  #使用presto包做极快的Wilcoxon测试

  ########################by seurat_clusters###############################
  combined_markers_dir = paste(workdir,"/","combined_markers", sep="")
  if(!file.exists(combined_markers_dir)){
    dir.create(combined_markers_dir)
  }
  setwd(combined_markers_dir)
  #运行
  all_markers <- wilcoxauc(
    sc.combined,
    group_by = "seurat_clusters",
    assay = "data",
    groups_use = NULL,
    seurat_assay = "RNA",
  ) %>% 
    dplyr::mutate(avg_log2FC = log(exp(logFC),base = 2), pct.1 = round(pct_in/100,3), pct.2 = round(pct_out/100,3)) %>%
    dplyr::select(group,feature,avg_log2FC,pct.1,pct.2,pval,padj) 
    
  colnames(all_markers) <- c("cluster","gene","avg_log2FC","pct.1","pct.2","p_val","p_val_adj") 

  all_markers <- dplyr::filter(all_markers, avg_log2FC > 0 & p_val_adj < 0.05) %>% dplyr::arrange(cluster,-avg_log2FC,-pct.1,p_val)
  all_markers <- all_markers[mixedorder(all_markers$cluster),] #按照自然排序
  write.table(all_markers,file = "combined_markers_export.xls",sep = "\t",quote = F,row.names = F,col.names = T)   

  #差异marker
  filter_markers <- dplyr::filter(all_markers,avg_log2FC > 0.5 & p_val_adj < 0.05 & pct.1 >= 0.2 & pct.1 >= pct.2)
  write.table(filter_markers,file = "deg_combined_markers_export.xls",sep = "\t",quote = F,row.names = F,col.names = T)

  #提取每个类群排名靠前的top marker
  top3_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
  write.table(top3_markers,file = "top3_combined_markers_export.xls",sep = "\t",quote = F,row.names = F,col.names = T)

  top10_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  write.table(top10_markers,file = "top10_combined_markers_export.xls",sep = "\t",quote = F,row.names = F,col.names = T)

  top50_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
  write.table(top50_markers,file = "top50_combined_markers_export.xls",sep = "\t",quote = F,row.names = F,col.names = T)  

  top100_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
  write.table(top100_markers,file = "top100_combined_markers_export.xls",sep = "\t",quote = F,row.names = F,col.names = T)

  #scanpy作图 很慢
  # if(opt$nmethod == "log" && opt$imethod == "cca"){
  #   system("p3 /users/chengc/work/pipeline/ablifegraphic/honeycomb/apps/single-cell/scanpy_plot.py -s ../sc.h5ad -m ./top3_combined_markers_export.xls -O top3_scanpy")
  # }else{
  #   system("p3 /users/chengc/work/pipeline/ablifegraphic/honeycomb/apps/single-cell/scanpy_plot2.py -s ../sc.h5ad -m ./top3_combined_markers_export.xls -O top3_scanpy")
  # }

  #作图
  Idents(sc.combined) <- "seurat_clusters"
  tryCatch(
    {
      #按照top3作图
      #heatmap
      sc.combined_sub <- ScaleData(sc.combined, features = top3_markers$gene)
      p1 <- DoHeatmap(sc.combined_sub, features = top3_markers$gene, size = 5) + 
        theme(axis.text.y = element_text(size = 10,color = "black")) + NoLegend()
      ggsave("top3_makrer_Heatmap.png", p1, width = 15, height = 15)
      ggsave("top3_makrer_Heatmap.pdf", p1, width = 15, height = 15)
      
      #VlnPlot
      p1 <- VlnPlot(
        object = sc.combined,
        features = top3_markers$gene,
        pt.size = 0.1,
        stack = T,
        flip = T
      ) + NoLegend()
      ggsave("top3_makrer_VlnPlot.png",p1,width = 15,height = 15,dpi = 300)
      ggsave("top3_makrer_VlnPlot.pdf",p1,width = 15,height = 15)

      #DotPlot
      p1 <- DotPlot(
        object = sc.combined,
        features = unique(top3_markers$gene)
      ) + RotatedAxis()
      ggsave("top3_makrer_DotPlot.png",p1,width = 15,height = 12,dpi = 300)
      ggsave("top3_makrer_DotPlot.pdf",p1,width = 15,height = 12)
            
    },
    error = function(e){message("markers不存在")}
  )

  setwd(workdir)

  ########################by customclassif###############################
  if(opt$skipanno == "false"){
    combined_markers_dir = paste(workdir,"/","combined_markers_bycelltype", sep="")
    if(!file.exists(combined_markers_dir)){
      dir.create(combined_markers_dir)
    }
    setwd(combined_markers_dir)
    #运行
    all_markers <- wilcoxauc(
      sc.combined,
      group_by = "customclassif",
      assay = "data",
      groups_use = NULL,
      seurat_assay = "RNA",
    ) %>% 
      dplyr::mutate(avg_log2FC = log(exp(logFC),base = 2), pct.1 = round(pct_in/100,3), pct.2 = round(pct_out/100,3)) %>%
      dplyr::select(group,feature,avg_log2FC,pct.1,pct.2,pval,padj) 
      
    colnames(all_markers) <- c("cluster","gene","avg_log2FC","pct.1","pct.2","p_val","p_val_adj") 

    all_markers <- dplyr::filter(all_markers, avg_log2FC > 0 & p_val_adj < 0.05) %>% dplyr::arrange(cluster,-avg_log2FC,-pct.1,p_val)
    all_markers <- all_markers[mixedorder(all_markers$cluster),] #按照自然排序
    write.table(all_markers,file = "combined_markers_export.xls",sep = "\t",quote = F,row.names = F,col.names = T)   

    #差异marker
    filter_markers <- dplyr::filter(all_markers,avg_log2FC > 0.5 & p_val_adj < 0.05 & pct.1 >= 0.2 & pct.1 >= pct.2)
    write.table(filter_markers,file = "deg_combined_markers_export.xls",sep = "\t",quote = F,row.names = F,col.names = T)

    #提取每个类群排名靠前的top marker
    top3_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
    write.table(top3_markers,file = "top3_combined_markers_export.xls",sep = "\t",quote = F,row.names = F,col.names = T)

    top10_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    write.table(top10_markers,file = "top10_combined_markers_export.xls",sep = "\t",quote = F,row.names = F,col.names = T)

    top50_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
    write.table(top50_markers,file = "top50_combined_markers_export.xls",sep = "\t",quote = F,row.names = F,col.names = T)  

    top100_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
    write.table(top100_markers,file = "top100_combined_markers_export.xls",sep = "\t",quote = F,row.names = F,col.names = T)

    # #scanpy作图
    # if(opt$nmethod == "log" && opt$imethod == "cca"){
    #   system("p3 /users/chengc/work/pipeline/ablifegraphic/honeycomb/apps/single-cell/scanpy_plot.py -s ../sc.h5ad -m ./top3_combined_markers_export.xls -O top3_scanpy")
    # }else{
    #   system("p3 /users/chengc/work/pipeline/ablifegraphic/honeycomb/apps/single-cell/scanpy_plot2.py -s ../sc.h5ad -m ./top3_combined_markers_export.xls -O top3_scanpy")
    # }
 
    #作图
    sc.combined$customclassif <- factor(sc.combined$customclassif,levels=mixedsort(unique(sc.combined$customclassif)))    
    Idents(sc.combined) <- "customclassif"
    tryCatch(
      {
        #按照top3作图
        #heatmap
        sc.combined_sub <- ScaleData(sc.combined, features = top3_markers$gene)
        p1 <- DoHeatmap(sc.combined_sub, features = top3_markers$gene, size = 5) + 
          theme(axis.text.y = element_text(size = 10,color = "black")) + NoLegend()
        ggsave("top3_makrer_Heatmap.png", p1, width = 15, height = 15)
        ggsave("top3_makrer_Heatmap.pdf", p1, width = 15, height = 15)
        
        #VlnPlot
        p1 <- VlnPlot(
          object = sc.combined,
          features = top3_markers$gene,
          pt.size = 0.1,
          stack = T,
          flip = T
        ) + NoLegend()
        ggsave("top3_makrer_VlnPlot.png",p1,width = 15,height = 15,dpi = 300)
        ggsave("top3_makrer_VlnPlot.pdf",p1,width = 15,height = 15)

        #DotPlot
        p1 <- DotPlot(
          object = sc.combined,
          features = unique(top3_markers$gene),
        ) + RotatedAxis()
        ggsave("top3_makrer_DotPlot.png",p1,width = 15,height = 12,dpi = 300)
        ggsave("top3_makrer_DotPlot.pdf",p1,width = 15,height = 12)
              
      },
      error = function(e){message("markers不存在")}
    ) 

    setwd(workdir)    

  }

}
