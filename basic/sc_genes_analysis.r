#!/usr/bin/env Rscript

###Set VERSION
VERSION = "3.0.0"

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1 && (args[1] == "-v" | args[1] == "--version")){
  message("Version: \n\t", VERSION, "\n")
  quit(status = 0)
}

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
# https://satijalab.org/seurat/pbmc3k_tutorial.html
# https://github.com/CostaLab/scrna_seurat_pipeline
#####################################################################################

#####################################################################################
#####参数获取
#####################################################################################
# windows系统检查当前程序路径
# script.dir <- dirname(sys.frame(1)$ofile)

library(tidyverse)
get_this_file <- function(){
    commandArgs() %>% 
       tibble::enframe(name=NULL) %>%
       tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
       dplyr::filter(key == "--file") %>%
       dplyr::pull(value)
}
script.dir <- get_this_file()

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(optparse))      ## Options
suppressPackageStartupMessages(library(futile.logger)) ## logger
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(doParallel))
library(SingleCellExperiment)
library(scater)
library(ggrepel)
library(grDevices)
library(gplots)
library(Seurat)
library(patchwork)
library(cowplot)
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
        c("-c","--comparelist"),
        action = "store",
        type = "character",
        default = "",
        help = "group1:group2,group2:group3,group3:group4"
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
        c("-i", "--infile"),
        action = "store",
        type = "character",
        help = "The cell analysis rds file"
    ),
    # make_option(
    #     c("-d", "--dims"),
    #     action = "store",
    #     type = "integer",
    #     default = 30,
    #     help = "PCA dims number"
    # ),
    make_option(
        c("--defaultassay"),
        action = "store",
        type = "character",
        default = "RNA",
        help = "default assay"
    ),
    make_option(
        c("-a", "--annofile"),
        action = "store",
        type = "character",
        help = "The Annotation file"
    ),
    make_option(
        c("--skipanno"),
        action = "store",
        type = "character",
        default = "false",
        help = "skip cell anno, true or false"
    ),
    make_option(
        c("--mmethod"),
        action = "store",
        type = "character",
        default = "seurat",
        help = "find marker methods, seurat or presto, default is seurat"
    ),
    make_option(
        c("--logfc"),
        action = "store",
        type = "double",
        default = 0.5,
        help = "deg logfc threshold"
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

# opt$dims <- 20
# opt$logfc <- 0.5
# opt$infile <- "/users/wuyf/public/Singlecell/pipeline/harmony/2_cells_analysis/cells_analysis.rds"
# opt$samplelist <- "ITP_1,HC_1"
# opt$listname <- "ITP,HC"
# opt$comparelist <- "ITP:HC"

listnames <- strsplit(opt$listname,",")[[1]]
groupnum = length(unique(listnames))
#####################################################################################
## preprocessing数据读取，数据对象建立
# We start by reading in the data. The Read10X function reads in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).
#####################################################################################

# Load the PBMC dataset

sc.combined <- readRDS(file = opt$infile)

workdir <- getwd()

#
if(opt$mmethod == "seurat"){

  #####################################################################################
  ## Identify differential expressed genes across conditions 差异基因分析
  # Now that we've aligned the stimulated and control cells, we can start to do comparative analyses and look at the differences induced by stimulation. One way to look broadly at these changes is to plot the average expression of both the stimulated and control cells and look for genes that are visual outliers on a scatter plot. Here, we take the average expression of both the stimulated and control naive T cells and CD14 monocyte populations and generate the scatter plots, highlighting genes that exhibit dramatic responses to interferon stimulation.
  #####################################################################################

  ## 鉴定差异基因
  clusternum <- length(unique(factor(Idents(sc.combined))))
  DefaultAssay(sc.combined) <- "RNA"
  sc.idents.bak <- Idents(sc.combined)
  # Idents(sc.combined) <- sc.idents.bak
  maxIdents<- unique(factor(Idents(sc.combined)))
  sc.combined$celltype.group <- paste(Idents(sc.combined), sc.combined$group, sep = "_")
  sc.combined$celltype <- Idents(sc.combined)
  Idents(sc.combined) <- "celltype.group"

  compares <- strsplit(opt$comparelist,",")[[1]]

  workdir <- getwd()

  deg_dir3 = paste(workdir,"/","deg_bycluster", sep="")
  if(!file.exists(deg_dir3)){
    dir.create(deg_dir3)
  }
  setwd(deg_dir3)

  source("/public/pipeline/ablifepipe/latest/src/SingleCell/allsample/sc_genes_visualizations_function.r")

  get_diff <- function(cid,test_group,ctrl_group){
    tryCatch(
        {
          FindMarkers(sc.combined,
                        ident.1 = paste(cid, test_group, sep = "_"),
                        ident.2 = paste(cid, ctrl_group, sep = "_"),
                        min.cells.group = 3,
                        logfc.threshold = 0,
                        min.pct = 0,
                        only.pos = FALSE) %>%
          rownames_to_column(var = "gene") %>%
          cbind(cluster = cid, .)
        },
        error=function(cond) {
          ret_code <<- -1
          data.frame()
        },
        finally={
        })
    
    }

  plot_scatter_by_group <- function(cid,test_group,ctrl_group){
      # cluster <- "CD4 Naive T"
      Idents(sc.combined) <- sc.idents.bak
      library(ggplot2)
      library(cowplot)
      theme_set(theme_cowplot())

      sub.cells <- subset(sc.combined, idents = cid)
      # seurat.subset <- subset(x = sub.cells, subset = (group == test_group & group == ctrl_group) & (another_condition == "Ambigious"))
      sub.cells <- subset(sub.cells,group == c(test_group,ctrl_group))
      # sub.cells <- sub.cells@meta.data
      # table(sub.cells$group)

      degdir = getwd()

      top5 = rbind(top5_up,top5_down)

      expscatter_dir = paste("","expscatter", sep="")
      if(!file.exists(expscatter_dir)){
        dir.create(expscatter_dir)
      }

      setwd(expscatter_dir)

      tryCatch(
        {
          Idents(sub.cells) <- "group"
          avg.sub.cells <- log1p(AverageExpression(sub.cells, verbose = FALSE)$RNA)
          avg.sub.cells.data=data.frame(avg.sub.cells)
          avg.sub.cells.data$gene <- rownames(avg.sub.cells.data)
          diffgenes  = subset(top5, cluster == cid)
          genes.to.label = unique(diffgenes$gene)
          print(genes.to.label)
          p1 <- ggplot(avg.sub.cells.data,aes_string(x=test_group,y=ctrl_group)) + 
              geom_point(
          aes(
            fill = ifelse(gene %in% genes.to.label, "blue", "black")
            ),
            shape = 21,
            alpha = 0.6,
          ) +

          scale_fill_manual(values = c("grey", "red"))+

          ggtitle(paste("cluster ",cid,sep=""))
          p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)+NoLegend()
          ggsave(file=paste("cluster",cid,"_top5_expscatter",".pdf",sep=""),width = 100,height = 100,units = "mm")
          ggsave(file=paste("cluster",cid,"_top5_expscatter",".png",sep=""),width = 100,height = 100,units = "mm", dpi=300)
        },
        error=function(cond) {
          ret_code <<- -1
        },
        finally={
          print(paste("cluster",cid,sep=""))
        })
      
      Idents(sc.combined) <- sc.idents.bak
      setwd(degdir)

      # showgene_dir = paste(degdir,"/top5_geneplot", sep="")
      # if(!file.exists(showgene_dir)){
      #   dir.create(showgene_dir)
      # }
      # setwd(showgene_dir)

      # # top5 单个基因的cluster方式展示
      # # top5 = rbind(top5,top5_down)
      # diffgenes2  = subset(top5_up, cluster == cid)
      # showgene = unique(diffgenes2$gene)
      # showgenenum = length(showgene)
      # prefix = paste("cluster",cid,"_up",sep="")
      # nowdir = getwd()
      # # if(showgenenum>0){
      # #   tryCatch(
      # #   {
      # #     show_genes(showgene,showgenenum,groupnum,clusternum,prefix,cid)
      # #   },
      # #   error=function(cond) {
      # #     ret_code <<- -1
      # #   },
      # #   finally={
      # #     print(paste("warning: cluster",cid,"_up",sep=""))
      # #   })
        
      # # }
      # # system("py3 /users/chengc/work/pipeline/ablifegraphic/honeycomb/apps/single-cell/scanpy_plot.py -s ../../sc.h5ad -m ../top3_combined_markers_export.xls")
      # setwd(nowdir)

      # diffgenes2  = subset(top5_down, cluster == cid)
      # showgene = unique(diffgenes2$gene)
      # showgenenum = length(showgene)
      # prefix = paste("cluster",cid,"_down",sep="")
      # nowdir = getwd()

      # # if(showgenenum>0){
      # #   tryCatch(
      # #   {
      # #     show_genes(showgene,showgenenum,groupnum,clusternum,prefix,cid)
      # #   },
      # #   error=function(cond) {
      # #     ret_code <<- -1
      # #   },
      # #   finally={
      # #     print(paste("warning: cluster",cid,"_down",sep=""))
      # #   })
      # # }
      # setwd(nowdir)
      
      # setwd(degdir)

  }

  compare_group <- function(gid){
      
      groups <- strsplit(compares[gid],":")[[1]]
      test_group = groups[1]
      ctrl_group = groups[2]

      deg_dir = paste(deg_dir3,"/",test_group,"_vs_", ctrl_group,"_deg", sep="")

      if(!file.exists(deg_dir)){
        dir.create(deg_dir)
      }

      setwd(deg_dir)

      # https://blog.csdn.net/qq_18055167/article/details/104437236
      # diff_genes <- map_dfr(c(maxIdents,test_group,ctrl_group), get_diff)
      Idents(sc.combined) <- "celltype.group"
      diff_genes <- pmap_dfr(list(maxIdents,test_group,ctrl_group), get_diff)
      # diff_genes <- future_mapply(get_diff,maxIdents,test_group,ctrl_group,future.seed = TRUE)

      diff_genes <- diff_genes[ , c(1:2,4:6, 3, 7)]
      diff_genes <- diff_genes %>%
              dplyr::arrange(cluster, p_val)

      cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_all_genes_export','.xls',sep=''))
      write.table(diff_genes,file=paste(test_group,'_vs_',ctrl_group,'_all_genes_export','.xls',sep=''),append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

      # saveRDS(diff_genes, file = "./diff_genes.rds")
      # saveRDS(conserved_markers, file = "./conserved_markers.rds")
      # saveRDS(combined_markers, file = "./combined_markers.rds")
      # 提取每个类群排名靠前的3个或10个保守标记物
      #上调

      
      diff_genes <- diff_genes %>%
              dplyr::arrange(cluster, -avg_log2FC, p_val)
      deg <- subset(diff_genes,avg_log2FC>=opt$logfc)
      deg <- subset(deg,pct.1>=0.2|pct.2>=0.2)
      deg <- subset(deg,p_val_adj<=0.05)

      cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_filter_up_deg','.xls',sep=''))
      write.table(deg,file=paste(test_group,'_vs_',ctrl_group,'_filter_up_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

      # deg_up <- subset(diff_genes,avg_log2FC>=0)
      deg_up <- deg
      top3_up <<- deg_up %>% 
      group_by(cluster) %>% 
      top_n(n = 3, 
              wt = avg_log2FC)
      cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top3_up_deg','.xls',sep=''))
      write.table(top3_up,file=paste(test_group,'_vs_',ctrl_group,'_top3_up_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)
      
      top5_up <<- deg_up %>% 
      group_by(cluster) %>% 
      top_n(n = 5, 
              wt = avg_log2FC)
      cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top5_up_deg','.xls',sep=''))
      write.table(top5_up,file=paste(test_group,'_vs_',ctrl_group,'_top5_up_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

      top10_up <<- deg_up %>% 
      group_by(cluster) %>% 
      top_n(n = 10, 
              wt = avg_log2FC)
      cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top10_up_deg','.xls',sep=''))
      write.table(top10_up,file=paste(test_group,'_vs_',ctrl_group,'_top10_up_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

      top50_up <<- deg_up %>% 
      group_by(cluster) %>% 
      top_n(n = 50, 
              wt = avg_log2FC)
      cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top50_up_deg','.xls',sep=''))
      write.table(top50_up,file=paste(test_group,'_vs_',ctrl_group,'_top50_up_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

      top100_up <<- deg_up %>% 
      group_by(cluster) %>% 
      top_n(n = 100, 
              wt = avg_log2FC)
      cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top100_up_deg','.xls',sep=''))
      write.table(top100_up,file=paste(test_group,'_vs_',ctrl_group,'_top100_up_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)



      diff_genes <- diff_genes %>%
              dplyr::arrange(cluster, avg_log2FC, p_val)
      deg <- subset(diff_genes,avg_log2FC<=-opt$logfc)
      deg <- subset(deg,pct.1>=0.2|pct.2>=0.2)
      deg <- subset(deg,p_val_adj<=0.05)

      cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_filter_down_deg','.xls',sep=''))
      write.table(deg,file=paste(test_group,'_vs_',ctrl_group,'_filter_down_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

      # deg_down <- subset(diff_genes,avg_log2FC<0)
      deg_down <- deg
      top3_down <<- deg_down %>% 
      group_by(cluster) %>% 
      top_n(n = 3, 
              wt = -avg_log2FC)
      cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top3_down_deg','.xls',sep=''))
      write.table(top3_down,file=paste(test_group,'_vs_',ctrl_group,'_top3_down_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

      top5_down <<- deg_down %>% 
      group_by(cluster) %>% 
      top_n(n = 5, 
              wt = -avg_log2FC)
      cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top5_down_deg','.xls',sep=''))
      write.table(top5_down,file=paste(test_group,'_vs_',ctrl_group,'_top5_down_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

      top10_down <<- deg_down %>% 
      group_by(cluster) %>% 
      top_n(n = 10, 
              wt = -avg_log2FC)
      cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top10_down_deg','.xls',sep=''))
      write.table(top10_down,file=paste(test_group,'_vs_',ctrl_group,'_top10_down_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

      top50_down <<- deg_down %>% 
      group_by(cluster) %>% 
      top_n(n = 50, 
              wt = -avg_log2FC)
      cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top50_down_deg','.xls',sep=''))
      write.table(top50_down,file=paste(test_group,'_vs_',ctrl_group,'_top50_down_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

      top100_down <<- deg_down %>% 
      group_by(cluster) %>% 
      top_n(n = 100, 
              wt = -avg_log2FC)
      cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top100_down_deg','.xls',sep=''))
      write.table(top100_down,file=paste(test_group,'_vs_',ctrl_group,'_top100_down_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

      # # pmap_dfr(list(maxIdents,test_group,ctrl_group), plot_scatter_by_group)
      # expscatter_dir = paste("","expscatter", sep="")
      # if(!file.exists(expscatter_dir)){
      #   dir.create(expscatter_dir)
      # }

      # # showgene_dir = paste("","top5_geneplot", sep="")
      # # if(!file.exists(showgene_dir)){
      # #   dir.create(showgene_dir)
      # # }

      # future_mapply(plot_scatter_by_group,maxIdents,test_group,ctrl_group,future.seed = TRUE)
      # # plot_scatter_by_group(0,test_group,ctrl_group)
      setwd(deg_dir3)
    }

  # future_lapply(1:length(compares), compare_group)

  for(gid in 1:length(compares)){
    tryCatch(
        {
          compare_group(gid)
        },
        error=function(cond) {
          ret_code <<- -1
        },
        finally={
          print(paste("done compare: ",compares[gid],sep=""))
        })
  }

  ## 保存RDS

  # saveRDS(sc.combined, file = "./genes_analysis.rds")

  ########################################################################

  ## 按照细胞类型做差异分析

  if(opt$skipanno == "false"){

    ## 鉴定差异基因，按照大的细胞类型分别做
    Idents(sc.combined) <- sc.idents.bak
    clusternum <- length(unique(factor(Idents(sc.combined))))
    DefaultAssay(sc.combined) <- "RNA"
    sc.idents.bak <- Idents(sc.combined)
    # Idents(sc.combined) <- sc.idents.bak

    sc.combined$cellclass.group <- paste(sc.combined$customclassif, sc.combined$group, sep = "_")
    sc.combined$cellclass <- sc.combined$customclassif
    Idents(sc.combined) <- "cellclass.group"
    maxIdents<- unique(factor(sc.combined$customclassif))

    compares <- strsplit(opt$comparelist,",")[[1]]

    deg_dir2 = paste(workdir,"/","deg_bycelltype", sep="")
    if(!file.exists(deg_dir2)){
      dir.create(deg_dir2)
    }
    setwd(deg_dir2)


    # workdir <- getwd()
    # setwd(workdir)

    # source("/users/chengc/work/pipeline/ablifepipe/src/SingleCell/allsample/sc_rbp_visualizations_function.r")


    get_diff <- function(cid,test_group,ctrl_group){
      tryCatch(
          {
            FindMarkers(sc.combined,
                          ident.1 = paste(cid, test_group, sep = "_"),
                          ident.2 = paste(cid, ctrl_group, sep = "_"),
                          min.cells.group = 3,
                          logfc.threshold = 0,
                          min.pct = 0,
                          only.pos = FALSE) %>%
            rownames_to_column(var = "gene") %>%
            cbind(cluster = cid, .)
          },
          error=function(cond) {
            ret_code <<- -1
            data.frame()
          },
          finally={
          })
      
      }

    plot_scatter_by_group <- function(cid,test_group,ctrl_group,deg_up,deg_down){
        # cluster <- "CD4 Naive T"
        Idents(sc.combined) <- sc.idents.bak
        library(ggplot2)
        library(cowplot)
        theme_set(theme_cowplot())

        sub.cells <- subset(sc.combined, customclassif == cid)
        # seurat.subset <- subset(x = sub.cells, subset = (group == test_group & group == ctrl_group) & (another_condition == "Ambigious"))
        sub.cells <- subset(sub.cells,group == c(test_group,ctrl_group))
        # sub.cells <- sub.cells@meta.data
        # table(sub.cells$group)

        nowdir = getwd()

        degs = rbind(deg_up,deg_down)

        expscatter_dir = paste("","expscatter", sep="")
        if(!file.exists(expscatter_dir)){
          dir.create(expscatter_dir)
        }

        setwd(expscatter_dir)

        tryCatch(
          {
            Idents(sub.cells) <- "group"
            avg.sub.cells <- log1p(AverageExpression(sub.cells, verbose = FALSE)$RNA)
            avg.sub.cells.data=data.frame(avg.sub.cells)
            avg.sub.cells.data$gene <- rownames(avg.sub.cells.data)
            diffgenes  = subset(degs, cluster == cid)
            genes.to.label = unique(diffgenes$gene)
            print(genes.to.label)
            p1 <- ggplot(avg.sub.cells.data,aes_string(x=test_group,y=ctrl_group)) + 
                geom_point(
            aes(
              fill = ifelse(gene %in% genes.to.label, "blue", "black")
              ),
              shape = 21,
              alpha = 0.6,
            ) +

            scale_fill_manual(values = c("black", "red"))+

            ggtitle(paste("C",cid,sep=""))
            p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)+NoLegend()
            ggsave(file=paste("C",cid,"_deg_expscatter",".pdf",sep=""),width = 100,height = 100,units = "mm")
            ggsave(file=paste("C",cid,"_deg_expscatter",".png",sep=""),width = 100,height = 100,units = "mm", dpi=300)
          },
          error=function(cond) {
            ret_code <<- -1
          },
          finally={
            print(paste("C",cid,sep=""))
          })
          setwd(nowdir)
    }


    compare_group <- function(gid){
        
        groups <- strsplit(compares[gid],":")[[1]]
        test_group = groups[1]
        ctrl_group = groups[2]

        deg_dir = paste(deg_dir2,"/",test_group,"_vs_", ctrl_group,"_deg", sep="")

        if(!file.exists(deg_dir)){
          dir.create(deg_dir)
        }

        setwd(deg_dir)

        # https://blog.csdn.net/qq_18055167/article/details/104437236
        # diff_genes <- map_dfr(c(maxIdents,test_group,ctrl_group), get_diff)
        Idents(sc.combined) <- "cellclass.group"
        diff_genes <- pmap_dfr(list(maxIdents,test_group,ctrl_group), get_diff)
        # diff_genes <- future_mapply(get_diff,maxIdents,test_group,ctrl_group,future.seed = TRUE)

        diff_genes <- diff_genes[ , c(1:2,4:6, 3, 7)]
        diff_genes <- diff_genes %>%
                dplyr::arrange(cluster, p_val)

        cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_all_genes_export','.xls',sep=''))
        write.table(diff_genes,file=paste(test_group,'_vs_',ctrl_group,'_all_genes_export','.xls',sep=''),append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

        # saveRDS(diff_genes, file = "./diff_genes.rds")
        # saveRDS(conserved_markers, file = "./conserved_markers.rds")
        # saveRDS(combined_markers, file = "./combined_markers.rds")
        # 提取每个类群排名靠前的3个或10个保守标记物
        #上调

        
        diff_genes <- diff_genes %>%
                dplyr::arrange(cluster, -avg_log2FC, p_val)
        deg <- subset(diff_genes,avg_log2FC>=opt$logfc)
        deg <- subset(deg,pct.1>=0.2|pct.2>=0.2)
        deg <- subset(deg,p_val_adj<=0.05)

        cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_filter_up_deg','.xls',sep=''))
        write.table(deg,file=paste(test_group,'_vs_',ctrl_group,'_filter_up_deg','.xls',sep=''),
                    append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

        # deg_up <- subset(diff_genes,avg_log2FC>=0)
        deg_up <- deg
        top3_up <<- deg_up %>% 
        group_by(cluster) %>% 
        top_n(n = 3, 
                wt = avg_log2FC)
        cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top3_up_deg','.xls',sep=''))
        write.table(top3_up,file=paste(test_group,'_vs_',ctrl_group,'_top3_up_deg','.xls',sep=''),
                    append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)
        
        top5_up <<- deg_up %>% 
        group_by(cluster) %>% 
        top_n(n = 5, 
                wt = avg_log2FC)
        cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top5_up_deg','.xls',sep=''))
        write.table(top5_up,file=paste(test_group,'_vs_',ctrl_group,'_top5_up_deg','.xls',sep=''),
                    append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

        top10_up <<- deg_up %>% 
        group_by(cluster) %>% 
        top_n(n = 10, 
                wt = avg_log2FC)
        cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top10_up_deg','.xls',sep=''))
        write.table(top10_up,file=paste(test_group,'_vs_',ctrl_group,'_top10_up_deg','.xls',sep=''),
                    append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

        top50_up <<- deg_up %>% 
        group_by(cluster) %>% 
        top_n(n = 50, 
                wt = avg_log2FC)
        cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top50_up_deg','.xls',sep=''))
        write.table(top50_up,file=paste(test_group,'_vs_',ctrl_group,'_top50_up_deg','.xls',sep=''),
                    append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

        top100_up <<- deg_up %>% 
        group_by(cluster) %>% 
        top_n(n = 100, 
                wt = avg_log2FC)
        cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top100_up_deg','.xls',sep=''))
        write.table(top100_up,file=paste(test_group,'_vs_',ctrl_group,'_top100_up_deg','.xls',sep=''),
                    append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)



        diff_genes <- diff_genes %>%
                dplyr::arrange(cluster, avg_log2FC, p_val)
        deg <- subset(diff_genes,avg_log2FC<=-opt$logfc)
        deg <- subset(deg,pct.1>=0.2|pct.2>=0.2)
        deg <- subset(deg,p_val_adj<=0.05)

        cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_filter_down_deg','.xls',sep=''))
        write.table(deg,file=paste(test_group,'_vs_',ctrl_group,'_filter_down_deg','.xls',sep=''),
                    append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

        # deg_down <- subset(diff_genes,avg_log2FC<0)
        deg_down <- deg
        top3_down <<- deg_down %>% 
        group_by(cluster) %>% 
        top_n(n = 3, 
                wt = -avg_log2FC)
        cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top3_down_deg','.xls',sep=''))
        write.table(top3_down,file=paste(test_group,'_vs_',ctrl_group,'_top3_down_deg','.xls',sep=''),
                    append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

        top5_down <<- deg_down %>% 
        group_by(cluster) %>% 
        top_n(n = 5, 
                wt = -avg_log2FC)
        cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top5_down_deg','.xls',sep=''))
        write.table(top5_down,file=paste(test_group,'_vs_',ctrl_group,'_top5_down_deg','.xls',sep=''),
                    append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

        top10_down <<- deg_down %>% 
        group_by(cluster) %>% 
        top_n(n = 10, 
                wt = -avg_log2FC)
        cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top10_down_deg','.xls',sep=''))
        write.table(top10_down,file=paste(test_group,'_vs_',ctrl_group,'_top10_down_deg','.xls',sep=''),
                    append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

        top50_down <<- deg_down %>% 
        group_by(cluster) %>% 
        top_n(n = 50, 
                wt = -avg_log2FC)
        cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top50_down_deg','.xls',sep=''))
        write.table(top50_down,file=paste(test_group,'_vs_',ctrl_group,'_top50_down_deg','.xls',sep=''),
                    append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)

        top100_down <<- deg_down %>% 
        group_by(cluster) %>% 
        top_n(n = 100, 
                wt = -avg_log2FC)
        cat(c(''),file=paste(test_group,'_vs_',ctrl_group,'_top100_down_deg','.xls',sep=''))
        write.table(top100_down,file=paste(test_group,'_vs_',ctrl_group,'_top100_down_deg','.xls',sep=''),
                    append=T,quote=F,sep='\t',row.names = FALSE,col.names = TRUE)


        # # pmap_dfr(list(maxIdents,test_group,ctrl_group), plot_scatter_by_group)
        # expscatter_dir = paste("","expscatter", sep="")
        # if(!file.exists(expscatter_dir)){
        #   dir.create(expscatter_dir)
        # }

        # # showgene_dir = paste("","deg_geneplot", sep="")
        # # if(!file.exists(showgene_dir)){
        # #   dir.create(showgene_dir)
        # # }

        # # future_mapply(plot_scatter_by_group,maxIdents,test_group,ctrl_group,deg_up,deg_down,future.seed = TRUE)

        # for(cid in maxIdents){
        #   tryCatch(
        #       {
        #         plot_scatter_by_group(cid,test_group,ctrl_group,deg_up,deg_down)
        #       },
        #       error=function(cond) {
        #         ret_code <<- -1
        #       },
        #       finally={
        #         print(paste("done exp scatter: ",cid,sep=""))
        #       })
        # }
        
        # setwd(deg_dir)

        # # plot_scatter_by_group(0,test_group,ctrl_group)
        setwd(deg_dir2)
      }

    # future_lapply(1:length(compares), compare_group)

    for(gid in 1:length(compares)){
      tryCatch(
          {
            compare_group(gid)
          },
          error=function(cond) {
            ret_code <<- -1
          },
          finally={
            print(paste("done compare: ",compares[gid],sep=""))
          })
    }

  }

  ########################################################################

  ## 使用所有细胞做差异分析,by all

  setwd(workdir)

  deg_dir = paste(workdir,"/","deg_byall", sep="")
  if(!file.exists(deg_dir)){
    dir.create(deg_dir)
  }
  setwd(deg_dir)
  sc.idents.bak <- Idents(sc.combined)


  get_diff_all <- function(test_group,ctrl_group){
    tryCatch(
        {
          FindMarkers(sc.combined,
                        ident.1 = test_group,
                        ident.2 = ctrl_group,
                        min.cells.group = 3,
                        logfc.threshold = 0,
                        min.pct = 0,
                        only.pos = FALSE)
        },
        error=function(cond) {
          ret_code <<- -1
          data.frame()
        },
        finally={
        })
    
    }

  plot_scatter_by_group <- function(test_group,ctrl_group){
      # cluster <- "CD4 Naive T"
      Idents(sc.combined) <- sc.idents.bak
      # Idents(sc.combined) <- "RBP"
      library(ggplot2)
      library(cowplot)
      theme_set(theme_cowplot())

      # sub.cells <- subset(sc.combined, idents = cid)
      # seurat.subset <- subset(x = sub.cells, subset = (group == test_group & group == ctrl_group) & (another_condition == "Ambigious"))
      # sub.cells <- subset(sc.combined,group == c(test_group,ctrl_group))
      sub.cells <- subset(sc.combined,group == c(test_group,ctrl_group))
      # sub.cells <- sub.cells@meta.data
      # table(sub.cells$group)

      degdir = getwd()

      top10 = rbind(top10_up,top10_down)

      expscatter_dir = paste("","expscatter", sep="")
      if(!file.exists(expscatter_dir)){
        dir.create(expscatter_dir)
      }

      setwd(expscatter_dir)

      tryCatch(
        {
          Idents(sub.cells) <- "group"
          avg.sub.cells <- log1p(AverageExpression(sub.cells, verbose = FALSE)$RNA)
          avg.sub.cells.data=data.frame(avg.sub.cells)
          avg.sub.cells.data$gene <- rownames(avg.sub.cells.data)
          # diffgenes  = subset(top10, cluster == cid)
          diffgenes  = top10
          genes.to.label = unique(rownames(diffgenes))
          print(genes.to.label)
          p1 <- ggplot(avg.sub.cells.data,aes_string(x=test_group,y=ctrl_group)) + 
              geom_point(
          aes(
            fill = ifelse(gene %in% genes.to.label, "blue", "black")
            ),
            shape = 21,
            alpha = 0.6,
          ) +

          scale_fill_manual(values = c("black", "red"))+

          ggtitle(paste("top10",sep=""))
          p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, max.overlaps=Inf, xnudge=0, ynudge=0)+NoLegend()
          ggsave(file=paste("top10_expscatter",".pdf",sep=""),width = 100,height = 100,units = "mm")
          ggsave(file=paste("top10_expscatter",".png",sep=""),width = 100,height = 100,units = "mm", dpi=300)
        },
        error=function(cond) {
          ret_code <<- -1
        },
        finally={
          print(paste("top10",sep=""))
        })
      

      setwd(degdir)

  }

  compare_group <- function(gid){
      
      groups <- strsplit(compares[gid],":")[[1]]
      test_group = groups[1]
      ctrl_group = groups[2]

      deg_dir2 = paste(deg_dir,"/",test_group,"_vs_", ctrl_group,"_deg", sep="")

      if(!file.exists(deg_dir2)){
        dir.create(deg_dir2)
      }

      setwd(deg_dir2)

      # https://blog.csdn.net/qq_18055167/article/details/104437236
      # diff_genes <- map_dfr(c(maxIdents,test_group,ctrl_group), get_diff)
      # Idents(sc.combined) <- "celltype.group"
      # diff_genes <- pmap_dfr(list(maxIdents,test_group,ctrl_group), get_diff)
      
      Idents(sc.combined) <- "group"
      # Idents(sc.combined) <- "grade"
      # Idents(sc.combined) <- "cnv"
      diff_genes<-""
      diff_genes <- get_diff_all(test_group,ctrl_group)
      # diff_genes <- future_mapply(get_diff,maxIdents,test_group,ctrl_group,future.seed = TRUE)

      # diff_genes <- diff_genes[ , c(1:2,4:6, 3, 7)]
      diff_genes <- diff_genes %>%
              dplyr::arrange(p_val)

      cat(c('Gene\t'),file=paste(test_group,'_vs_',ctrl_group,'_all_genes_export','.xls',sep=''))
      write.table(diff_genes,file=paste(test_group,'_vs_',ctrl_group,'_all_genes_export','.xls',sep=''),append=T,quote=F,sep='\t',row.names = TRUE,col.names = TRUE)

      # saveRDS(diff_genes, file = "./diff_genes.rds")
      # saveRDS(conserved_markers, file = "./conserved_markers.rds")
      # saveRDS(combined_markers, file = "./combined_markers.rds")
      # 提取每个类群排名靠前的3个或10个保守标记物
      #上调

      
      diff_genes <- diff_genes %>%
              dplyr::arrange(-avg_log2FC, p_val)
      deg <- subset(diff_genes,avg_log2FC>=opt$logfc)
      deg <- subset(deg,pct.1>=0.2|pct.2>=0.2)
      deg <- subset(deg,p_val_adj<=0.05)

      cat(c('Gene\t'),file=paste(test_group,'_vs_',ctrl_group,'_filter_up_deg','.xls',sep=''))
      write.table(deg,file=paste(test_group,'_vs_',ctrl_group,'_filter_up_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = TRUE,col.names = TRUE)

      # deg_up <- subset(diff_genes,avg_log2FC>=0)
      deg_up <- deg
      top3_up <<- deg_up %>% 
      top_n(n = 3, 
              wt = avg_log2FC)
      cat(c('Gene\t'),file=paste(test_group,'_vs_',ctrl_group,'_top3_up_deg','.xls',sep=''))
      write.table(top3_up,file=paste(test_group,'_vs_',ctrl_group,'_top3_up_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = TRUE,col.names = TRUE)
      
      top5_up <<- deg_up %>% 
      top_n(n = 5, 
              wt = avg_log2FC)
      cat(c('Gene\t'),file=paste(test_group,'_vs_',ctrl_group,'_top5_up_deg','.xls',sep=''))
      write.table(top5_up,file=paste(test_group,'_vs_',ctrl_group,'_top5_up_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = TRUE,col.names = TRUE)

      top10_up <<- deg_up %>% 
      top_n(n = 10, 
              wt = avg_log2FC)
      cat(c('Gene\t'),file=paste(test_group,'_vs_',ctrl_group,'_top10_up_deg','.xls',sep=''))
      write.table(top10_up,file=paste(test_group,'_vs_',ctrl_group,'_top10_up_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = TRUE,col.names = TRUE)

      top50_up <<- deg_up %>% 
      top_n(n = 50, 
              wt = avg_log2FC)
      cat(c('Gene\t'),file=paste(test_group,'_vs_',ctrl_group,'_top50_up_deg','.xls',sep=''))
      write.table(top50_up,file=paste(test_group,'_vs_',ctrl_group,'_top50_up_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = TRUE,col.names = TRUE)

      top100_up <<- deg_up %>% 
      top_n(n = 100, 
              wt = avg_log2FC)
      cat(c('Gene\t'),file=paste(test_group,'_vs_',ctrl_group,'_top100_up_deg','.xls',sep=''))
      write.table(top100_up,file=paste(test_group,'_vs_',ctrl_group,'_top100_up_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = TRUE,col.names = TRUE)



      diff_genes <- diff_genes %>%
              dplyr::arrange(avg_log2FC, p_val)
      deg <- subset(diff_genes,avg_log2FC<=-opt$logfc)
      deg <- subset(deg,pct.1>=0.2|pct.2>=0.2)
      deg <- subset(deg,p_val_adj<=0.05)

      cat(c('Gene\t'),file=paste(test_group,'_vs_',ctrl_group,'_filter_down_deg','.xls',sep=''))
      write.table(deg,file=paste(test_group,'_vs_',ctrl_group,'_filter_down_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = TRUE,col.names = TRUE)

      # deg_down <- subset(diff_genes,avg_log2FC<0)
      deg_down <- deg
      top3_down <<- deg_down %>% 
      top_n(n = 3, 
              wt = -avg_log2FC)
      cat(c('Gene\t'),file=paste(test_group,'_vs_',ctrl_group,'_top3_down_deg','.xls',sep=''))
      write.table(top3_down,file=paste(test_group,'_vs_',ctrl_group,'_top3_down_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = TRUE,col.names = TRUE)

      top5_down <<- deg_down %>% 
      top_n(n = 5, 
              wt = -avg_log2FC)
      cat(c('Gene\t'),file=paste(test_group,'_vs_',ctrl_group,'_top5_down_deg','.xls',sep=''))
      write.table(top5_down,file=paste(test_group,'_vs_',ctrl_group,'_top5_down_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = TRUE,col.names = TRUE)

      top10_down <<- deg_down %>% 
      top_n(n = 10, 
              wt = -avg_log2FC)
      cat(c('Gene\t'),file=paste(test_group,'_vs_',ctrl_group,'_top10_down_deg','.xls',sep=''))
      write.table(top10_down,file=paste(test_group,'_vs_',ctrl_group,'_top10_down_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = TRUE,col.names = TRUE)

      top50_down <<- deg_down %>% 
      top_n(n = 50, 
              wt = -avg_log2FC)
      cat(c('Gene\t'),file=paste(test_group,'_vs_',ctrl_group,'_top50_down_deg','.xls',sep=''))
      write.table(top50_down,file=paste(test_group,'_vs_',ctrl_group,'_top50_down_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = TRUE,col.names = TRUE)

      top100_down <<- deg_down %>% 
      top_n(n = 100, 
              wt = -avg_log2FC)
      cat(c('Gene\t'),file=paste(test_group,'_vs_',ctrl_group,'_top100_down_deg','.xls',sep=''))
      write.table(top100_down,file=paste(test_group,'_vs_',ctrl_group,'_top100_down_deg','.xls',sep=''),
                  append=T,quote=F,sep='\t',row.names = TRUE,col.names = TRUE)

      # # pmap_dfr(list(maxIdents,test_group,ctrl_group), plot_scatter_by_group)
      # expscatter_dir = paste("","expscatter", sep="")
      # if(!file.exists(expscatter_dir)){
      # 	dir.create(expscatter_dir)
      # }

      # # showgene_dir = paste("","top5_geneplot", sep="")
      # # if(!file.exists(showgene_dir)){
      # # 	dir.create(showgene_dir)
      # # }

      # plot_scatter_by_group(test_group,ctrl_group)
      # # plot_scatter_by_group(0,test_group,ctrl_group)
      setwd("../")
    }


  # opt$comparelist = "HighGrade:LowGrade"
  compares <- strsplit(opt$comparelist,",")[[1]]
  Idents(sc.combined) <- "group"

  # future_lapply(1:length(compares), compare_group)

  for(gid in 1:length(compares)){
    tryCatch(
        {
          compare_group(gid)
        },
        error=function(cond) {
          ret_code <<- -1
        },
        finally={
          print(paste("done compare: ",compares[gid],sep=""))
        })
  }
    
}else if(opt$mmethod == "presto"){

  #差异分析的函数
  presto_de <- function(ass="RNA",ide="seurat_clusters"){

    DefaultAssay(sc.combined) <- ass
    Idents(sc.combined) <- ide
    sc.combined$group_clusters <- paste0(sc.combined$group,"_",Idents(sc.combined))

    if(opt$comparelist != ""){
    comparison <- list()
    compare <- strsplit(opt$comparelist,",")[[1]]
    for(i in 1:length(compare)){
      comparison[i] <- list(c(strsplit(compare[i],":")[[1]]))
      }
    }

    #统计
    deg_stat <- data.frame()

    #差异分析
    for(i in 1:length(comparison)){
      #创建目录
      temp_outdir <- paste0(deg_dir,"/",comparison[[i]][1],"_vs_",comparison[[i]][2],"_deg")
      if(!file.exists(temp_outdir)){
        dir.create(temp_outdir)
      }
      setwd(temp_outdir) 

      #使用FindAllMarkers非常慢，即便是开了多线程也非常慢
      #使用presto包做极快的Wilcoxon测试
      all_markers <- list()

      cellgroup <- Idents(sc.combined) %>% as.character() %>% unique() %>% mixedsort()

      for(j in cellgroup){
        
        #必须保证两组间都有细胞
        n1 <- sum(sc.combined$group_clusters == paste0(comparison[[i]][1],"_",j))
        n2 <- sum(sc.combined$group_clusters == paste0(comparison[[i]][2],"_",j))
        
        if(n1 > 1 && n2 > 1){
          
          #由于列表的索引是1开头的
          all_markers[[j]]<- wilcoxauc(
            sc.combined,
            group_by = "group_clusters",
            assay = "data",
            groups_use = c(paste0(comparison[[i]][1],"_",j),paste0(comparison[[i]][2],"_",j)),
            seurat_assay = ass
          ) %>% dplyr::select(group,feature,logFC,pct_in,pct_out,pval,padj) 
          
          all_markers[[j]] <- dplyr::filter(all_markers[[j]], group==paste0(comparison[[i]][1],"_",j))
          all_markers[[j]]$group <- j 
          
        }
        
      }

      all_markers <- do.call(rbind,all_markers) #合并所有的结果成为一个数据框

      all_markers <- dplyr::mutate(all_markers,avg_log2FC = log(exp(logFC),base = 2), pct.1 = round(pct_in/100,3), pct.2 = round(pct_out/100,3)) %>% dplyr::select(group,feature,avg_log2FC,pct.1,pct.2,pval,padj) 

      colnames(all_markers) <- c("cluster","gene","avg_log2FC","pct.1","pct.2","p_val","p_val_adj") 
      
      all_markers <- all_markers[mixedorder(all_markers$cluster),] #按照自然排序

      write.table(all_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_all_genes_export.xls"),sep = "\t",quote = F,row.names = F,col.names = T)

      #####################################################################################################
      #down差异
      filter_down_markers <- dplyr::filter(all_markers,avg_log2FC < -opt$logfc & p_val_adj < 0.05 & (pct.1 > 0.2 | pct.2 > 0.2)) %>% dplyr::arrange(cluster,avg_log2FC,-pct.1,p_val)
      filter_down_markers <- filter_down_markers[mixedorder(filter_down_markers$cluster),] #按照自然排序
      write.table(filter_down_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_filter_down_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)

      #提取每个类群排名靠前的top marker
      top3_markers <- filter_down_markers %>% group_by(cluster) %>% top_n(n = -3, wt = avg_log2FC)
      write.table(top3_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top3_down_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)

      top10_markers <- filter_down_markers %>% group_by(cluster) %>% top_n(n = -10, wt = avg_log2FC)
      write.table(top10_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top10_down_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)

      top50_markers <- filter_down_markers %>% group_by(cluster) %>% top_n(n = -50, wt = avg_log2FC)
      write.table(top50_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top50_down_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)    

      top100_markers <- filter_down_markers %>% group_by(cluster) %>% top_n(n = -100, wt = avg_log2FC)
      write.table(top100_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top100_down_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)

      tryCatch(
      {
        #差异marker数量统计
        deg_stat_tmp <- filter_down_markers %>% 
        group_by(cluster) %>% 
        summarise(count = n()) %>%
        mutate(com = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_down")) %>%
        as.data.frame()  
        deg_stat <- rbind(deg_stat,deg_stat_tmp)      
      },
      error = function(e){message("差异marker数量统计error")}
      )

      #####################################################################################################
      #up差异
      filter_up_markers <- dplyr::filter(all_markers,avg_log2FC > opt$logfc & p_val_adj < 0.05 & (pct.1 > 0.2 | pct.2 > 0.2)) %>% dplyr::arrange(cluster,-avg_log2FC,-pct.1,p_val)
      filter_up_markers <- filter_up_markers[mixedorder(filter_up_markers$cluster),] #按照自然排序
      write.table(filter_up_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_filter_up_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)

      #提取每个类群排名靠前的top marker
      top3_markers <- filter_up_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
      write.table(top3_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top3_up_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)

      top10_markers <- filter_up_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
      write.table(top10_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top10_up_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)

      top50_markers <- filter_up_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
      write.table(top50_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top50_up_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)    

      top100_markers <- filter_up_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
      write.table(top100_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top100_up_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)

      tryCatch(
      {
        #差异marker数量统计
        deg_stat_tmp <- filter_up_markers %>% 
        group_by(cluster) %>% 
        summarise(count = n()) %>%
        mutate(com = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_up")) %>%
        as.data.frame()  
        deg_stat <- rbind(deg_stat,deg_stat_tmp) %>% unique()   
      },
      error = function(e){message("差异marker数量统计error")}
      )
      
    }

    #统计表格
    setwd(deg_dir)
    #长转宽
    deg_stat_wide <- deg_stat %>% pivot_wider(names_from = com, values_from = count) %>% as.data.frame()
    deg_stat_wide[is.na(deg_stat_wide)] <- 0 
    deg_stat_wide <- deg_stat_wide[mixedorder(deg_stat_wide$cluster),]
    write.table(deg_stat_wide,file = "deg_stat.xls",sep = "\t",quote = F,row.names = F,col.names = T)

  }

  ####################deg by seurat clusters#############################
  deg_dir = paste(workdir,"/","deg_bycluster", sep="")
  if(!file.exists(deg_dir)){
    dir.create(deg_dir)
  }
  setwd(deg_dir)  

  presto_de("RNA","seurat_clusters")

  ####################deg by celltype#############################
  if(opt$skipanno == "false"){
    deg_dir = paste(workdir,"/","deg_bycelltype", sep="")
    if(!file.exists(deg_dir)){
      dir.create(deg_dir)
    }
    setwd(deg_dir)  

    presto_de("RNA","customclassif") 
  }

  ####################deg by all#############################
  deg_dir = paste(workdir,"/","deg_byall", sep="")
  if(!file.exists(deg_dir)){
    dir.create(deg_dir)
  }
  setwd(deg_dir)  

  DefaultAssay(sc.combined) <- "RNA"
  Idents(sc.combined) <- "group"

  if(opt$comparelist != ""){
    comparison <- list()
    compare <- strsplit(opt$comparelist,",")[[1]]
    for(i in 1:length(compare)){
      comparison[i] <- list(c(strsplit(compare[i],":")[[1]]))
    }
  }

  #统计
  deg_stat <- data.frame()

  #差异分析
  for(i in 1:length(comparison)){
    #创建目录
    temp_outdir <- paste0(deg_dir,"/",comparison[[i]][1],"_vs_",comparison[[i]][2],"_deg")
    if(!file.exists(temp_outdir)){
      dir.create(temp_outdir)
    }
    setwd(temp_outdir) 
    
    #必须保证两组间都有细胞
    n1 <- sum(sc.combined$group == comparison[[i]][1])
    n2 <- sum(sc.combined$group == comparison[[i]][2])
    
    if(n1 > 5 && n2 > 5){
      
      #由于列表的索引是1开头的
      all_markers <- wilcoxauc(
        sc.combined,
        group_by = "group",
        assay = "data",
        groups_use = c(comparison[[i]][1],comparison[[i]][2]),
        seurat_assay = "RNA"
      ) %>% dplyr::filter(group==comparison[[i]][1]) %>%
        dplyr::select(feature,logFC,pct_in,pct_out,pval,padj)
      
    }
    
    all_markers <- dplyr::mutate(all_markers,avg_log2FC = log(exp(logFC),base = 2), pct.1 = round(pct_in/100,3), pct.2 = round(pct_out/100,3)) %>% dplyr::select(feature,avg_log2FC,pct.1,pct.2,pval,padj) 
    
    colnames(all_markers) <- c("gene","avg_log2FC","pct.1","pct.2","p_val","p_val_adj") 
    
    write.table(all_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_all_genes_export.xls"),sep = "\t",quote = F,row.names = F,col.names = T)
  
    #####################################################################################################
    #down差异
    filter_down_markers <- dplyr::filter(all_markers,avg_log2FC < -opt$logfc & p_val_adj < 0.05 & (pct.1 > 0.2 | pct.2 > 0.2)) %>% dplyr::arrange(avg_log2FC,-pct.1,p_val)
    write.table(filter_down_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_filter_down_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)
    
    #提取每个类群排名靠前的top marker
    top3_markers <- head(filter_down_markers,n = 3)
    write.table(top3_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top3_down_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)
    
    top10_markers <- head(filter_down_markers,n = 10)
    write.table(top10_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top10_down_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)
    
    top50_markers <- head(filter_down_markers,n = 50)
    write.table(top50_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top50_down_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)    
    
    top100_markers <- head(filter_down_markers,n = 100)
    write.table(top100_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top100_down_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)
    
    tryCatch(
      {
        #差异marker数量统计
        deg_stat_tmp <- filter_down_markers %>% 
          summarise(count = n()) %>%
          mutate(group = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2])) %>%
          mutate(com = "down") %>%
          as.data.frame()  
        deg_stat <- rbind(deg_stat,deg_stat_tmp)      
      },
      error = function(e){message("差异marker数量统计error")}
    )
    
    #####################################################################################################
    #up差异
    filter_up_markers <- dplyr::filter(all_markers,avg_log2FC > opt$logfc & p_val_adj < 0.05 & (pct.1 > 0.2 | pct.2 > 0.2)) %>% dplyr::arrange(-avg_log2FC,-pct.1,p_val)
    write.table(filter_up_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_filter_up_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)
    
    #提取每个类群排名靠前的top marker
    top3_markers <- head(filter_up_markers,n = 3)
    write.table(top3_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top3_up_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)
    
    top10_markers <- head(filter_up_markers,n = 10)
    write.table(top10_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top10_up_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)
    
    top50_markers <- head(filter_up_markers,n = 50)
    write.table(top50_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top50_up_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)    
    
    top100_markers <- head(filter_up_markers,n = 100)
    write.table(top100_markers,file = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2],"_top100_up_deg.xls"),sep = "\t",quote = F,row.names = F,col.names = T)
    
    tryCatch(
      {
        #差异marker数量统计
        deg_stat_tmp <- filter_up_markers %>% 
          summarise(count = n()) %>%
          mutate(group = paste0(comparison[[i]][1],"_vs_",comparison[[i]][2])) %>%
          mutate(com = "up") %>%
          as.data.frame()  
        deg_stat <- rbind(deg_stat,deg_stat_tmp) %>% unique()    
      },
      error = function(e){message("差异marker数量统计error")}
    ) 
    
  }

  #统计表格
  setwd(deg_dir)
  #长转宽
  deg_stat_wide <- deg_stat %>% pivot_wider(names_from = com, values_from = count) %>% as.data.frame()
  deg_stat_wide[is.na(deg_stat_wide)] <- 0 
  write.table(deg_stat_wide,file = "deg_stat.xls",sep = "\t",quote = F,row.names = F,col.names = T)



}


