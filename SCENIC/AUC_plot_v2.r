rm(list = ls())
#setwd("D:/Workspace/R_learning/sc_analysis/SCENIC/RX-7710652-N_Polydactyly/")

#传参
library(optparse)
option_list <- list(
  make_option(c("-l", "--loom"),action = "store",type = "character",
              help = "The Input loom file"),
  make_option(c("-i", "--rds"),action = "store",type = "character",
              help = "The Input rds file"),
  make_option(c("-g", "--genelist"),action = "store",type = "character",
              help = "The Input rds file"),
  make_option(c("-d", "--idents"),action = "store",type = "character",default = "seurat_clusters",
              help = "The active idents"),              
  make_option(c("-s", "--splitby"),action = "store",type = "character",default = "",
              help = "The split by"),              
  make_option(c("-o", "--outdir"),action = "store",type = "character",default = "./",
              help = "The outdir,the default is './'")
)

opt <- parse_args(OptionParser(option_list = option_list)) ####解释器

##加载包
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(stats)
library(ggplot2)
#library(ggheatmap)
library(reshape2)
library(ggsci)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(limma)
library(stringr)
library(gtools)

#创建文件夹
if(!file.exists(opt$outdir)){
  dir.create(opt$outdir)
}
setwd(opt$outdir) 
workdir <- getwd()

#test
# opt$loom <- "/data1/projects/single-cell/RX-7710668-N_SAR/figure/20230823/Fig2/tf/auc_mtx.loom"
# opt$rds <- "/data1/projects/single-cell/RX-7710668-N_SAR/figure/20230823/Fig2/tf/cells_analysis_sub.rds"
# opt$genelist <- "/users/wuyf/public/Singlecell/pyscenic/v2/bycluster/3_Regulons_plot/AUC_plot_RSS_TF/RSS_TF"

#
loom <- open_loom(opt$loom) 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
#regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
#tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

#查看富集到的全部的TF
#读入单细胞数据
#opt$rds <- "filter.bycells.rds"
sc <- readRDS(opt$rds)
#将每一个细胞的每一个TF的表达量与单细胞的每一个细胞匹配
sub_regulonAUC <- regulonAUC[,match(colnames(sc),colnames(regulonAUC))]
dim(sub_regulonAUC)
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(sc))
#[1] TRUE  一定要为TRUE

##==============================TF_AUC与seurat结合===========================
#普通展示
next_regulonAUC <- regulonAUC[,match(colnames(sc),colnames(regulonAUC))]
dim(next_regulonAUC)#将AUC结果于seurat对象结合

regulon_AUC <- regulonAUC@NAMES
sc@meta.data = cbind(sc@meta.data ,t(assay(next_regulonAUC[regulon_AUC,])))

#自己选定感兴趣的或者比较重要的转录因子
#opt$genelist <- "keratinocytes_fibroblasts_RSS_TF"
TF_plot <- read.table(opt$genelist,quote = "",header = F)
TF_plot <- as.character(TF_plot$V1)

#dotplot
# p1 <- DotPlot(sc, features = TF_plot,group.by = "customclassif")+
#   theme_bw()+
#   theme(panel.grid = element_blank(), 
#         axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
#   labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
# ggsave("AUC_dotplot.png",p1,width = 8,height = 5)
# ggsave("AUC_dotplot.pdf",p1,width = 8,height = 5)

Idents(sc) <- opt$idents
Idents(sc) <- factor(Idents(sc),levels=mixedsort(levels(Idents(sc))))

#设定图形宽度
mywidth <- length(levels(Idents(sc)))
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

if(opt$splitby != ""){
  splitnum <- length(unique(sc@meta.data[,opt$splitby]))
}

#作图
if(opt$splitby == ""){

  #featureplot
  for(i in TF_plot){
    p1 <- FeaturePlot(sc, features = i,reduction = 'umap')
    ggsave(paste0("AUC_",i,"_featureplot.png"),p1,width = 5,height = 4)
    ggsave(paste0("AUC_",i,"_featureplot.pdf"),p1,width = 5,height = 4)
  }

  #violin plot 小提琴图
  for(i in TF_plot){
    p1 <- VlnPlot(sc, features = i,pt.size = 0) + labs(y = "AUC Score") + NoLegend()
    ggsave(paste0("AUC_",i,"_vlnplot.png"),p1,width = mywidth,height = 4)
    ggsave(paste0("AUC_",i,"_vlnplot.pdf"),p1,width = mywidth,height = 4)
  }

}else{

  #featureplot split by group
  for(i in TF_plot){
    p1 <- FeaturePlot(sc, features = i, reduction = 'umap', split.by = opt$splitby)
    ggsave(paste0("AUC_",i,"_featureplot_split.png"),p1,width = 5*splitnum,height = 4)
    ggsave(paste0("AUC_",i,"_featureplot_split.pdf"),p1,width = 5*splitnum,height = 4)
  }

  #violin plot 小提琴图 split by group
  for(i in TF_plot){
    p1 <- VlnPlot(sc, features = i, pt.size = 0, split.by = opt$splitby) + labs(y = "AUC Score")
    ggsave(paste0("AUC_",i,"_vlnplot_split.png"),p1,width = mywidth,height = 4)
    ggsave(paste0("AUC_",i,"_vlnplot_split.pdf"),p1,width = mywidth,height = 4)
  }

}
