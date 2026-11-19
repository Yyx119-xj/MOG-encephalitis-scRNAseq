# setwd("D:/workspace/R/sc_analysis/pct_plot/")

rm(list = ls())

#加载包
library(ggpubr)
library(dplyr)
library(ggsci)
library(RColorBrewer)
library(patchwork)
# library(Seurat)
library(gtools)
library(optparse)

#参数设置
option_list <- list(
  make_option(c("-i", "--inputfile"),action = "store",type = "character",
              help = "The Input file"),
  make_option(c("-c", "--comparelist"),action = "store",type = "character",default = "",
              help = "The comparelist"),
  make_option(c("-h","--height"),action = "store",type = "numeric",default = "5",
              help = "The height"),
  make_option(c("-w","--width"),action = "store",type = "numeric",default = "12",
              help = "The width"),
  make_option(c("-o", "--outdir"),action = "store",type = "character",default = "./",
              help = "The outdir,the default is './'")
)

opt <- parse_args(OptionParser(option_list = option_list,add_help_option=FALSE)) ####解释器

#创建文件夹
if(!file.exists(opt$outdir)){
  dir.create(opt$outdir)
}
setwd(opt$outdir) 
workdir <- getwd()

##测试
# opt$inputfile <- "stat_celltype_fraction_by_sample.ed.xls"
# opt$comparelist <- "LTBI:HC,TB:HC"

##数据预处理
mydata <- read.table(opt$inputfile,header = T,
                     sep = "\t",check.names = F,com = "",quote = "")

if(opt$comparelist != ""){
  comparison <- list()
  compare <- strsplit(opt$comparelist,",")[[1]]
  for(i in 1:length(compare)){
    comparison[i] <- list(c(strsplit(compare[i],":")[[1]]))
  }
}

#排序并转成因子
mydata$newgroup <- factor(mydata$newgroup,levels = mydata$newgroup %>% unique() %>% mixedsort())

if(opt$comparelist != ""){

  #作图并标注显著性
  ##boxplot
  p1 <- ggboxplot(mydata, x = "samplegroup",y = "Percentage", 
                  fill = "samplegroup",add = "jitter",
                  palette = "npg")
  p2 <- facet(p1, facet.by = "newgroup",nrow = 1,scales = "fixed",strip.position = "bottom",panel.labs.background = list(color = NA)) +
    theme(legend.position="right",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(color = "black"),
          strip.placement = "outside",
          strip.background.x = element_blank(),
          strip.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,color = "black")) + 
    scale_y_continuous(expand = c(0,0)) + # y轴从0開始
    stat_compare_means(
      method = "t.test", #t.test or wilcox.test
      paired = F,hide.ns = F,
      comparisons = comparison,
      label = "p.signif",
      tip.length = c(0,0,0,0),
      bracket.size = 0.5,
      step.increase = 0.05,
      vjust = 1.5,
      size = 4
    )
  ggsave("cell_percentage_boxplot.png",p2,height = opt$height,width = opt$width,dpi = 300)
  ggsave("cell_percentage_boxplot.pdf",p2,height = opt$height,width = opt$width)

  ##barplot
  p3 <- ggbarplot(mydata, x = "samplegroup",y = "Percentage", 
                  fill = "samplegroup",facet.by = "newgroup",
                  add = "mean_se",error.plot = "upper_errorbar",
                  palette = "npg",ggtheme = theme_classic())
  p4 <- facet(p3, facet.by = "newgroup",nrow = 1,scales = "fixed",strip.position = "bottom",panel.labs.background = list(color = NA)) +
    theme(legend.position="right",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(color = "black"),
          strip.placement = "outside",
          strip.background.x = element_blank(),
          strip.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,color = "black")) + 
    scale_y_continuous(expand = c(0,0)) + # y轴从0開始
    stat_compare_means(
      method = "t.test", #t.test or wilcox.test
      paired = F,hide.ns = F,
      comparisons = comparison,
      label = "p.signif",
      tip.length = c(0,0,0,0),
      bracket.size = 0.5,
      step.increase = 0.05,
      vjust = 1.5,
      size = 4
    )
  ggsave("cell_percentage_barplot.png",p4,height = opt$height,width = opt$width,dpi = 300)
  ggsave("cell_percentage_barplot.pdf",p4,height = opt$height,width = opt$width)
}else{
  #作图不标注显著性
  ##boxplot
  p1 <- ggboxplot(mydata, x = "samplegroup",y = "Percentage", 
                  fill = "samplegroup",add = "jitter",
                  palette = "npg")
  p2 <- facet(p1, facet.by = "newgroup",nrow = 1,scales = "fixed",strip.position = "bottom",panel.labs.background = list(color = NA)) +
    theme(legend.position="right",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(color = "black"),
          strip.placement = "outside",
          strip.background.x = element_blank(),
          strip.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,color = "black")) + 
    scale_y_continuous(expand = c(0,0)) # y轴从0開始
  ggsave("cell_percentage_boxplot.png",p2,height = opt$height,width = opt$width,dpi = 300)
  ggsave("cell_percentage_boxplot.pdf",p2,height = opt$height,width = opt$width)

  ##barplot
  p3 <- ggbarplot(mydata, x = "samplegroup",y = "Percentage", 
                  fill = "samplegroup",facet.by = "newgroup",
                  add = "mean_se",error.plot = "upper_errorbar",
                  palette = "npg",ggtheme = theme_classic())
  p4 <- facet(p3, facet.by = "newgroup",nrow = 1,scales = "fixed",strip.position = "bottom",panel.labs.background = list(color = NA)) +
    theme(legend.position="right",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(color = "black"),
          strip.placement = "outside",
          strip.background.x = element_blank(),
          strip.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,color = "black")) + 
    scale_y_continuous(expand = c(0,0)) # y轴从0開始

  ggsave("cell_percentage_barplot.png",p4,height = opt$height,width = opt$width,dpi = 300)
  ggsave("cell_percentage_barplot.pdf",p4,height = opt$height,width = opt$width)  
}
