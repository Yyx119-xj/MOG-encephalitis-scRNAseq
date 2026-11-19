# setwd("D:/workspace/R/sc_analysis/scRNAtoolVis/")

rm(list = ls())

suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option(c("-r", "--rdsfile"),action = "store",type = "character",default = "",
              help = "The Input rdsfile"),
  make_option(c("-m", "--markerstab"),action = "store",type = "character",default = "",
              help = "The Input markers tab"),
  make_option(c("--width"),action = "store",type = "numeric",default = "10",
              help = "The width"),
  make_option(c("--height"),action = "store",type = "numeric",default = "6",
              help = "The height")
)

opt <- parse_args(OptionParser(option_list = option_list)) ####解释器

#测试
# opt$rdsfile <- "../cells_analysis.rds"
# opt$markerstab <- "cellmarkers_export.xls"

#创建文件夹
workdir <- getwd()

#加载包
library(scRNAtoolVis)
library(Seurat)
library(patchwork)
library(scales)
library(ggsci)
library(gtools)

#读取数据
sc <- readRDS(opt$rdsfile)
DefaultAssay(sc) <- "RNA"
sc$cellcluster <- factor(sc$cellcluster,levels=mixedsort(unique(sc$cellcluster)))
sc$customclassif <- factor(sc$customclassif,levels=mixedsort(unique(sc$customclassif)))

markers <- read.table(opt$markerstab,header = T,sep = "\t",
                    check.names = F,com = "",quote = "")
markers <- markers[mixedorder(markers$cluster),]
markers <- markers[!duplicated(markers$gene),]

#dotplot
#by customclassif
Idents(sc) <- "customclassif"
p1 <- jjDotPlot(object = sc,
                id = "customclassif",
                markerGene = markers,
                plot.margin = c(4,1,1,1),
                xtree = F,
                ytree = F,
                anno = T,
                aesGroName = "cluster",
                yPosition = length(unique(Idents(sc))) + 1,
                x.text.angle = 45,
                x.text.hjust = 1,
                x.text.vjust = 1,
                gene.order = markers$gene
                )
ggsave("cellmarkers_dotplot_bycelltype.png",p1,width = opt$width,height = opt$height,dpi = 300)
ggsave("cellmarkers_dotplot_bycelltype.pdf",p1,width = opt$width,height = opt$height)

#by cellcluster
Idents(sc) <- "cellcluster"
p2 <- jjDotPlot(object = sc,
                id = "cellcluster",
                markerGene = markers,
                plot.margin = c(4,1,1,1),
                xtree = F,
                ytree = F,
                anno = T,
                aesGroName = "cluster",
                yPosition = length(unique(Idents(sc))) + 1,
                x.text.angle = 45,
                x.text.hjust = 1,
                x.text.vjust = 1,
                gene.order = markers$gene
                )
ggsave("cellmarkers_dotplot_bycluster.png",p2,width = opt$width,height = opt$height,dpi = 300)
ggsave("cellmarkers_dotplot_bycluster.pdf",p2,width = opt$width,height = opt$height)

#heatmap
#by customclassif
tryCatch(
    {
      Idents(sc) <- "customclassif"
      p3 <- AverageHeatmap(object = sc,
                            # group.by = "customclassif",
                            markerGene = markers$gene)
      png("cellmarkers_average_heatmap_bycelltype.png",width = opt$height * 30,height = opt$width * 20,units = "mm",res = 300)
      print(p3)
      dev.off()
      pdf("cellmarkers_average_heatmap_bycelltype.pdf",width = opt$height,height = opt$width*2/3)
      print(p3)
      dev.off()

      #by cellcluster
      Idents(sc) <- "cellcluster"
      p4 <- AverageHeatmap(object = sc,
                            # group.by = "cellcluster",
                            markerGene = markers$gene)
      png("cellmarkers_average_heatmap_bycluster.png",width = opt$height * 30,height = opt$width * 20,units = "mm",res = 300)
      print(p4)
      dev.off()
      pdf("cellmarkers_average_heatmap_bycluster.pdf",width = opt$height,height = opt$width*2/3)
      print(p4)
      dev.off()
    },
    # warning = function(w){message("warning")},
    error = function(e){message("error")}
)

tryCatch(
    {
      Idents(sc) <- "customclassif"
      p3 <- averageHeatmap(object = sc,
                            # group.by = "customclassif",
                            markerGene = markers$gene)
      png("cellmarkers_average_heatmap_bycelltype.png",width = opt$height * 30,height = opt$width * 20,units = "mm",res = 300)
      print(p3)
      dev.off()
      pdf("cellmarkers_average_heatmap_bycelltype.pdf",width = opt$height,height = opt$width*2/3)
      print(p3)
      dev.off()

      #by cellcluster
      Idents(sc) <- "cellcluster"
      p4 <- averageHeatmap(object = sc,
                            # group.by = "cellcluster",
                            markerGene = markers$gene)
      png("cellmarkers_average_heatmap_bycluster.png",width = opt$height * 30,height = opt$width * 20,units = "mm",res = 300)
      print(p4)
      dev.off()
      pdf("cellmarkers_average_heatmap_bycluster.pdf",width = opt$height,height = opt$width*2/3)
      print(p4)
      dev.off()
    },
    # warning = function(w){message("warning")},
    error = function(e){message("error")}
)

