# setwd("D:/workspace/R/sc_analysis/correlation/")
rm(list = ls())

####加载包
library(Seurat)
library(tidyverse)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(reshape2)
library(optparse)

option_list <- list(
  make_option(c("-r", "--rdsfile"),action = "store",type = "character",
              help = "The Input rdsfile"),
  make_option(c("--ident1"),action = "store",type = "character",default = "seurat_clusters",
              help = "The first ident"),
  make_option(c("--ident2"),action = "store",type = "character",default = "",
              help = "The second ident"),
  make_option(c("-m", "--method"),action = "store",type = "character",default = "spearman",
              help = "pearson or spearman, default is spearman"),
  make_option(c("-o", "--outdir"),action = "store",type = "character",default = "./",
              help = "The outdir,the default is './'")
)

opt <- parse_args(OptionParser(option_list = option_list)) ####解释器

#创建文件夹
if(!file.exists(opt$outdir)){
  dir.create(opt$outdir)
}
setwd(opt$outdir) 
workdir <- getwd()

#测试
# opt$rdsfile <- "D:/workspace/R/sc_analysis/CD4+T.rds"
# opt$ident1 <- "seurat_clusters"
# opt$ident2 <- "group"
# opt$method <- "spearman"

#读数据
scRNA <- readRDS(opt$rdsfile)

#colnames(scRNA@meta.data)

if(opt$ident2 == ""){
  if(opt$ident1 == "seurat_clusters"){
    Idents(scRNA) <- paste0("C",scRNA$seurat_clusters)
  }else{
    Idents(scRNA) <- scRNA@meta.data[,opt$ident1]
  }
}else{
  if(opt$ident1 == "seurat_clusters"){
    Idents(scRNA) <- paste0("C",scRNA$seurat_clusters,"_",scRNA@meta.data[,opt$ident2])
  }else{
    Idents(scRNA) <- paste0(scRNA@meta.data[,opt$ident1],"_",scRNA@meta.data[,opt$ident2])
  }
}

exp <- AverageExpression(scRNA)$RNA
features <- names(tail(sort(apply(exp, 1, sd)),2000))#提取2000个高变基因，也就是标准差最大的基因
exp_sub <- exp[which(row.names(exp)%in% features),]
exp_cor <- cor(exp_sub, method = opt$method)
cat('cell\t',file = paste(opt$method,'_cor','.txt',sep=''))
write.table(exp_cor, file = paste(opt$method,'_cor','.txt',sep=''), 
            quote = F, sep = '\t', col.names = T,row.names = T,append = T)

exp_cor_melt <- melt(exp_cor,measure.vars=colnames(exp_cor))
colnames(exp_cor_melt) <- c("cell1","cell2","cor")
write.table(exp_cor_melt, file = paste0(opt$method,"_cor_melt.txt"),
            quote = F, sep = '\t', col.names = T,row.names = F)

#做热图
pheatmap::pheatmap(exp_cor, 
                   color = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")),
                   filename = "spearman_cor.pdf",
                   width = 8,
                   height = 7
)
pheatmap::pheatmap(exp_cor, 
                   color = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")),
                   filename = "spearman_cor.png",
                   width = 8,
                   height = 7
)




