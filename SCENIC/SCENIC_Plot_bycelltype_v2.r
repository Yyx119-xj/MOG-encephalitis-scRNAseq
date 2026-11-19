rm(list = ls())
#setwd("D:/Workspace/R_learning/sc_analysis/SCENIC/")

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
library(optparse)
#library(ggheatmap)
library(reshape2)
library(ggsci)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(limma)
#packageVersion("SCENIC")

option_list <- list(
  make_option(c("-l", "--loom"),action = "store",type = "character",
              help = "The Input loom file"),
  make_option(c("-i", "--rds"),action = "store",type = "character",
              help = "The Input rds file"),
  make_option(c("-r", "--regulons"),action = "store",type = "character",default = "regulons.csv",
              help = "The Input rds file"),
  make_option(c("-c", "--comparelist"),action = "store",type = "character",
              help = "The compare list"),
  make_option(c("-q","--qvalue"),action = "store",type = "numeric",default = "0.05",
              help = "The threshold to filtering qvalue"),
  make_option(c("-f","--fc"),action = "store",type = "numeric",default = "0.5",
              help = "The threshold to filtering foldchange"),
  make_option(c("-o", "--outdir"),action = "store",type = "character",default = "./",
              help = "The outdir,the default is './'")
)

opt <- parse_args(OptionParser(option_list = option_list)) ####解释器

work_dir <- getwd()

# # test
# opt$loom <- "/public/home/wuyf/project/singlecell/RX2025-01094-nT/analysis/analysis_GSE242414/scenic/auc_mtx.loom"
# opt$rds <- "/public/home/wuyf/project/singlecell/RX2025-01094-nT/analysis/analysis_GSE242414/scenic/cells_analysis.rds"
# opt$regulons <- "/public/home/wuyf/project/singlecell/RX2025-01094-nT/analysis/analysis_GSE242414/scenic/regulons.csv"
# opt$comparelist <- "Male:Female"

## 提取 auc_mtx.loom 信息
#inputDir='./outputs/'
#scenicLoomPath=file.path(inputDir,'auc_mtx.loom')
#opt$loom <- "auc_mtx.loom"
loom <- open_loom(opt$loom) 

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
#regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
#tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

##
cat("Cell\t",file = "regulonAUC.txt")
write.table(as.matrix(regulonAUC@assays@data@listData[["AUC"]]),
            "regulonAUC.txt",sep = "\t",quote = F,
            col.names = T,append = T)
#write.csv(as.matrix(regulonAUC@assays@data@listData[["AUC"]]), file = "regulonAUC.csv")

#1.查看富集到的全部的TF
#读入单细胞数据
#opt$rds <- "filter.bycells.rds"
scRNA <- readRDS(opt$rds)
#将每一个细胞的每一个TF的表达量与单细胞的每一个细胞匹配
sub_regulonAUC <- regulonAUC[,match(colnames(scRNA),colnames(regulonAUC))]
dim(sub_regulonAUC)
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(scRNA))
#[1] TRUE  一定要为TRUE

##增加TF插槽并输出rds
regulon_AUC <- regulonAUC@NAMES
scRNA[['TF']] <- CreateAssayObject(data = assay(sub_regulonAUC[regulon_AUC,]))
saveRDS(scRNA,"scRNA_TF.rds")
# DefaultAssay(scRNA) <- 'TF'
# ave_exp <- AverageExpression(scRNA, group.by=c("group"))[['TF']]
# cat('regulons\t',file = paste('regulons_ave_auc_bygroup','.txt',sep=''))
# write.table(ave_exp, file = "regulons_ave_auc_bygroup.txt", quote = F, sep = '\t', col.names = T,row.names = T,append = T)

##
meta <- data.frame(row.names = colnames(scRNA), 
                        celltype = scRNA$customclassif,
                        sample = scRNA$replicate,
                        group = scRNA$group)
#meta$newgroup <- paste(cellTypes$celltype,cellTypes$group,sep = "-")

##########bycelltype
selectedResolution <- "celltype"
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=meta[colnames(sub_regulonAUC), 
                                        selectedResolution]) 
rss=na.omit(rss)

#####
#创建文件夹
rss_dir <- "1_Regulons_RSS"
if(!file.exists(rss_dir)){
  dir.create(rss_dir)
}
setwd(rss_dir)

cat("Cell\t",file = "RSS_bycelltype_all.txt")
write.table(rss,"RSS_bycelltype_all.txt",sep = "\t",quote = F,
             col.names = T,append = T)
##heatmap图
rssPlot <- 
  plotRSS(
    rss,
    zThreshold = 2,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr = 0.1,
    varName = "cellType",
    #col.low = '#330066',
    #col.mid = '#66CC66',
    #col.high = '#FFCC33'
  )
plotly::ggplotly(rssPlot$plot)
ggsave("RSS_bycelltype_filter.png",width = 5,height = 10,dpi = 300)
ggsave("RSS_bycelltype_filter.pdf",width = 5,height = 10)

##pheatmap展示
rss_data <- rssPlot$plot$data
rss_data<-dcast(rss_data, 
                Topic~rss_data$cellType,
                value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]

cat("Cell\t",file = "RSS_bycelltype_filter.txt")
write.table(rss_data,"RSS_bycelltype_filter.txt",sep = "\t",quote = F,
            col.names = T,append = T)

pheatmap::pheatmap(rss_data,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                   scale = "row",filename = "RSS_bycelltype_filter2.pdf",
                   width = 6,height = 10)
pheatmap::pheatmap(rss_data,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                   scale = "row",filename = "RSS_bycelltype_filter2.png",
                   width = 6,height = 10)

####rss特异性TF结果
B_rss <- as.data.frame(rss)
#需要作图的细胞类型
celltype <- unique(meta$celltype) %>% as.character()
rssRanklist <- list()

for(i in 1:length(celltype)) {
  
  data_rank_plot <- cbind(as.data.frame(rownames(B_rss)),
                          as.data.frame(B_rss[,celltype[i]]))#提取数据
  
  colnames(data_rank_plot) <- c("TF", "celltype")
  data_rank_plot=na.omit(data_rank_plot)#去除NA
  data_rank_plot <- data_rank_plot[order(data_rank_plot$celltype,decreasing=T),]#降序排列
  data_rank_plot$rank <- seq(1, nrow(data_rank_plot))#添加排序
  
  p <- ggplot(data_rank_plot, aes(x=rank, y=celltype)) + 
    geom_point(size=3, shape=16, color="#1F77B4",alpha =0.4)+
    geom_point(data = data_rank_plot[1:6,],
               size=3, color='#DC050C')+ #选择前6个标记，自行按照需求选择
    theme_bw()+
    theme(axis.title = element_text(colour = 'black', size = 12),
          axis.text = element_text(colour = 'black', size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x='Regulons Rank', y='Specificity Score',title =celltype[i])+
    geom_text_repel(data= data_rank_plot[1:6,],
                    aes(label=TF), color="black", size=3, fontface="italic", 
                    arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                    point.padding = 0.3, segment.color = 'black', 
                    segment.size = 0.3, force = 1, max.iter = 3e3)
  #rssRanklist[[i]] <- p
  ggsave(filename = paste(celltype[i],"rssRank.png",sep = "_"),p,
         width = 4,height = 7,dpi = 300)
  ggsave(filename = paste(celltype[i],"rssRank.pdf",sep = "_"),p,
         width = 4,height = 7)  
  
}
setwd(work_dir)

#####SCENIC的差异分析
#创建文件夹
de_dir <- "2_Regulons_DEG"
if(!file.exists(de_dir)){
  dir.create(de_dir)
}
setwd(de_dir)

#opt$comparelist <- "ARDS:Healthy,ARDS_Sepsis:Sepsis"
compare <- strsplit(opt$comparelist,",")[[1]]

anaAUC <- regulonAUC
anaAUC <- anaAUC[onlyNonDuplicatedExtended(rownames(anaAUC)),]

for (i in 1:length(compare)){
  com <- strsplit(compare[i],":")[[1]]
  deg_dir <- paste(com[1],"vs",com[2],sep = "_")
  if(!file.exists(deg_dir)){
    dir.create(deg_dir)
  }
  setwd(deg_dir)
  
  for (j in celltype){
    #选取需要分析的细胞
    cellinfo <- subset(meta,celltype==j & group%in%com)
    targets <- data.table(FileName=rownames(cellinfo),Target=cellinfo$group)

    #任何一组的细胞数小于5就会跳过
    if(sum(targets$Target==com[1]) <= 5 || sum(targets$Target==com[2]) <= 5){next}
    
    #接下来其实就是常规的差异分析了
    lev<-unique(targets$Target)
    f <- factor(targets$Target, levels=lev) 
    design <- model.matrix(~0+f) #样本矩阵
    colnames(design) <- lev #design矩阵重命名
    
    eset=getAUC(anaAUC) #获取所有细胞TF的AUC值
    eset=eset[,targets$FileName]#我们只选取要分析的T细胞
    eset<-t(scale(t(eset))) #scale
    
    #比较分析limma
    tmp <- paste(com[1],com[2],sep = "-")
    cont.wt <- makeContrasts(tmp,levels=design)
    fit <- lmFit(eset, design)
    fit2 <- contrasts.fit(fit, cont.wt) 
    fit2 <- eBayes(fit2) 
    tT=topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)#差异分析结果
    
    tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
    colnames(tT)=c("FDR","P.Value","logFC")#选择关键的三列即可
    rownames(tT)<-gsub("\\(","",rownames(tT))
    rownames(tT)<-gsub("\\)","",rownames(tT))
    rownames(tT)<-gsub("\\+","",rownames(tT))
    
    #保存结果
    name <- paste(j,"_",com[1],"_vs_",com[2],".txt",sep = "")
    cat("Symbol\t",file = name)
    write.table(tT,file=name,sep="\t",quote=F,append = T)
    
    #接下来就可以按照不同的阈值筛选显著的TF了，其实这样的差异分析对于结果的解释和
    #应用上更有说服力，这样我们可以直观的看到某些转录因子在不同细胞状态下的差异
    #可视化的话可以用热图、火山图等展示。
    logFoldChange <- opt$fc
    adjustP <- opt$qvalue
    diffSig <- tT[with(tT, (abs(logFC)>logFoldChange & FDR < adjustP )), ]#显著筛选
    diffSig=na.omit(diffSig)
    
    #保存结果
    name <- paste(j,"_",com[1],"_vs_",com[2],"_diffSig.txt",sep = "")
    cat("Symbol\t",file = name)
    write.table(diffSig,file=name,sep="\t",quote=F,append = T)
    
    ##火山图
    #提取top5的差异基因在图上展示
    top1 <- arrange(diffSig,logFC) %>% filter(logFC<0) %>% head(n=5)
    top2 <- arrange(diffSig,desc(logFC)) %>% filter(logFC>0) %>% head(n=5)
    top <- rbind(top1,top2) 
    top$Symbol <- rownames(top)
    if(length(rownames(top))==0){next}
    
    tT$Significant <- ifelse(tT$FDR < 0.05 & abs(tT$logFC) >= 0.5, 
                             ifelse(tT$logFC >= 0.5, "Up", "Down"), "Other")
    ggplot(
      # 数据、映射、颜色
      tT, aes(x = logFC, y = -log10(FDR))) +
      geom_point(aes(color = Significant), size=2) +
      scale_color_manual(values = c("blue","grey", "red")) +
      # 注释
      geom_text_repel(
        data = top,
        aes(label = Symbol),
        size = 3,
        force = 1,
        nudge_x = -0.2
        #nudge_y = 0.2
        # box.padding = unit(0.35, "lines"),
        # point.padding = unit(0.3, "lines")
      ) +
      # 辅助线
      #geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      #geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
      # 坐标轴
      #scale_y_continuous(limits = c(0, 20))+
      labs(x="log2(fold change)",
           y="-log10(FDR)") +
      # 图例
      theme_bw(base_size = 12)+
      theme(legend.position = "bottom")+
      ggtitle(paste0(deg_dir," ",j))
    #保存结果
    name <- paste(j,"_",com[1],"_vs_",com[2],"_top5_volcano.png",sep = "")
    ggsave(filename = name,width = 6,height = 7,dpi = 300)
    name <- paste(j,"_",com[1],"_vs_",com[2],"_top5_volcano.pdf",sep = "")
    ggsave(filename = name,width = 6,height = 7)
    
  }
  setwd("../")
}

setwd(work_dir)

#=======================================================================================
#pySCENIC 转录因子与靶基因
#=======================================================================================
#我们还可以可视化感兴趣的TF调控的靶基因
#TF与靶基因的关系在pyscenic分析得到的第二个文件

sce_regulons <- read.csv(opt$regulons)
sce_regulons <- sce_regulons[-2, ]
colnames(sce_regulons) <- sce_regulons[1,]
sce_regulons <- sce_regulons[-1, ]
colnames(sce_regulons) <- c("TF","ID","AUC","NES","MotifSimilarityQvalue","OrthologousIdentity",
                            "Annotation","Context","TargetGenes","RankAtMax")

#sce_regulons <- sce_regulons[which(sce_regulons$AUC>0.1),]
sce_regulons <- select(sce_regulons,TF,TargetGenes,ID,AUC,NES)

sce_regulons$TargetGenes <-gsub("\\[","",sce_regulons$TargetGenes)
sce_regulons$TargetGenes <-gsub("\\]","",sce_regulons$TargetGenes)
sce_regulons$TargetGenes <-gsub("\\(","",sce_regulons$TargetGenes)
sce_regulons$TargetGenes <-gsub("\\)","",sce_regulons$TargetGenes)
sce_regulons$TargetGenes <-gsub("\\'","",sce_regulons$TargetGenes)

sce_regulons_top5 <- sce_regulons %>% group_by(TF) %>% 
  slice_max(AUC,n = 5)
sce_regulons_top10 <- sce_regulons %>% group_by(TF) %>% 
  slice_max(AUC,n = 10)

#保存文件
reg_dir <- "3_Regulons_plot"
if(!file.exists(reg_dir)){
  dir.create(reg_dir)
}
setwd(reg_dir)

write.table(sce_regulons,file = "regulons.txt",sep = "\t",
            quote = F,row.names = F,col.names = T)
#write.table(sce_regulons_top5,file = "regulons_AUC_top5.txt",sep = "\t",
#            quote = F,row.names = F,col.names = T)
#write.table(sce_regulons_top10,file = "regulons_AUC_top10.txt",sep = "\t",
#            quote = F,row.names = F,col.names = T)
setwd(work_dir)
#后续使用perl在linux上处理更方便，然后使用cytoscape画图


