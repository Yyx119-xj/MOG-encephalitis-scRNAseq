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
# https://satijalab.org/seurat/pbmc3k_tutorial.html
# https://github.com/CostaLab/sc.combined_seurat_pipeline
#####################################################################################

#####################################################################################
#####参数获取
#####################################################################################
suppressPackageStartupMessages(library(optparse))
## 参数读取
option_list <- list(
    make_option(
        c("-n","--ncpus"),
        action = "store",
        type = "integer",
        default = 6,
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
        c("--MinTotalUMI"),
        action = "store",
        type = "integer",
        default = 1000,
        help = "minimal total UMI count"
    ),
    make_option(
        c("--MinGenes"),
        action = "store",
        type = "integer",
        default = 500,
        help = "minimal genes"
    ),
    make_option(
        c("--MaxMT"),
        action = "store",
        type = "integer",
        default = 15,
        help = "maximal percentage of mitochondria"
    ),
    make_option(
        c("--MinCellsInGene"),
        action = "store",
        type = "integer",
        default = 5,
        help = "minimal cells with UMI>0 in one gene"
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
        c("-f","--format"),
        action = "store",
        type = "character",
        default = "",
        help = "The count matrix file format, 10x or count or h5"
    ),
    make_option(
        c("-b","--bdfiles"),
        action = "store",
        type = "character",
        default = "",
        help = "The count matrix file of sample"
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
        help = "integrate methods, cca , harmony or false, default is cca"
    ),  
    # make_option(
    #     c("-a","--annofiles"),
    #     action = "store",
    #     type = "character",
    #     help = "The Annotation file of samples"
    # ),
    make_option(
        c("--skipfilter"),
        action = "store",
        type = "character",
        default = "false",
        help = "false or true"
    ),
    make_option(
        c("--rdsfile"),
        action = "store",
        type = "character",
        default = "",
        help = "The data object rds file"
    ),
    make_option(
        c("--keepclusters"),
        action = "store",
        type = "character",
        default = "",
        help = "The keep cluster"
    ),
    make_option(
        c("--skipcellcycle"),
        action = "store",
        type = "character",
        default = "true",
        help = "false or true"
    ),
    make_option(
        c("--sgenelist"),
        action = "store",
        type = "character",
        default = "",
        help = "The keep cluster"
    ),
    make_option(
        c("--g2mgenelist"),
        action = "store",
        type = "character",
        default = "",
        help = "The keep cluster"
    ),
    make_option(
        c("--rmbgRNA"),
        action = "store",
        type = "character",
        default = "false",
        help = "remove background RNA, default false"
    )                     
)
opt <- parse_args(OptionParser(option_list = option_list))

#加载包
options(stringsAsFactors = FALSE)
`%ni%` <- Negate(`%in%`)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(scater)
library(grDevices)
library(gplots)
library(Seurat)
library(patchwork)
library(hdf5r)
library(DoubletFinder)
library(stats)
library(ggthemes)
library(reshape2)
library(harmony)
library(future.apply)
library(doParallel)
library(SoupX)

#
registerDoParallel(cores=opt$ncpus)
plan("multicore", workers = opt$ncpus)
options(future.globals.maxSize = opt$MaxMemMega * 1024^2)

#
min_total_UMI = opt$MinTotalUMI
min_genes = opt$MinGenes
max_MT = opt$MaxMT
min_cells_in_gene = opt$MinCellsInGene

# 测试参数
# opt$samplelist <- "ITP_1,HC_1"
# opt$listname <- "ITP,HC"

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

#####################################################################################
## 数据读取，数据对象建立
#####################################################################################

countfiles <- strsplit(opt$bdfiles,",")[[1]]
# annofiles <- strsplit(opt$annofiles,",")[[1]]
samplenames <- strsplit(opt$samplelist,",")[[1]]
listnames <- strsplit(opt$listname,",")[[1]]
umiformat <- strsplit(opt$format,",")[[1]]
groupnum = length(unique(listnames))

sc.list <- list()

cell_filter <- function(i){
    testcountfile <- ""
    if(opt$rmbgRNA == "false"){
        if(umiformat[i] == "10x"){
            #10x标准输出，三个文件
            testcountfile <- Read10X(data.dir = countfiles[i])
        }
        if(umiformat[i] == "BGI"){
            #BGI标准输出，三个文件
            testcountfile <- Read10X(data.dir = countfiles[i],gene.column = 1)
        }         
        if(umiformat[i] == "count"){
            #row为gene, col为cell
            testcountfile <- read.table(file = countfiles[i],header = T,com = '',sep = "\t",quote = NULL,check.names = F,row.names = 1)
            # testcountfile <- t(testcountfile)
        }
        if(umiformat[i] == "h5"){
            testcountfile= Read10X_h5(countfiles[i])
        }  
        if(umiformat[i] == "rds"){
            #这里的输入rds必须是名字已经修改过的
            tmp = readRDS(countfiles[i])
            testcountfile <- GetAssayData(tmp,assay = "RNA",slot = "count")
        }        
    }

    if(opt$rmbgRNA == "true"){
        if(umiformat[i] == "10x"){
            raw_dir <- file.path(countfiles[i], "raw_feature_bc_matrix") 
            filt_dir <- file.path(countfiles[i], "filtered_feature_bc_matrix") 
            raw_matrix <- Read10X(data.dir = raw_dir)
            filt_matrix <- Read10X(data.dir = filt_dir)
            # 找到 raw 和 filtered 矩阵中共同存在的基因
            common_genes <- intersect(rownames(raw_matrix), rownames(filt_matrix))
            # 将两个矩阵按照相同的交集基因、相同的顺序进行过滤
            raw_matrix <- raw_matrix[common_genes, ]
            filt_matrix <- filt_matrix[common_genes, ]
            # 3. 创建 SoupChannel 对象 (此时绝对不会再报错了)
            sc <- SoupChannel(raw_matrix, filt_matrix)
            # 创建临时 Seurat 对象来获取聚类
            seu_tmp <- CreateSeuratObject(counts = filt_matrix)
            seu_tmp <- NormalizeData(seu_tmp, verbose = FALSE)
            seu_tmp <- FindVariableFeatures(seu_tmp, verbose = FALSE)
            seu_tmp <- ScaleData(seu_tmp, verbose = FALSE)
            seu_tmp <- RunPCA(seu_tmp, verbose = FALSE)
            seu_tmp <- RunUMAP(seu_tmp, dims = 1:20, verbose = FALSE)
            seu_tmp <- FindNeighbors(seu_tmp, dims = 1:20, verbose = FALSE)
            seu_tmp <- FindClusters(seu_tmp, resolution = 0.5, verbose = FALSE)
            # 将聚类信息传递给 SoupX
            meta <- seu_tmp@meta.data
            umap <- seu_tmp@reductions$umap@cell.embeddings
            sc <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
            sc <- setDR(sc, umap)
            # 估算污染并校正
            sc <- autoEstCont(sc,forceAccept = T)
            testcountfile <- adjustCounts(sc, roundToInt = TRUE)
        }
        if(umiformat[i] == "BGI"){
            raw_dir <- file.path(countfiles[i], "raw_matrix") 
            filt_dir <- file.path(countfiles[i], "filter_matrix") 
            raw_matrix <- Read10X(data.dir = raw_dir,gene.column = 1)
            filt_matrix <- Read10X(data.dir = filt_dir,gene.column = 1)
            # 找到 raw 和 filtered 矩阵中共同存在的基因
            common_genes <- intersect(rownames(raw_matrix), rownames(filt_matrix))
            # 将两个矩阵按照相同的交集基因、相同的顺序进行过滤
            raw_matrix <- raw_matrix[common_genes, ]
            filt_matrix <- filt_matrix[common_genes, ]
            # 3. 创建 SoupChannel 对象 (此时绝对不会再报错了)
            sc <- SoupChannel(raw_matrix, filt_matrix)
            # 创建临时 Seurat 对象来获取聚类
            seu_tmp <- CreateSeuratObject(counts = filt_matrix)
            seu_tmp <- NormalizeData(seu_tmp, verbose = FALSE)
            seu_tmp <- FindVariableFeatures(seu_tmp, verbose = FALSE)
            seu_tmp <- ScaleData(seu_tmp, verbose = FALSE)
            seu_tmp <- RunPCA(seu_tmp, verbose = FALSE)
            seu_tmp <- RunUMAP(seu_tmp, dims = 1:20, verbose = FALSE)
            seu_tmp <- FindNeighbors(seu_tmp, dims = 1:20, verbose = FALSE)
            seu_tmp <- FindClusters(seu_tmp, resolution = 0.5, verbose = FALSE)
            # 将聚类信息传递给 SoupX
            meta <- seu_tmp@meta.data
            umap <- seu_tmp@reductions$umap@cell.embeddings
            sc <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
            sc <- setDR(sc, umap)
            # 估算污染并校正
            sc <- autoEstCont(sc,forceAccept = T)
            testcountfile <- adjustCounts(sc, roundToInt = TRUE)
        }        
    }

    colnames(testcountfile) = paste(samplenames[i],colnames(testcountfile),sep='_')
    test <- CreateSeuratObject(counts = testcountfile, project = "test",min.cells = 1, min.features = 1) 

    ## 过滤检出基因细胞数过少的基因
    # 提取计数
    counts <- GetAssayData(object = test, slot = "counts")
    # 根据在每个细胞的计数是否大于0为每个基因输出一个逻辑向量
    nonzero <- counts > 0
    # 将所有TRUE值相加，如果每个基因的TRUE值超过10个，则返回TRUE。
    keep_genes <- Matrix::rowSums(nonzero) >= min_cells_in_gene
    test_discard_genes <- Matrix::rowSums(nonzero) < min_cells_in_gene
    # 仅保留那些在10个以上细胞中表达的基因
    filtered_counts <- counts[keep_genes, ]
    # 重新赋值给经过过滤的Seurat对象
    # test <- CreateSeuratObject(filtered_counts, meta.data = test@meta.data)
    test <- CreateSeuratObject(counts = filtered_counts, data = filtered_counts, meta.data = test@meta.data)
    test[["group"]] <- listnames[i]
    test[["replicate"]] <- samplenames[i]
        
    #####################################################################################
    ## cell筛选
    #####################################################################################

    # 线粒体基因，人小鼠和大鼠，其它物种需要自行修改
    # 大鼠NCBI refseq基因中的线粒体基因是|^ATP6$|^ATP8$|^COX1$|^COX2$|^COX3$|^CYTB$|^ND1$|^ND2$|^ND3$|^ND4$|^ND4L$|^ND5$|^ND6$
    test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^MT-|^mt-|^Mt-|^ATP6$|^ATP8$|^COX1$|^COX2$|^COX3$|^CYTB$|^ND1$|^ND2$|^ND3$|^ND4$|^ND4L$|^ND5$|^ND6$")    

    # 核糖体基因
    test[["percent.ribo"]] <- PercentageFeatureSet(test, "^RP[SL]|^Rp[sl]")

    # 红细胞
    # Percentage hemoglobin genes - includes all genes starting with HB except HBP.
    # And finally, with the same method we will calculate proportion hemoglobin genes, which can give an indication of red blood cell contamination.
    test[["percent.hb"]] <- PercentageFeatureSet(test, "^HB[^(P)]|^Hb[^(p)]")
    # test[["percent.plat"]] <- PercentageFeatureSet(test, "^PECAM1$|^PF4$|^Pecam1$|^Pf4$")

    #################################################
    ## Visualize QC metrics as a violin plot, test
    ################################################
    # QC
    p1 <- VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 5, pt.size = 0.1)
    # p1 <- VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.1)
    pdf(file=paste("QC_all_",samplenames[i],".pdf",sep=""),width = 20,height = 4)
    print(p1)
    dev.off()

    png(file = paste("QC_all_",samplenames[i],".png",sep=""),width = 600,height = 120,units = "mm",res = 300,pointsize = 2)
    par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
    print(p1)
    dev.off()

    # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
    # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
    p1 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.1)
    p2 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1)
    pdf(file=paste("QC_FeatureScatter_",samplenames[i],".pdf",sep=""),width = 10,height = 4)
    print(p1 + p2)
    dev.off()

    png(file = paste("QC_FeatureScatter_",samplenames[i],".png",sep=""),width = 250,height = 100,units = "mm",res = 300,pointsize = 2)
    par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
    print(p1 + p2)
    dev.off()

    qc.total <- test$nCount_RNA > 0
    # qc.lib <- test$sum < 4e5
    qc.lib <- (test$nCount_RNA < min_total_UMI & test$nCount_RNA > 0)
    # qc.nexprs <- test$detected < 5e3
    qc.nexprs <- (test$nFeature_RNA < min_genes & test$nCount_RNA > 0)
    # qc.spike <- test$altexps_ERCC_percent > 10
    qc.mito <- (test$percent.mt > max_MT & test$nCount_RNA > 0)

    # test <- subset(test, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    test <- subset(test, subset = nFeature_RNA >= min_genes & nCount_RNA >= min_total_UMI & percent.mt <= max_MT)

    # 标准化
    test <- NormalizeData(test) %>% FindVariableFeatures(nfeatures = 2000)
    if(opt$nmethod == "log"){
        test <- test %>% ScaleData()
    }else if(opt$nmethod == "sct"){
        test <- SCTransform(test, variable.features.n = 2000)
    }
    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(test), 10)

    # plot variable features with and without labels
    pdf(file=paste("QC_highly_variable_features_",samplenames[i],".pdf",sep=""),width = 6,height = 4)
    plot1 <- VariableFeaturePlot(test)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    # plot1 + plot2
    print(plot2)
    dev.off()

    png(file = paste("QC_highly_variable_features_",samplenames[i],".png",sep=""),width = 150,height = 100,units = "mm",res = 300,pointsize = 2)
    par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
    plot1 <- VariableFeaturePlot(test)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    # plot1 + plot2
    print(plot2)
    dev.off()

    # 双细胞过滤 
    # remove doublet using DoubletFinder ; 
    # https://nbisweden.github.io/workshop-sc.combinedseq/labs/compiled/seurat/seurat_01_qc.html  
    test <-  RunPCA(test, verbose = F, npcs = 20)
    pc.num <- 1:15
    test <- RunUMAP(test, dims = pc.num, verbose = F)
    test <- FindNeighbors(test, dims = pc.num) %>% FindClusters(resolution = 0.3)
    ## 寻找最优pK值
    sweep.res.list <- paramSweep_v3(test, PCs = pc.num, sct = T)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    # pK_bcmvn <- 0.2
    print("pK_bcmvn:")
    print(pK_bcmvn)

    ## 排除不能检出的同源doublets，优化期望的doublets数量
    # qc.total <- test$nCount_RNA > 0
    DoubletRate = ncol(test)*8*1e-6                     # 5000细胞对应的doublets rate是3.9%
    homotypic.prop <- modelHomotypic(test$seurat_clusters)   # 最好提供celltype
    nExp_poi <- round(DoubletRate*ncol(test)) 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    ## 使用确定好的参数鉴定doublets
    test <- doubletFinder_v3(test, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = T)                            
    ## 结果展示，分类结果在pbmc@meta.data中
    DF.name = colnames(test@meta.data)[grepl("DF.classification", colnames(test@meta.data))]
    p1 <- DimPlot(test, reduction = "umap", group.by = DF.name)
    pdf(file=paste("DimPlot_doublet_",samplenames[i],".pdf",sep=""),width = 5,height = 4)
    print(p1)
    dev.off()

    png(file = paste("DimPlot_doublet_",samplenames[i],".png",sep=""),width = 125,height = 100,units = "mm",res = 300,pointsize = 2)
    par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
    print(p1)
    dev.off()

    p1 <- VlnPlot(test, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
    pdf(file=paste("VlnPlot_doublet_",samplenames[i],".pdf",sep=""),width = 5,height = 4)
    print(p1)
    dev.off()

    png(file = paste("VlnPlot_doublet_",samplenames[i],".png",sep=""),width = 125,height = 100,units = "mm",res = 300,pointsize = 2)
    par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
    print(p1)
    dev.off()

    test.doublet <- test[, test@meta.data[, DF.name] == "Doublet"]
    qc.doublet <- test.doublet$nCount_RNA > 0
    test <- test[, test@meta.data[, DF.name] == "Singlet"]
    qc.keep <- test$nCount_RNA > 0

    # discard <- qc.lib | qc.nexprs | qc.spike | qc.mito
    # discard <- qc.lib | qc.nexprs | qc.mito | qc.doublet
    # discard <- qc.lib | qc.nexprs | qc.mito | test@meta.data[, DF.name] == "Doublet"
    # discard

    #增加过滤后每个样本的ncount中位数和nfeature中位数
    median_count <- test@meta.data$nCount_RNA[test$replicate == samplenames[i]] %>% median()
    median_feature <- test@meta.data$nFeature_RNA[test$replicate == samplenames[i]] %>% median()

    # DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs),
    #           SpikeProp=sum(qc.spike), MitoProp=sum(qc.mito), Total=sum(discard))
    dis <- DataFrame(Raw=sum(qc.total), LibSize=sum(qc.lib), NExprs=sum(qc.nexprs),
            MitoProp=sum(qc.mito),DoubletProp=sum(qc.doublet), Total=sum(qc.total)-sum(qc.keep), Keep=sum(qc.keep), ExpressGenes=sum(test_discard_genes), KeepGenes=sum(keep_genes), median_count=median_count, median_feature=median_feature)
    dis
    cat(paste('原始细胞数\t按比对总量过滤','(<',min_total_UMI,')','\t按检出基因数过滤','(<',min_genes,')','\t按线粒体比例过滤','(>',max_MT,'%)','\t按doublet过滤','\t总过滤细胞数\t最终保留的细胞数\t筛选掉少数细胞表达的基因','(<',min_cells_in_gene,')','\t最终保留的基因数\t过滤后count中位数\t过滤后基因中位数\n',sep=''),file=paste(samplenames[i],'_discard_stat','.xls',sep=''))
    write.table(dis,file=paste(samplenames[i],'_discard_stat','.xls',sep=''),
                append=T,quote=F,sep='\t',row.names = FALSE,col.names = FALSE)
    
    # test = test[, test@meta.data[, DF.name] == "Singlet"]

    # Compute the relative expression of each gene per cell
    coum <- test@assays$RNA@counts
    coum <- Matrix::t(Matrix::t(coum)/Matrix::colSums(coum)) * 100
    most_expressed <- order(apply(coum, 1, median), decreasing = T)[30:1]
    print(most_expressed)
    gg <- as.data.frame(t(coum[most_expressed, ]))
    dfs <- melt(gg)
    print(head(dfs))
    ggplot(dfs,aes(x=variable,y=value))+geom_boxplot()+coord_flip()+labs(title="Top30 most expressed genes",y="% total count per cell", x = "")+theme_bw() + theme_paper
    ggsave(
        paste("percentage_of_counts_per_gene_",samplenames[i],".pdf",sep=""),
        width = 200,
        height = 150,
        units = "mm"
    )
    ggsave(
        paste("percentage_of_counts_per_gene_",samplenames[i],".png",sep=""),
        width = 200,
        height = 150,
        units = "mm"
    )

    return(test)
}

#执行cell_filter
if(opt$skipfilter == "false"){
    sc.list <- future_lapply(1:length(countfiles), cell_filter)
    # saveRDS(sc.list,"sc.list.rds")    
}
if(opt$skipfilter == "true" & opt$rdsfile != ""){    
    if(opt$keepclusters == ""){
        sc <- readRDS(opt$rdsfile)
        sc.list <- SplitObject(sc, split.by = "replicate")
    }else{
        sc <- readRDS(opt$rdsfile)
        Idents(sc) <- sc$seurat_clusters
        keeps <- strsplit(opt$keepclusters,",")[[1]]
        sc_sub <- subset(sc, idents= keeps)
        sc.list <- SplitObject(sc_sub, split.by = "replicate")
        # saveRDS(sc.list,"sc.list.rds") 
    }
}

##细胞周期函数
cell_cycle_score <- function(obj = sc.combined){
    DefaultAssay(obj) <- "RNA"
    #读取g2m和s期基因列表
    #默认提供人和小鼠的，如果是其它物种需要传入相应的文件
    if(opt$g2mgenelist != ""){
        g2m_genes <- read.table(opt$g2mgenelist,header = F,check.names = F,com = "",quote = "")
        g2m_genes <- g2m_genes[,1]
    }else{
        g2m_genes <- c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA","Gtse1","Ccnb2","Ckap2","Tpx2","Ube2c","Ndc80","Aurka","Ttk","Tubb4b","Top2a","Gas2l3","Kif11","Hmgb2","G2e3","Aurkb","Kif23","Rangap1","Nusap1","Dlgap5","Nuf2","Nek2","Tmpo","Hjurp","Cdc25c","Birc5","Tacc3","Cdk1","Cks2","Cenpe","Kif2c","Ect2","Cbx5","Ncapd2","Hmmr","Cdca3","Mki67","Ctcf","Cdca2","Bub1","Kif20b","Cdc20","Anln","Lbr","Cks1b","Ckap5","Psrc1","Ckap2l","Smc4","Anp32e","Cenpa","Cenpf","Cdca8")
    }
    if(opt$sgenelist != ""){
        s_genes <- read.table(opt$sgenelist,header = F,check.names = F,com = "",quote = "")
        s_genes <- s_genes[,1]
    }else{
        s_genes <- c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8","Chaf1b","Rad51","Cdc45","Mcm5","Pcna","Gins2","Dscc1","Cdc6","Ubr7","E2f8","Tyms","Wdr76","Clspn","Usp1","Casp8ap2","Pola1","Slbp","Fen1","Blm","Dtl","Ccne2","Gmnn","Rfc2","Tipin","Hells","Rrm1","Mcm4","Rpa2","Rrm2","Brip1","Nasp","Mcm2","Uhrf1","Exo1","Rad51ap1","Ung","Prim1","Cdca7","Mcm6","Msh2")
    }
    g2m_genes <- g2m_genes[g2m_genes %in% rownames(obj)] #提取矩阵中的G2M期marker基因
    s_genes <- s_genes[s_genes %in% rownames(obj)] #提取矩阵中的S期marker基因
    #细胞周期评分用的是Normalize Data（经log转换后的矩阵)
    obj <- CellCycleScoring(obj, g2m.features=g2m_genes, s.features=s_genes) 

    # # Running a PCA on cell cycle genes reveals
    # obj <- RunPCA(obj, features = c(s_genes, g2m_genes)) 
    # p1 <- DimPlot(obj,group.by = 'Phase')
    # ggsave("Phase_PCA_plot_before_regression.png",width = 5,height = 4,dpi = 300)
    # ggsave("Phase_PCA_plot_before_regression.pdf",width = 5,height = 4)

    # if(opt$skipcellcycle == "false"){
    #     if(opt$nmethod == "log"){
    #         #对周期评分进行回归
    #         obj <- ScaleData(obj, vars.to.regress = c("S.Score", "G2M.Score"))
    #         # Running a PCA on cell cycle genes reveals
    #         obj <- RunPCA(obj, features = c(s_genes, g2m_genes)) 
    #         p1 <- DimPlot(obj,group.by = 'Phase')
    #         ggsave("Phase_PCA_plot_after_regression.png",width = 5,height = 4,dpi = 300)
    #         ggsave("Phase_PCA_plot_after_regression.pdf",width = 5,height = 4)
    #     }else if(opt$nmethod == "sct"){
    #         #对周期评分进行回归
    #         obj <- SCTransform(obj, variable.features.n = 2000 , vars.to.regress = c("S.Score", "G2M.Score"), verbose = T)
    #         # Running a PCA on cell cycle genes reveals
    #         obj <- RunPCA(obj, features = c(s_genes, g2m_genes)) 
    #         p1 <- DimPlot(obj,group.by = 'Phase')
    #         ggsave("Phase_PCA_plot_after_regression.png",width = 5,height = 4,dpi = 300)
    #         ggsave("Phase_PCA_plot_after_regression.pdf",width = 5,height = 4)
    #     }
    # }
    return(obj)    
}

##整合数据
if(length(sc.list)==1){
    sc.combined <- sc.list[[1]]
    # #细胞周期
    # tryCatch(
    #     {
    #     sc.combined <- cell_cycle_score(sc.combined)
    #     },
    #     # warning = function(w){message("warning")},
    #     error = function(e){message("error")}
    # )
    DefaultAssay(sc.combined) <- "RNA"   
    sc.combined <- DietSeurat(sc.combined,assays = "RNA") #精简rds，只保留RNA槽
    sc.combined <- cell_cycle_score(sc.combined) #细胞周期
    if(opt$nmethod=="log"){
        if(opt$skipcellcycle == "true"){
            sc.combined <- NormalizeData(sc.combined) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData(vars.to.regress = c("nCount_RNA","percent.mt")) #RNA槽标准化
        }else if(opt$skipcellcycle == "false"){
            sc.combined <- NormalizeData(sc.combined) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData(vars.to.regress = c("nCount_RNA","percent.mt","S.Score","G2M.Score")) #RNA槽标准化        
        }  
    }else if(opt$nmethod=="sct"){
        if(opt$skipcellcycle == "true"){
            sc.combined <- SCTransform(sc.combined, variable.features.n = 2000 , vars.to.regress = c("nCount_RNA","percent.mt"), verbose = T) #SCT标准化，仅仅是为了降维聚类使用  
        }else if(opt$skipcellcycle == "false"){
            sc.combined <- SCTransform(sc.combined, variable.features.n = 2000 , vars.to.regress = c("nCount_RNA","percent.mt","S.Score","G2M.Score"), verbose = T) #SCT标准化，仅仅是为了降维聚类使用        
        }          
    }           
    sc.combined <- RunPCA(sc.combined, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)    
}else{
    if(opt$nmethod=="log" & opt$imethod=="cca"){
        features <- SelectIntegrationFeatures(object.list = sc.list)
        sc.anchors <- FindIntegrationAnchors(object.list = sc.list, anchor.features = features)
        sc.combined <- IntegrateData(anchorset = sc.anchors)
        sc.combined <- cell_cycle_score(sc.combined) #细胞周期
        DefaultAssay(sc.combined) <- "integrated"
        if(opt$skipcellcycle == "true"){
            sc.combined <- sc.combined %>% ScaleData(vars.to.regress = c("nCount_RNA","percent.mt")) #RNA槽标准化
        }else if(opt$skipcellcycle == "false"){
            sc.combined <- sc.combined %>% ScaleData(vars.to.regress = c("nCount_RNA","percent.mt","S.Score","G2M.Score")) #RNA槽标准化        
        }         
        sc.combined <- RunPCA(sc.combined, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)        
    }
    if(opt$nmethod=="log" & opt$imethod=="harmony"){
        sc.combined <- merge(sc.list[[1]], sc.list[2:length(sc.list)]) #合并数据
        DefaultAssay(sc.combined) <- "RNA"
        sc.combined <- DietSeurat(sc.combined,assays = "RNA") #精简rds，只保留RNA槽
        sc.combined <- cell_cycle_score(sc.combined) #细胞周期
        if(opt$skipcellcycle == "true"){
            sc.combined <- NormalizeData(sc.combined) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData(vars.to.regress = c("nCount_RNA","percent.mt")) #RNA槽标准化
        }else if(opt$skipcellcycle == "false"){
            sc.combined <- NormalizeData(sc.combined) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData(vars.to.regress = c("nCount_RNA","percent.mt","S.Score","G2M.Score")) #RNA槽标准化        
        }                        
        sc.combined <- RunPCA(sc.combined, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
        sc.combined <- RunHarmony(sc.combined, group.by.vars = "replicate", max_iter = 20)
    }
    if(opt$nmethod=="sct" & opt$imethod=="harmony"){        
        sc.combined <- merge(sc.list[[1]], sc.list[2:length(sc.list)]) #合并数据
        DefaultAssay(sc.combined) <- "RNA"
        sc.combined <- DietSeurat(sc.combined,assays = "RNA") #精简rds，只保留RNA槽
        sc.combined <- cell_cycle_score(sc.combined) #细胞周期
        if(opt$skipcellcycle == "true"){
            sc.combined <- SCTransform(sc.combined, variable.features.n = 2000 , vars.to.regress = c("nCount_RNA","percent.mt"), verbose = T) #SCT标准化，仅仅是为了降维聚类使用  
        }else if(opt$skipcellcycle == "false"){
            sc.combined <- SCTransform(sc.combined, variable.features.n = 2000 , vars.to.regress = c("nCount_RNA","percent.mt","S.Score","G2M.Score"), verbose = T) #SCT标准化，仅仅是为了降维聚类使用        
        }                
        sc.combined <- RunPCA(sc.combined, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
        sc.combined <- RunHarmony(sc.combined, group.by.vars="replicate", max_iter = 20)
    }
    if(opt$nmethod=="log" & opt$imethod=="false"){        
        sc.combined <- merge(sc.list[[1]], sc.list[2:length(sc.list)]) #合并数据
        DefaultAssay(sc.combined) <- "RNA"
        sc.combined <- DietSeurat(sc.combined,assays = "RNA") #精简rds，只保留RNA槽  
        sc.combined <- cell_cycle_score(sc.combined) #细胞周期
        if(opt$skipcellcycle == "true"){
            sc.combined <- NormalizeData(sc.combined) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData(vars.to.regress = c("nCount_RNA","percent.mt")) #RNA槽标准化
        }else if(opt$skipcellcycle == "false"){
            sc.combined <- NormalizeData(sc.combined) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData(vars.to.regress = c("nCount_RNA","percent.mt","S.Score","G2M.Score")) #RNA槽标准化        
        }        
        sc.combined <- RunPCA(sc.combined, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
    }        
    if(opt$nmethod=="sct" & opt$imethod=="false"){        
        sc.combined <- merge(sc.list[[1]], sc.list[2:length(sc.list)]) #合并数据
        DefaultAssay(sc.combined) <- "RNA"
        sc.combined <- DietSeurat(sc.combined,assays = "RNA") #精简rds，只保留RNA槽
        sc.combined <- cell_cycle_score(sc.combined) #细胞周期
        if(opt$skipcellcycle == "true"){
            sc.combined <- SCTransform(sc.combined, variable.features.n = 2000 , vars.to.regress = c("nCount_RNA","percent.mt"), verbose = T) #SCT标准化，仅仅是为了降维聚类使用  
        }else if(opt$skipcellcycle == "false"){
            sc.combined <- SCTransform(sc.combined, variable.features.n = 2000 , vars.to.regress = c("nCount_RNA","percent.mt","S.Score","G2M.Score"), verbose = T) #SCT标准化，仅仅是为了降维聚类使用        
        }                   
        sc.combined <- RunPCA(sc.combined, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
    }    
}

# saveRDS(sc.combined, file = "./sc_preprocessing.rds")

#作图
p1 <- VizDimLoadings(sc.combined, dims = 1:2, reduction = "pca")
pdf(file=paste("VizPlot_PCA",".pdf",sep=""),width = 8,height = 5)
p1
dev.off()

png(file = paste("VizPlot_PCA",".png",sep=""),width = 200,height = 125,units = "mm",res = 300,pointsize = 1.5)
par(mar = c(5, 5, 2, 2), cex.axis = 1.5, cex.lab = 2)
p1
dev.off()

p1 <- DimPlot(sc.combined, reduction = "pca", group.by = "replicate", split.by = "group",pt.size = 0.01, raster = F )
pdf(file=paste("DimPlot_PCA",".pdf",sep=""),width = 4 * groupnum,height = 4)
p1
dev.off()

png(file = paste("DimPlot_PCA",".png",sep=""),width = 100 * groupnum,height = 100,units = "mm",res = 300,pointsize = 2)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
p1
dev.off()

#
pdf(file=paste("DimHeatmap_PCA",".pdf",sep=""),width = 15,height = 10)
DimHeatmap(sc.combined, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()

png(file = paste("DimHeatmap_PCA",".png",sep=""),width = 450,height = 300,units = "mm",res = 300,pointsize = 15)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
DimHeatmap(sc.combined, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()

#
pdf(file=paste("ElbowPlot",".pdf",sep=""),width = 10,height = 4)
p1 <- ElbowPlot(sc.combined, ndims = 25)
p2 <- ElbowPlot(sc.combined, ndims = 100)
p1+p2
dev.off()

png(file = paste("ElbowPlot",".png",sep=""),width = 250,height = 100,units = "mm",res = 300,pointsize = 2)
par(mar = c(5, 5, 2, 2), cex.axis = 2, cex.lab = 2)
p1 <- ElbowPlot(sc.combined, ndims = 25)
p2 <- ElbowPlot(sc.combined, ndims = 100)
p1+p2
dev.off()

#降维聚类函数
dim_reduction_plot <- function(red){
    sc.combined <- FindNeighbors(sc.combined, reduction = red, dims = 1:30)
    sc.combined <- FindClusters(sc.combined, resolution = 0.4)
    sc.combined <- RunUMAP(sc.combined, dims = 1:30, reduction = red)
    return(sc.combined)   
}

#降维聚类
if(opt$imethod == "harmony"){
    tryCatch(
        {
            sc.combined <- dim_reduction_plot("harmony")
        },
        # warning = function(w){message("warning")},
        error = function(e){message("error")}
    ) 
}else{
    tryCatch(
        {
            sc.combined <- dim_reduction_plot("pca")
        },
        # warning = function(w){message("warning")},
        error = function(e){message("error")}
    )    
}

#保存
saveRDS(sc.combined, file = "./sc_preprocessing.rds")

#Dimplot
tryCatch(
    {
        ## umap
        p1 <- DimPlot(sc.combined, reduction = "umap", group.by = "replicate", pt.size = 0.1, raster = F )
        p2 <- DimPlot(sc.combined, reduction = "umap", group.by = "group", pt.size = 0.1, raster = F )
        p3 <- DimPlot(sc.combined, reduction = "umap", group.by = "Phase", pt.size = 0.1, raster = F )    
        p4 <- DimPlot(sc.combined, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, raster = F, label = T )   

        pdf(file=paste("umap_cluster","pdf",sep="."),width = 14,height = 10)
        print((p1|p2)/(p3|p4) + theme_paper)
        dev.off()
        png(file = paste("umap_cluster",".png",sep=""),width = 350,height = 250,units = "mm",res = 300,pointsize = 10)
        print((p1|p2)/(p3|p4))
        dev.off() 
    },
    # warning = function(w){message("warning")},
    error = function(e){message("error")}
)  

