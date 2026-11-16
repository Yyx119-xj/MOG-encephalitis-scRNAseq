rm(list = ls())
#setwd("D:/Workspace/R_learning/sc_analysis/pseudotime")

#https://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247488352&idx=2&sn=c08863668cbde5d66fd5d3e592f3aac9&scene=21#wechat_redirect
#http://cole-trapnell-lab.github.io/monocle-release/

#参数获取
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-r", "--rdsfile"),action = "store",type = "character",default = "",
              help = "The Input rds file"),
  make_option(c("-e","--mean_exp"),action = "store",type = "numeric",default = "0.1",
              help = "The threshold of mean exp"),  
  make_option(c("-g", "--genelist"),action = "store",type = "character",default = "",
              help = "The Input genelist file"),
  make_option(c("-o", "--outdir"),action = "store",type = "character",default = "./",
              help = "The outdir,the default is './'")                            
)
opt <- parse_args(OptionParser(option_list = option_list)) ####解释器

#创建文件夹
if(!file.exists(opt$outdir)){
  dir.create(opt$outdir)
}
setwd(opt$outdir) 
work_dir <- getwd() #设定工作目录

##加载包
library(monocle)
library(RColorBrewer)
library(Seurat)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(gtools)
library(reshape2)

#设定一个全局的主题
mytheme <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1,color = "black"),
  axis.text.y = element_text(color = "black"),
  axis.title.x = element_text(color = "black"),  # 设置 x 轴标题颜色
  axis.title.y = element_text(color = "black")   # 设置 y 轴标题颜色
)

#设定颜色
mycolor3 <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#8DD3C7",
  "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", 
  "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99", "#999999",
  "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
  "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#DCF0B9"
)

scPalette <- function (n) 
{
  colorSpace <- mycolor3
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}
# mycolor <- scPalette(10)

# 测试
# opt$rdsfile <- "../subset/2_cells_analysis/cells_analysis.rds"
# opt$mean_exp <- 0.12
# opt$genelist <- "/users/wuyf/public/RBP/human/RBP_list_human.txt"

##读入已经注释的rds文件
scRNAsub <- readRDS(opt$rdsfile)

##读入genelist文件
# gene <- read.table(opt$genelist,header = F,quote = "")
# gene <- as.character(gene$V1)

##获取counts文件
data <- as(as.matrix(scRNAsub@assays$RNA@counts),'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#创建CellDataSet(cds)对象
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
##expressionFamily参数用于指定表达矩阵的数据类型，有几个选项可以选择：
##稀疏矩阵用negbinomial.size();FPKM值用tobit();logFPKM值用gaussianff()

##mycds是Monocle为我们的数据生成的对象，相当于我们在seurat使用的scRNA对象。
##数据导入后需要进行标准化和其他预处理：

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
#mycds <- detectGenes(mycds, min_expr = 0.5) #不用筛选，因为前期的分析已经筛选过了
##Warning messages:`group_by_()` was deprecated in dplyr 0.7.0.
##与seurat把标准化后的表达矩阵保存在对象中不同，
##monocle只保存一些中间结果在对象中，需要用时再用这些中间结果转化。
##经过上面三个函数的计算，mycds对象中多了SizeFactors、Dipersions、num_cells_expressed和num_genes_expressed等信息。

##完成数据导入和预处理后，就可以考虑选择哪些基因代表细胞的发育特征
##Monocle官网教程提供了4个选择方法：
##选择发育差异表达基因
##选择clusters差异表达基因
##选择离散程度高的基因
##自定义发育marker基因

##使用monocle选择的高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= opt$mean_exp & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cat("高变基因数量：",length(disp.genes),"\n")
mycds <- setOrderingFilter(mycds, disp.genes)
p1 <- plot_ordering_genes(mycds)
ggsave("disp_genes_bymonocle.png",p1,width = 5,height = 5)
ggsave("disp_genes_bymonocle.pdf",p1,width = 5,height = 5)

###使用seurat选择的高变基因
#disp.genes <- VariableFeatures(scRNAsub)
#mycds <- setOrderingFilter(mycds, disp.genes)
#p1 <- plot_ordering_genes(mycds)
#ggsave("disp_genes_byseurat.png",p1,width = 5,height = 5)
#ggsave("disp_genes_byseurat.pdf",p1,width = 5,height = 5)

##使用clusters差异表达基因
#disp_table <- FindAllMarkers(scRNAsub)
#disp.genes <- subset(disp_table,p_val_adj<0.05)$gene
#mycds <- setOrderingFilter(mycds, disp.genes)
#p1 <- plot_ordering_genes(mycds)
#ggsave("disp_genes_byseuratdeg.png",p1,width = 5,height = 5)
#ggsave("disp_genes_byseuratdeg.pdf",p1,width = 5,height = 5)


##降维及细胞排序
##使用disp.genes开展后续分析

#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
# mycds <- orderCells(mycds,root_state = NULL)
# mycds <- orderCells(mycds,reverse = T)

#保存mycds变量
save(mycds,file = "mycds.RData")

#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
ggsave("State.pdf", plot = plot1, width = 6, height = 5)
ggsave("State.png", plot = plot1, width = 6, height = 5)
##Cluster轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
ggsave("Cluster.pdf", plot = plot2, width = 6, height = 5)
ggsave("Cluster.png", plot = plot2, width = 6, height = 5)
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime") + scale_color_viridis_c()
ggsave("Pseudotime.pdf", plot = plot3, width = 6, height = 5)
ggsave("Pseudotime.png", plot = plot3, width = 6, height = 5)
##Celltype轨迹图，细胞注释结果
plot4 <- plot_cell_trajectory(mycds, color_by = "customclassif")
ggsave("Celltype.pdf", plot = plot4, width = 6, height = 5)
ggsave("Celltype.png", plot = plot4, width = 6, height = 5)
##Celltype树形图
plot5 <- plot_complex_cell_trajectory(mycds, x = 1, y = 2,color_by = "customclassif")
ggsave("Celltype_tree.pdf", plot = plot5, width = 6, height = 5)
ggsave("Celltype_tree.png", plot = plot5, width = 6, height = 5)
##Group轨迹图，细胞注释结果
plot6 <- plot_cell_trajectory(mycds, color_by = "group")
ggsave("Group.pdf", plot = plot6, width = 6, height = 5)
ggsave("Group.png", plot = plot6, width = 6, height = 5)


##合并作图
plotc <- plot1|plot2|plot3
ggsave("Combination1.pdf", plot = plotc, width = 10, height = 3.5)
ggsave("Combination1.png", plot = plotc, width = 10, height = 3.5)
plotc <- plot4|plot5
ggsave("Combination2.pdf", plot = plotc, width = 10, height = 5)
ggsave("Combination2.png", plot = plotc, width = 10, height = 5)
##保存结果
write.csv(pData(mycds), "pseudotime.csv")


##轨迹图分面显示
p1 <- plot_cell_trajectory(mycds, color_by = "State") + facet_wrap(~State, nrow = 1)
p2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters") + facet_wrap(~seurat_clusters, nrow = 1)
p3 <- plot_cell_trajectory(mycds, color_by = "group") + facet_wrap(~group, nrow = 1)

ggsave("trajectory_facet_byState.png", plot = p1, width = 8, height = 5)
ggsave("trajectory_facet_byState.pdf", plot = p1, width = 8, height = 5)

ggsave("trajectory_facet_byseurat_clusters.png", plot = p2, width = 8, height = 5)
ggsave("trajectory_facet_byseurat_clusters.pdf", plot = p2, width = 8, height = 5)

ggsave("trajectory_facet_bygroup.png", plot = p3, width = 8, height = 3)
ggsave("trajectory_facet_bygroup.pdf", plot = p3, width = 8, height = 3)

##细胞密度图
df <- pData(mycds) 
## pData(cds)取出的是cds对象中cds@phenoData@data的内容
##pseudotime state density
p1 <- ggplot(df, aes(Pseudotime, colour = State, fill=State)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
ggsave("pseudotime_density_state.png", plot = p1, width = 6, height = 5)
ggsave("pseudotime_density_state.pdf", plot = p1, width = 6, height = 5)
##pseudotime seurat_clusters density
p2 <- ggplot(df, aes(Pseudotime, colour = seurat_clusters, fill=seurat_clusters)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
ggsave("pseudotime_density_seurat_clusters.png", plot = p2, width = 6, height = 5)
ggsave("pseudotime_density_seurat_clusters.pdf", plot = p2, width = 6, height = 5)

##拟时相关基因聚类热图
##
#高变基因
#disp_table <- dispersionTable(mycds)
#disp.genes <- subset(disp_table, mean_expression >= 0.1&dispersion_empirical >= 1*dispersion_fit)
#mean_expression >= 1 这里可以根据情况设置
#disp.genes <- as.character(disp.genes$gene_id)
#差异基因bygroup
# diff_test <- differentialGeneTest(mycds[disp.genes,], 
#                                   fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_bygroup <- differentialGeneTest(mycds[disp.genes,],
                                fullModelFormulaStr = "~group")
diff_test_bygroup <- subset(diff_test_bygroup,qval<1e-04)
diff_test_bygroup <- diff_test_bygroup[order(diff_test_bygroup$qval),]
write.table(diff_test_bygroup,"diff_genes_bygroup.xls",quote = F,row.names = F,sep = "\t")

# target_diff_gene <- diff_test[row.names(diff_test) %in% gene,]
# write.table(target_diff_gene,"diff_genes_bygroup_rbp.xls",quote=F,row.names=F,col.names=T,sep="\t")


#差异基因bypseudotime
diff_test <- differentialGeneTest(mycds[disp.genes,], 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
#diff_test <- differentialGeneTest(mycds[disp.genes,],
#                                 fullModelFormulaStr = "~group")
diff_test <- subset(diff_test,qval<1e-04)
diff_test <- diff_test[order(diff_test$qval),]
write.table(diff_test,"diff_genes_bypseudotime.xls",quote = F,row.names = F,sep = "\t")

#使用差异基因做cluster结果,cluster方法是kmeans
sig_gene_names <- row.names(diff_test)
p1 <- plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=4,
                              show_rownames=F, return_heatmap=T)
ggsave("pseudotime_heatmap_all.png", plot = p1, width = 5, height = 10)
ggsave("pseudotime_heatmap_all.pdf", plot = p1, width = 5, height = 10)

## 输出结果cluster结果,cluster方法是kmeans
clusters = cutree(p1$tree_row,k=4)
clustering = data.frame(clusters)
clustering[,2] = clustering[,1]
clustering[,1]  = rownames(clustering)
clustering = clustering[order(clustering[,2]),]
colnames(clustering)[1]="Symbol"
colnames(clustering)[2]="Cluster"
write.table(clustering,"cluster_list.xls",quote=F,row.names=F,col.names=T,sep="\t")

#取top100差异基因做cluster并画heatmap，cluster num 默认是4
sig_gene_names <- row.names(diff_test) %>% head(100) 
p1 <- plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=4,
                             show_rownames=T, return_heatmap=T)
ggsave("pseudotime_heatmap_top100.png", plot = p1, width = 5, height = 10)
ggsave("pseudotime_heatmap_top100.pdf", plot = p1, width = 5, height = 10)

##基因可视化,选择diff top10 qval最小
top10_diff_genes <- head(sig_gene_names,n = 10)

p1 <- plot_genes_violin(mycds[top10_diff_genes,], grouping = "State", color_by = "State")
ggsave("top10_diff_genes_state_violin.png", plot = p1, width = 3, height = 10)
ggsave("top10_diff_genes_state_violin.pdf", plot = p1, width = 3, height = 10)

p2 <- plot_genes_in_pseudotime(mycds[top10_diff_genes,], color_by = "State")
ggsave("top10_diff_genes_state_pseudotime.png", plot = p2, width = 3, height = 10)
ggsave("top10_diff_genes_state_pseudotime.pdf", plot = p2, width = 3, height = 10)

p3 <- plot_genes_in_pseudotime(mycds[top10_diff_genes,], color_by = "seurat_clusters")
ggsave("top10_diff_genes_cluster_pseudotime.png", plot = p3, width = 3, height = 10)
ggsave("top10_diff_genes_cluster_pseudotime.pdf", plot = p3, width = 3, height = 10)

plotc <- p1|p2|p3
ggsave("top10_diff_genes_combination.png", plot = plotc, width = 8, height = 10)
ggsave("top10_diff_genes_combination.pdf", plot = plotc, width = 8, height = 10)

##使用输入基因画heatmap图，cluster num 默认是4
# sig_gene_names_rbp <- row.names(diff_test)[row.names(diff_test) %in% gene] %>% head(50)
# target_diff_gene <- diff_test[row.names(diff_test) %in% gene,]
# write.table(target_diff_gene,"diff_genes_bypseudotime_rbp.xls",quote=F,row.names=F,col.names=T,sep="\t")

# p1 <- plot_pseudotime_heatmap(mycds[sig_gene_names_rbp,], num_clusters=4,
#                              show_rownames=T, return_heatmap=T)
# ggsave("pseudotime_heatmap_rbp50.png", plot = p1, width = 5, height = 10)
# ggsave("pseudotime_heatmap_rbp50.pdf", plot = p1, width = 5, height = 10)

# ##基因可视化,选择diff top10 qval最小的输入基因
# top10_diff_genes_rbp <- head(sig_gene_names_rbp,n = 10)

# p1 <- plot_genes_violin(mycds[top10_diff_genes_rbp,], grouping = "State", color_by = "State")
# ggsave("top10_diff_genes_rbp_state_violin.png", plot = p1, width = 3, height = 10)
# ggsave("top10_diff_genes_rbp_state_violin.pdf", plot = p1, width = 3, height = 10)

# p2 <- plot_genes_in_pseudotime(mycds[top10_diff_genes_rbp,], color_by = "State")
# ggsave("top10_diff_genes_rbp_state_pseudotime.png", plot = p2, width = 3, height = 10)
# ggsave("top10_diff_genes_rbp_state_pseudotime.pdf", plot = p2, width = 3, height = 10)

# p3 <- plot_genes_in_pseudotime(mycds[top10_diff_genes_rbp,], color_by = "seurat_clusters")
# ggsave("top10_diff_genes_rbp_cluster_pseudotime.png", plot = p3, width = 3, height = 10)
# ggsave("top10_diff_genes_rbp_cluster_pseudotime.pdf", plot = p3, width = 3, height = 10)

# plotc <- p1|p2|p3
# ggsave("top10_diff_genes_rbp_combination.png", plot = plotc, width = 8, height = 10)
# ggsave("top10_diff_genes_rbp_combination.pdf", plot = plotc, width = 8, height = 10)

##BEAM分析
# disp_table <- dispersionTable(mycds)
# disp.genes <- subset(disp_table, mean_expression >= 0.2&dispersion_empirical >= 1*dispersion_fit)
# disp.genes <- as.character(disp.genes$gene_id)
mycds_sub <- mycds[disp.genes,]
#plot_cell_trajectory(mycds_sub, color_by = "State")
beam_res <- BEAM(mycds_sub, branch_point = 1,progenitor_method = "duplicate")
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
write.table(beam_res,"beam_result.xls",quote = F,row.names = F,sep = "\t")

#使用qval<0.01基因画heatmap，并提取cluster基因
sig_beam_gene <- filter(beam_res,qval<0.01) %>% row.names()
mycds_sub_beam <- mycds_sub[sig_beam_gene,]
p1<-plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters = 3, show_rownames = F,return_heatmap = T)
ggsave("beam_heatmap_all.png", plot = p1$ph_res, width = 8, height = 10)
ggsave("beam_heatmap_all.pdf", plot = p1$ph_res, width = 8, height = 10)

## 输出结果cluster结果,cluster方法是kmeans
clusters = cutree(p1$ph_res$tree_row,k=3)
clustering = data.frame(clusters)
clustering[,2] = clustering[,1]
clustering[,1]  = rownames(clustering)
clustering = clustering[order(clustering[,2]),]
colnames(clustering)[1]="Symbol"
colnames(clustering)[2]="Cluster"
write.table(clustering,"beam_result_cluster.xls",quote=F,row.names=F,col.names=T,sep="\t")

##取top100画图
top100_beam_gene <- top_n(beam_res,-100,qval) %>% row.names()
mycds_sub_beam <- mycds_sub[top100_beam_gene,]
p1<-plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters = 3, show_rownames = T,return_heatmap = T)
ggsave("beam_heatmap_top100.png", plot = p1$ph_res, width = 8, height = 10)
ggsave("beam_heatmap_top100.pdf", plot = p1$ph_res, width = 8, height = 10)

# ##取target_gene 画图
# target_beam_gene <- top_n(beam_res[row.names(beam_res) %in% gene,],-50,qval) %>% row.names()
# mycds_sub_beam <- mycds_sub[target_beam_gene,]
# p1<-plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters = 4, show_rownames = T,return_heatmap = T)
# ggsave("beam_heatmap_rbp50.png", plot = p1$ph_res, width = 8, height = 10)
# ggsave("beam_heatmap_rbp50.pdf", plot = p1$ph_res, width = 8, height = 10)

#给seurat对象增加"Size_Factor"，"Pseudotime"，"State"三个属性
scRNAsub@meta.data <- mycds@phenoData@data
saveRDS(scRNAsub,"scRNA_anno.rds")
#导出增加meta信息
#mymeta <- select(scRNAsub@meta.data,replicate,group,seurat_clusters,customclassif,cellcluster,cellgroup,Pseudotime,State)
#cat('Cell\t',file = paste('mymeta','.txt',sep=''))
#write.table(mymeta, file = "mymeta.txt", quote = F, sep = '\t', col.names = TRUE,row.names = T,append = T)

tryCatch(
  {
    # p1 <- FeaturePlot(scRNAsub,reduction = "umap",features = c("Pseudotime"),cols = c("#0000ff","#ff0000")) & coord_fixed() & theme(aspect.ratio = 1)
    p1 <- FeaturePlot(scRNAsub,reduction = "umap",features = c("Pseudotime"),cols = c("#450c59","#f6e629")) & coord_fixed() & theme(aspect.ratio = 1)
    ggsave("Pseudotime_umap.png",width = 6,height = 6)
    ggsave("Pseudotime_umap.pdf",width = 6,height = 6)    
  },
  error = function(e){message("FeaturePlot error")}
)

############################################################################################################
#细胞比例分析 by State
############################################################################################################

meta_all <- dplyr::select(scRNAsub@meta.data,replicate,group,State) 
# meta_all <- scRNAsub@meta.data[,colnames(scRNAsub@meta.data) %in% c("replicate","group","State")]
meta_all$State <- as.numeric(as.character(meta_all$State)) #因子转数值

stat_cluster_fraction_by_sample <- meta_all %>% 
  group_by(replicate,group,State) %>% 
  summarise(count = n()) %>% 
  mutate(percentage = round(count/sum(count)*100,digits = 2)) %>% as.data.frame()

write.table(stat_cluster_fraction_by_sample, file = "State_fraction_by_sample.xls", 
            quote = F, sep = '\t', col.names = T,row.names = F)

#作图
#排序并转成因子
stat_cluster_fraction_by_sample$State <- factor(stat_cluster_fraction_by_sample$State,
                                                          levels = stat_cluster_fraction_by_sample$State %>% unique() %>% mixedsort())
#barplot bygroup
p1 <- ggbarplot(stat_cluster_fraction_by_sample, x = "group",y = "percentage",
                fill = "group",facet.by = "State",
                add = "mean_se",error.plot = "upper_errorbar",
                palette = "npg",ggtheme = theme_classic())
p2 <- facet(p1, facet.by = "State",nrow = 1,scales = "fixed",strip.position = "bottom") +
  theme(legend.position="right",
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        #axis.text.x=element_text(color = "black",size = 12),
        axis.text.y=element_text(color = "black",size = 12)) + 
  scale_y_continuous(expand = c(0,0))# y轴从0開始
ggsave("State_percentage_barplot_bygroup.png",p2,height = 4,width = 9,dpi = 300)
ggsave("State_percentage_barplot_bygroup.pdf",p2,height = 4,width = 9)

#barplot bysample
p1 <- ggbarplot(
  data = stat_cluster_fraction_by_sample,
  x = "State",
  y = "percentage",
  fill = "replicate",
  position = position_dodge(),
  palette = "npg",
  ggtheme = theme_classic()
) +
  theme(legend.position="right",
        #axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.text.x=element_text(color = "black",size = 12),
        axis.text.y=element_text(color = "black",size = 12)) +
  scale_y_continuous(expand = c(0,0))# y轴从0開始
ggsave("State_percentage_barplot_byreplicate.png",p1,height = 5,width = 12,dpi = 300)
ggsave("State_percentage_barplot_byreplicate.pdf",p1,height = 5,width = 12)

############################################################################################################
#取top30的差异基因画一下curve图
############################################################################################################

genes <- row.names(diff_test_bygroup) %>% head(30)

genes_exp <- list()

for(i in 1:length(genes)){
  A <- log2(exprs(mycds)[genes[i],]+1)
  A <- as.data.frame(A)
  genes_exp[[i]] <- A
}

gene_exp <- do.call(cbind, genes_exp)
colnames(gene_exp) <- genes
#将上述几个基因的拟时表达添加到monocle
pData(mycds) = cbind(pData(mycds), gene_exp)
#提取作图数据，只需要有基因表达和拟时即可
data <- pData(mycds)
colnames(data)
#选择需要的列即可，我这里的origin.ident就是分组
data <- data %>% select("group","Pseudotime":ncol(data)) %>% select(-"State")

#ggplot作图
#使用分屏，应该就是文献种的办法
#首先将data宽数据转化为长数据
data_long_m<-melt(data, id.vars = c("group", "Pseudotime"), #需保留的不参与聚合的变量列名
                  measure.vars = 3:ncol(data),#选择需要转化的列
                  variable.name = c('gene'),#聚合变量的新列名
                  value.name = 'value')#聚合值的新列名
colnames(data_long_m)

p1 <- ggplot(data_long_m, aes(x=Pseudotime, y=value, color=group))+
  geom_smooth(aes(fill=group))+ #平滑的填充
  xlab('pseudotime') + 
  ylab('Relative expression') +
  facet_wrap(~gene, scales = "free_y", ncol = 6)+ #分面，y轴用各自数据
  theme(axis.text = element_text(color = 'black',size = 12),
        axis.title = element_text(color = 'black',size = 14),
        strip.text = element_text(color = 'black',size = 14)) #分面标题
 # scale_color_manual(name=NULL, values = c("#089E86","#3D5387"))+#修改颜色
 # scale_fill_manual(name=NULL, values = c("#089E86","#3D5387"))#修改颜色

ggsave("top30_diff_genes_curve_plot.png",p1,width=18,height=15,dpi=300)
ggsave("top30_diff_genes_curve_plot.pdf",p1,width=18,height=15)


