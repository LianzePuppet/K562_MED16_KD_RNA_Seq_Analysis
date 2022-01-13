setwd("E:/Desktop/k562/code/K562_MED16_KD_RNA_Seq_Analysis/")
library_packages<-function(){
  library('limma')
  library('gplots')
  library('ggplot2')#火山图
  library('edgeR')
  library('reshape2')
  library('DESeq2')
  library('limma')
  library('edgeR')
  library('reshape2')
  library('DESeq2')
  library('BiocGenerics')
  library('xlsx')#GSEA需要的表达矩阵
  library('factoextra')#Cluster Dendrogram
  library('gmodels')
  library('ggpubr')
  library('ggplot2')
  library("RColorBrewer")
  library('EnhancedVolcano')
  library('airway')
  library('pheatmap')
  library('edgeR')
  library('limma')
  library(ggthemes)
  library(ggplot2)
  library(ggridges)
}
library_packages()

###############Step1:读取数据
exprSet <- read.csv("../../data/readCount.K562.csv")
read_data<-function(exprSet){
  a1=exprSet
  ids=a1[,1:2] # 第1，2列是基因名字
  dat=a1[, 3:6] #选取的是样本列
  rownames(dat)=a1$geneID#更该行名为gene symbol
  #head(ids)#是geneID   geneName的对应表格
  colnames(ids)=c('probe_id','symbol')#修改列名为探针ID和symbol
  #dat[1:4,1:4]
  dat=dat[ids$probe_id,]
  ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
  ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
  dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
  rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
  #dat[1:4,1:4] #保留每个基因ID第一次出现的信息
  #去重完毕
  #综上，当有probe对应多个gene symbol时，只保留表达量中位数最大的那一个，并成功将列名保留为symbol
  return (dat)
}
dat<-read_data(exprSet)

##################GSEA需要的表达矩阵
filename1="./data/exprSet.xls"
save_GSEA<-function(save_gct,filename){
  DESCRIPTION<-rownames(save_gct)
  save_gct<-data.frame(DESCRIPTION,save_gct)
  write.xlsx(save_gct,file=filename)#保存GSEA所需要的表达矩阵
}
#save_GSEA(exprSet2,filename1)


###########################
###  Then  for DESeq2  ####
###########################
group_list<-c('untreat','untreat','treat','treat')#分组矩阵
get_dds<-function(exprSet2,group_list){
  suppressMessages(library(DESeq2))
  #construct dds matrix and save as Rdata objects
  (colData <- data.frame(row.names=colnames(exprSet2), group_list=group_list) )
  dds <- DESeqDataSetFromMatrix(countData = exprSet2,
                                colData = colData,
                                design = ~ group_list)
  #构建dds对象，需要一个表达矩阵和分组矩阵
  dds <- DESeq(dds)
  return(dds)
}
dds<-get_dds(dat,group_list)




#############################
###   sample-to-sample    ###
###   sample-to-sample    ###
###          PCA          ###
#############################
paint_sample_to_sample<-function(dds,exprSet2){
  vsdata <- vst(dds, blind=FALSE)
  #用vst来标化数据，实际上还有rlog方法，或者就是log2的方法，官网推荐< 30个样本用rlog，大于30个样本用vst，速度更快，这里我们不要计较那么多了，就用vst，因为真实的TCGA数据，样本往往大于30个。
  sampleDists <- dist(t(assay(vsdata)))
  
  #sample-to-sample& 的距离矩阵的热图
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsdata$group_list,colnames(dat), sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  pdf("./figure/sample_to_sample_distances_heatmap.pdf")
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors,
           show_colnames = T,)
  dev.off()
  return(sampleDists)
}
paint_sample_to_sample(dds,dat)
#使用dist函数计算样本间的距离并用hclust进行层次聚类
#这个图说明两个样本分开了不需要特殊处理
vsdata <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsdata )))
res1<- hcut(sampleDists, k = 2, stand = TRUE)
# Visualize
pdf("./figure/Cluster Dendrogram.pdf", pointsize=20)
fviz_dend(res1, rect = TRUE,rect_border="cluster",rect_lty=2,lwd=1.2,rect_fill = T,cex = 1.5,color_labels_by_k=T,horiz=T)
dev.off()
#PCA
pca.info <- fast.prcomp(dat)
#显示PCA计算结果
head(pca.info$rotation)
#计算treat和untreat之间的差别
pca.data <- data.frame(sample = rownames(pca.info$rotation),Type = c(rep("untreat",2),rep("treat",2)),pca.info$rotation)
#绘图
pca.data $name<-colnames(dat)
pdf("./figure/PCA.pdf", pointsize=20)
ggscatter(pca.data,x = "PC1",y = "PC2",color = "Type",label = "name") + theme_base()
dev.off()



#得到DEG
res <- results(dds, contrast=c("group_list","treat","untreat"))#差异分析结果的对象
res <- res[!is.na(res$padj),]
get_DEG<-function(res){
  res <- results(dds, contrast=c("group_list","treat","untreat"))#差异分析结果的对象
  res <- res[!is.na(res$padj),]
  resOrdered <- res[order(res$padj),]#是一个dataform的对象
  head(resOrdered)
  DEG=as.data.frame(resOrdered)
  DEG<-na.omit(DEG)   #这里不确定要不要去掉NA值
  return(DEG)
}
DEG <-get_DEG(res)

#设定logFC_cutoff
logFC_cutoff<-log2(1.2)



#draw_MA_plot
#对log2 fold change进行收缩一下，得到的MA图会好看一些。
res_order<-res[order(row.names(res)),]
res1 = res_order
resApe <-lfcShrink(dds, coef=2,type="apeglm")
pdf("./figure/MA_plot.pdf")
plotMA(resApe, ylim = c(-5,5))
topGene <- rownames(res1)[which.min(res1$padj)]
with(res1[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
dev.off()

#计算上下调标准
#save_DESeq2_DEG
save_DESeq2_DEG<-function(DEG){
  DEG$change = as.factor(
    ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
           ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  head(DEG)
  table(DEG$change)
  DESeq2_DEG <- DEG
  write.csv(DESeq2_DEG,file = "./data/DESeq2_DEG.csv")
  return(DESeq2_DEG)
}


this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe nxumber of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',]))

DESeq2_DEG<-save_DESeq2_DEG(DEG)

#火山图绘制
pdf("./figure/_volcano2.pdf")
EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue',title=this_tile)
dev.off()
#X轴是变化倍数，Y轴是p.value




###########################
###  Then  for edegR   ####
###########################
dge <- DGEList(counts=dat,group=group_list)
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group_list)
rownames(design)<-colnames(dge)
colnames(design)<-levels(group_list)
dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
fit2 <- glmLRT(fit, contrast=c(1,-1))
DEG=topTags(fit2, n=nrow(dat))
DEG=data.frame(DEG)
k1 = (DEG$PValue < 0.05)&(DEG$logFC < -logFC_cutoff)
k2 = (DEG$PValue < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
head(DEG)
table(DEG$change)
edgeR_DEG <- DEG
write.csv(edgeR_DEG,file = "./data/edgeR_DEG.csv")
#绘制分类图
draw_classification<-function(exprSet){
  #edgeR使用LogFC 观察样本分类情况
  ##载入数据及生成DGEList
  d <- DGEList(counts=dat,group=factor(group_list))#construct an object
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)#Normalization
  d$samples
  degs=d
  
  # cpm normalization
  countsPerMillion <- cpm(degs)
  countCheck <- countsPerMillion > 1
  keep <- which(rowSums(countCheck) >= 2)
  degs.keep <- degs[keep,]
  dim(degs.keep)
  
  #通过logFC查看样本的分组情况
  degs.norm <- calcNormFactors(degs.keep, method = 'TMM')
  pdf("./figure/classification.pdf")
  plotMDS(degs.norm, col=as.numeric(degs.norm$samples$group))
  legend("left",as.character(unique(degs.norm$samples$group)), col=1:3, pch=20)
  dev.off()
}
draw_classification(dat)
#############################




###########################
#  Then  for limma/voom   #
###########################
library(limma)
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(dat)
deg <- DGEList(counts=dat)
deg <- calcNormFactors(deg)
logCPM <- cpm(deg, log=TRUE, prior.count=3)
pdf("./figure/limma_linear_model.pdf", pointsize=20)
v <- voom(deg,design,plot=TRUE, normalize="quantile")
dev.off()

fit <- lmFit(v, design)
group_list
cont.matrix=makeContrasts(contrasts=c('treat-untreat'),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)
tempOutput=topTable(fit2,coef='treat-untreat',n=Inf)
DEG_limmma_voom=na.omit(tempOutput)
logFC_cutoff <- with(DEG_limmma_voom,mean(abs(logFC)) + 2*sd(abs(logFC)) )
logFC_cutoff=0.7613179
k1 = (DEG_limmma_voom$P.Value < 0.05)&(DEG_limmma_voom$logFC < -logFC_cutoff)
k2 = (DEG_limmma_voom$P.Value < 0.05)&(DEG_limmma_voom$logFC > logFC_cutoff)
DEG_limmma_voom$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG_limmma_voom$change)
head(DEG_limmma_voom)
limma_voom_DEG <- DEG_limmma_voom
write.csv(limma_voom_DEG,file = "./data/limma_voom_DEG.csv")
###############



######################################################################
###################     三种方法得到上下调基因统计   #################
######################################################################

#分别得到三款软件计算的上调、下调以及既不上调也不下调的基因数
tj = data.frame(deseq2 = as.integer(table(DESeq2_DEG$change)),
                edgeR = as.integer(table(edgeR_DEG$change)),
                limma_voom = as.integer(table(limma_voom_DEG$change)),
                row.names = c("down","not","up")
);
UP=function(df){
  rownames(df)[df$change=="UP"]
}
DOWN=function(df){
  rownames(df)[df$change=="DOWN"]
}
up = intersect(intersect(UP(DESeq2_DEG),UP(edgeR_DEG)),UP(limma_voom_DEG))
write.csv(up,file="./data/up.csv")
down= intersect(intersect(DOWN(DESeq2_DEG),DOWN(edgeR_DEG)),DOWN(limma_voom_DEG))
write.csv(down,file="./data/down.csv")

#上调基因画维恩图
up_DESeq2<-UP(DESeq2_DEG)
write.csv(up_DESeq2,file="./data/up_DESeq2.csv")
up_edgeR<-UP(edgeR_DEG)
write.csv(up_edgeR,file="./data/up_edgeR.csv")
up_limma_voom<-UP(limma_voom_DEG)
write.csv(up_limma_voom,file="./data/up_limma_voom.csv")

df1<-read.csv("./data/up_limma_voom.csv",header = T,stringsAsFactors = F)
df2<-read.csv("./data/up_edgeR.csv",header = T,stringsAsFactors = F)
df3<-read.csv("./data/up_DESeq2.csv",header = T,stringsAsFactors = F)
library(ggvenn)
pdf("./figure/upregulate_gene.pdf", pointsize=20)
x<-list(limma_voom=df1$x,edgeR=df2$x,DESeq2=df3$x)
ggvenn(x,c("limma_voom","edgeR","DESeq2"),set_name_size = 6,fill_alpha = 1,text_size =4,show_percentage = F,
       stroke_color = "white",
       fill_color = c("#ffb2b2","#b2e7cb","#b2d4ec"),
       set_name_color = c("#ff0000","#4a9b83","#1d6295"))
dev.off()

#下调基因画韦恩图
down_DESeq2<-DOWN(DESeq2_DEG)
write.csv(down_DESeq2,file="./data/down_DESeq2.csv")
down_edgeR<-DOWN(edgeR_DEG)
write.csv(down_edgeR,file="./data/down_edgeR.csv")
down_limma_voom<-DOWN(limma_voom_DEG)
write.csv(down_limma_voom,file="./data/down_limma_voom.csv")
df4<-read.csv("./data/down_limma_voom.csv",header = T,stringsAsFactors = F)
df5<-read.csv("./data/down_edgeR.csv",header = T,stringsAsFactors = F)
df6<-read.csv("./data/down_DESeq2.csv",header = T,stringsAsFactors = F)
library(ggvenn)
x<-list(limma_voom=df4$x,edgeR=df5$x,DESeq2=df6$x)
pdf("./figure/downregulate_gene.pdf",  pointsize=20)
ggvenn(x,c("limma_voom","edgeR","DESeq2"),set_name_size = 6,fill_alpha = 1,text_size =4,show_percentage = F,
       stroke_color = "white",
       fill_color = c("#ffb2b2","#b2e7cb","#b2d4ec"),
       set_name_color = c("#ff0000","#4a9b83","#1d6295"))
dev.off()
#set_name_size是图名称大小，
upregulate_gene_intersect<-intersect(intersect(df1$x,df2$x),df3$x)
downregulate_gene_intersect<-intersect(intersect(df4$x,df5$x),df6$x)
upregulate_gene_union<-union(union(df1$x,df2$x),df3$x)
downregulate_gene_union<-union(union(df4$x,df5$x),df6$x)
write.csv(upregulate_gene_intersect,file="./data/upregulate_gene_intersect.csv")
write.csv(downregulate_gene_intersect,file="./data/downregulate_gene_intersect.csv")
write.csv(upregulate_gene_union,file="./data/upregulate_gene_union.csv")
write.csv(downregulate_gene_union,file="./data/downregulate_gene_union.csv")


setwd("E:/Desktop/k562/code/finalfinal")

library_packages<-function(){
  library('limma')
  library('gplots')
  library('ggplot2')#火山图
  library('edgeR')
  library('reshape2')
  library('DESeq2')
  library('limma')
  library('edgeR')
  library('reshape2')
  library('DESeq2')
  library('BiocGenerics')
  library('xlsx')#GSEA需要的表达矩阵
  library('factoextra')#Cluster Dendrogram
  library('gmodels')
  library('ggpubr')
  library('ggplot2')
  library("RColorBrewer")
  library('EnhancedVolcano')
  library('airway')
  library('pheatmap')
  library('edgeR')
  library('limma')
  library(ggthemes)
  library(ggplot2)
  library(ggridges)
}
library_packages()

###############Step1:读取数据
exprSet <- read.csv("../../data/readCount.K562.csv")
read_data<-function(exprSet){
  a1=exprSet
  ids=a1[,1:2] # 第1，2列是基因名字
  dat=a1[, 3:6] #选取的是样本列
  rownames(dat)=a1$geneID#更该行名为gene symbol
  #head(ids)#是geneID   geneName的对应表格
  colnames(ids)=c('probe_id','symbol')#修改列名为探针ID和symbol
  #dat[1:4,1:4]
  dat=dat[ids$probe_id,]
  ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
  ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
  dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
  rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
  #dat[1:4,1:4] #保留每个基因ID第一次出现的信息
  #去重完毕
  #综上，当有probe对应多个gene symbol时，只保留表达量中位数最大的那一个，并成功将列名保留为symbol
  return (dat)
}
dat<-read_data(exprSet)

##################GSEA需要的表达矩阵
filename1="./data/exprSet.xls"
save_GSEA<-function(save_gct,filename){
  DESCRIPTION<-rownames(save_gct)
  save_gct<-data.frame(DESCRIPTION,save_gct)
  write.xlsx(save_gct,file=filename)#保存GSEA所需要的表达矩阵
}
#save_GSEA(exprSet2,filename1)


###########################
###  Then  for DESeq2  ####
###########################
group_list<-c('untreat','untreat','treat','treat')#分组矩阵
get_dds<-function(exprSet2,group_list){
  suppressMessages(library(DESeq2))
  #construct dds matrix and save as Rdata objects
  (colData <- data.frame(row.names=colnames(exprSet2), group_list=group_list) )
  dds <- DESeqDataSetFromMatrix(countData = exprSet2,
                                colData = colData,
                                design = ~ group_list)
  #构建dds对象，需要一个表达矩阵和分组矩阵
  dds <- DESeq(dds)
  return(dds)
}
dds<-get_dds(dat,group_list)




#############################
###   sample-to-sample    ###
###   sample-to-sample    ###
###          PCA          ###
#############################
paint_sample_to_sample<-function(dds,exprSet2){
  vsdata <- vst(dds, blind=FALSE)
  #用vst来标化数据，实际上还有rlog方法，或者就是log2的方法，官网推荐< 30个样本用rlog，大于30个样本用vst，速度更快，这里我们不要计较那么多了，就用vst，因为真实的TCGA数据，样本往往大于30个。
  sampleDists <- dist(t(assay(vsdata)))
  
  #sample-to-sample& 的距离矩阵的热图
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsdata$group_list,colnames(dat), sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  pdf("./figure/sample_to_sample_distances_heatmap.pdf")
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors,
           show_colnames = T,)
  dev.off()
  return(sampleDists)
}
paint_sample_to_sample(dds,dat)
#使用dist函数计算样本间的距离并用hclust进行层次聚类
#这个图说明两个样本分开了不需要特殊处理
vsdata <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsdata )))
res1<- hcut(sampleDists, k = 2, stand = TRUE)
# Visualize
pdf("./figure/Cluster Dendrogram.pdf", pointsize=20)
fviz_dend(res1, rect = TRUE,rect_border="cluster",rect_lty=2,lwd=1.2,rect_fill = T,cex = 1.5,color_labels_by_k=T,horiz=T)
dev.off()
#PCA
pca.info <- fast.prcomp(dat)
#显示PCA计算结果
head(pca.info$rotation)
#计算treat和untreat之间的差别
pca.data <- data.frame(sample = rownames(pca.info$rotation),Type = c(rep("untreat",2),rep("treat",2)),pca.info$rotation)
#绘图
pca.data $name<-colnames(dat)
pdf("./figure/PCA.pdf", pointsize=20)
ggscatter(pca.data,x = "PC1",y = "PC2",color = "Type",label = "name") + theme_base()
dev.off()



#得到DEG
res <- results(dds, contrast=c("group_list","treat","untreat"))#差异分析结果的对象
res <- res[!is.na(res$padj),]
get_DEG<-function(res){
  res <- results(dds, contrast=c("group_list","treat","untreat"))#差异分析结果的对象
  res <- res[!is.na(res$padj),]
  resOrdered <- res[order(res$padj),]#是一个dataform的对象
  head(resOrdered)
  DEG=as.data.frame(resOrdered)
  DEG<-na.omit(DEG)   #这里不确定要不要去掉NA值
  return(DEG)
}
DEG <-get_DEG(res)

#设定logFC_cutoff
logFC_cutoff<-log2(1.2)



#draw_MA_plot
#对log2 fold change进行收缩一下，得到的MA图会好看一些。
res_order<-res[order(row.names(res)),]
res1 = res_order
resApe <-lfcShrink(dds, coef=2,type="apeglm")
pdf("./figure/MA_plot.pdf")
plotMA(resApe, ylim = c(-5,5))
topGene <- rownames(res1)[which.min(res1$padj)]
with(res1[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
dev.off()

#计算上下调标准
#save_DESeq2_DEG
save_DESeq2_DEG<-function(DEG){
  DEG$change = as.factor(
    ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
           ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  head(DEG)
  table(DEG$change)
  DESeq2_DEG <- DEG
  write.csv(DESeq2_DEG,file = "./data/DESeq2_DEG.csv")
  return(DESeq2_DEG)
}

DESeq2_DEG<-save_DESeq2_DEG(DEG)

this_tile <- paste0('Cutoff for FC is ',round(1.2,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',]))

#火山图绘制
pdf("./figure/_volcano2.pdf")
EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue',title=this_tile)
dev.off()
#X轴是变化倍数，Y轴是p.value




###########################
###  Then  for edegR   ####
###########################
dge <- DGEList(counts=dat,group=group_list)
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group_list)
rownames(design)<-colnames(dge)
colnames(design)<-levels(group_list)
dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
fit2 <- glmLRT(fit, contrast=c(1,-1))
DEG=topTags(fit2, n=nrow(dat))
DEG=data.frame(DEG)
k1 = (DEG$PValue < 0.05)&(DEG$logFC < -logFC_cutoff)
k2 = (DEG$PValue < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
head(DEG)
table(DEG$change)
edgeR_DEG <- DEG
write.csv(edgeR_DEG,file = "./data/edgeR_DEG.csv")
#绘制分类图
draw_classification<-function(exprSet){
  #edgeR使用LogFC 观察样本分类情况
  ##载入数据及生成DGEList
  d <- DGEList(counts=dat,group=factor(group_list))#construct an object
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)#Normalization
  d$samples
  degs=d
  
  # cpm normalization
  countsPerMillion <- cpm(degs)
  countCheck <- countsPerMillion > 1
  keep <- which(rowSums(countCheck) >= 2)
  degs.keep <- degs[keep,]
  dim(degs.keep)
  
  #通过logFC查看样本的分组情况
  degs.norm <- calcNormFactors(degs.keep, method = 'TMM')
  pdf("./figure/classification.pdf")
  plotMDS(degs.norm, col=as.numeric(degs.norm$samples$group))
  legend("left",as.character(unique(degs.norm$samples$group)), col=1:3, pch=20)
  dev.off()
}
draw_classification(dat)
#############################




###########################
#  Then  for limma/voom   #
###########################
library(limma)
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(dat)
deg <- DGEList(counts=dat)
deg <- calcNormFactors(deg)
logCPM <- cpm(deg, log=TRUE, prior.count=3)
pdf("./figure/limma_linear_model.pdf", pointsize=20)
v <- voom(deg,design,plot=TRUE, normalize="quantile")
dev.off()

fit <- lmFit(v, design)
group_list
cont.matrix=makeContrasts(contrasts=c('treat-untreat'),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)
tempOutput=topTable(fit2,coef='treat-untreat',n=Inf)
DEG_limmma_voom=na.omit(tempOutput)
logFC_cutoff <- with(DEG_limmma_voom,mean(abs(logFC)) + 2*sd(abs(logFC)) )
logFC_cutoff=0.7613179
k1 = (DEG_limmma_voom$P.Value < 0.05)&(DEG_limmma_voom$logFC < -logFC_cutoff)
k2 = (DEG_limmma_voom$P.Value < 0.05)&(DEG_limmma_voom$logFC > logFC_cutoff)
DEG_limmma_voom$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG_limmma_voom$change)
head(DEG_limmma_voom)
limma_voom_DEG <- DEG_limmma_voom
write.csv(limma_voom_DEG,file = "./data/limma_voom_DEG.csv")
###############



######################################################################
###################     三种方法得到上下调基因统计   #################
######################################################################

#分别得到三款软件计算的上调、下调以及既不上调也不下调的基因数
tj = data.frame(deseq2 = as.integer(table(DESeq2_DEG$change)),
                edgeR = as.integer(table(edgeR_DEG$change)),
                limma_voom = as.integer(table(limma_voom_DEG$change)),
                row.names = c("down","not","up")
);
UP=function(df){
  rownames(df)[df$change=="UP"]
}
DOWN=function(df){
  rownames(df)[df$change=="DOWN"]
}
up = intersect(intersect(UP(DESeq2_DEG),UP(edgeR_DEG)),UP(limma_voom_DEG))
write.csv(up,file="./data/up.csv")
down= intersect(intersect(DOWN(DESeq2_DEG),DOWN(edgeR_DEG)),DOWN(limma_voom_DEG))
write.csv(down,file="./data/down.csv")

#上调基因画维恩图
up_DESeq2<-UP(DESeq2_DEG)
write.csv(up_DESeq2,file="./data/up_DESeq2.csv")
up_edgeR<-UP(edgeR_DEG)
write.csv(up_edgeR,file="./data/up_edgeR.csv")
up_limma_voom<-UP(limma_voom_DEG)
write.csv(up_limma_voom,file="./data/up_limma_voom.csv")

df1<-read.csv("./data/up_limma_voom.csv",header = T,stringsAsFactors = F)
df2<-read.csv("./data/up_edgeR.csv",header = T,stringsAsFactors = F)
df3<-read.csv("./data/up_DESeq2.csv",header = T,stringsAsFactors = F)
library(ggvenn)
pdf("./figure/upregulate_gene.pdf", pointsize=20)
x<-list(limma_voom=df1$x,edgeR=df2$x,DESeq2=df3$x)
ggvenn(x,c("limma_voom","edgeR","DESeq2"),set_name_size = 6,fill_alpha = 1,text_size =4,show_percentage = F,
       stroke_color = "white",
       fill_color = c("#ffb2b2","#b2e7cb","#b2d4ec"),
       set_name_color = c("#ff0000","#4a9b83","#1d6295"))
dev.off()

#下调基因画韦恩图
down_DESeq2<-DOWN(DESeq2_DEG)
write.csv(down_DESeq2,file="./data/down_DESeq2.csv")
down_edgeR<-DOWN(edgeR_DEG)
write.csv(down_edgeR,file="./data/down_edgeR.csv")
down_limma_voom<-DOWN(limma_voom_DEG)
write.csv(down_limma_voom,file="./data/down_limma_voom.csv")
df4<-read.csv("./data/down_limma_voom.csv",header = T,stringsAsFactors = F)
df5<-read.csv("./data/down_edgeR.csv",header = T,stringsAsFactors = F)
df6<-read.csv("./data/down_DESeq2.csv",header = T,stringsAsFactors = F)
library(ggvenn)
x<-list(limma_voom=df4$x,edgeR=df5$x,DESeq2=df6$x)
pdf("./figure/downregulate_gene.pdf",  pointsize=20)
ggvenn(x,c("limma_voom","edgeR","DESeq2"),set_name_size = 6,fill_alpha = 1,text_size =4,show_percentage = F,
       stroke_color = "white",
       fill_color = c("#ffb2b2","#b2e7cb","#b2d4ec"),
       set_name_color = c("#ff0000","#4a9b83","#1d6295"))
dev.off()
#set_name_size是图名称大小，
upregulate_gene_intersect<-intersect(intersect(df1$x,df2$x),df3$x)
downregulate_gene_intersect<-intersect(intersect(df4$x,df5$x),df6$x)
upregulate_gene_union<-union(union(df1$x,df2$x),df3$x)
downregulate_gene_union<-union(union(df4$x,df5$x),df6$x)
write.csv(upregulate_gene_intersect,file="./data/upregulate_gene_intersect.csv")
write.csv(downregulate_gene_intersect,file="./data/downregulate_gene_intersect.csv")
write.csv(upregulate_gene_union,file="./data/upregulate_gene_union.csv")
write.csv(downregulate_gene_union,file="./data/downregulate_gene_union.csv")


