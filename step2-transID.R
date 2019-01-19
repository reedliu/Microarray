rm(list=ls())
load("GSE62832.eSet.Rdata")

library(GEOquery)
e <- exprs(eSet[[1]])
####################################
# 3.将表达矩阵的探针ID转换为gene ID
####################################
#3.1首先需要知道GSE62832对应平台是GPL6244 
# eSet

#3.2然后需要知道GPL6244对应哪个注释包
#先获取GPL平台和Bioc注释包的对应关系【https://support.bioconductor.org/p/22784/】
#或者直接在这里寻找【https://www.jianshu.com/p/f6906ba703a0】
#下载注释包
if(F){source("https://bioconductor.org/biocLite.R")
  options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
  biocLite("hugene10sttranscriptcluster.db")
}
#3.3一切就绪，开始探索、过滤、整合
library(hugene10sttranscriptcluster.db) 
ls("package:hugene10sttranscriptcluster.db") #查看所有包含的对象
s <- toTable(hugene10sttranscriptclusterSYMBOL) #将SYMBOL对象转为数据框
#【补充：这个注释包是别人上传的，我们只是拿来用。其中的表达矩阵中原来有3w多探针，作者过滤后只留下接近2w的探针。这仅仅是作者过滤后的结果，并非原始数据结果】
#一些简单的探索#
library(magrittr)
if(F){
  s$symbol %>% unique() %>% length() #看一下有多少基因【人类正常蛋白编码基因就2w左右】
  s$symbol %>% table() %>% sort() %>% tail() #基因使用探针最多的前6名【发现有两个基因用了10个探针】
  s$symbol %>% table() %>% sort() %>% table() #基因与探针对应关系
}
# 过滤+整合
if(T){
  dim(e) #过滤前
  #过滤【只保留和注释文件探针id相同的探针】
  efilt <- e[rownames(e)%in%s$probe_id,]
  dim(efilt)#过滤后
  maxp = by(efilt,s$symbol,function(x) rownames(x)[which.max(rowMeans(x))]) 
  uniprobes = as.character(maxp)
  efilt=efilt[rownames(efilt)%in%uniprobes,]
  #整合2【目的：将我们表达矩阵的行名换成刚才一对一的基因名，并且match这个函数保证了表达矩阵和注释包的顺序是一致的】
  rownames(efilt)=s[match(rownames(efilt),s$probe_id),2]
}

####################################
# 4.对表达矩阵进行一些检验
####################################
#表达矩阵不局限于GEO数据库的芯片分析，转录组及其他涉及基因、样本、分组关系的都会有一个表达量矩阵，就是一个基因在不同样本中（对照、处理；是否患病等）的表达差异。拿到表达矩阵，在进行后续分析之前，首先要检测这个矩阵是不是合理的，比如看管家基因是否突出、一致；样本分组是不是和实验设计一致，用PCA、hclust检验
#4.1 检测一些管家基因表达量
boxplot(efilt[,1])#看看第一个样本中总体基因表达量分布，可以看到基本为5左右
efilt['ACTB',] #激动蛋白Beta-actin的基因名是ACTB，管家基因
efilt['GAPDH',] #也是管家基因
#4.2 看表达矩阵的整体分布
#先把表达矩阵=》tidy data【四列：基因名、样本、表达量、表型分组(看文献按MAO、MNO分组)】
if(T){
  library(reshape2)
  pdata=pData(eSet[[1]]) #将样本表型信息从数据框中提取出来【取出来的是表型、样本的数据框】
  group_list=as.character(pdata$`metabolic status:ch1`) 
  m_efilt = melt(efilt) #先将原来矩阵“融化
  colnames(m_efilt)=c('symbol','sample','value') #重新命名三列
  m_efilt$group=rep(group_list,each=nrow(efilt)) 
}
#boxplot
if(T){
  library(ggpubr)
  library(ggplot2)
  
  ggplot(data = m_efilt, aes(x= sample, y = value, fill = group))+
    geom_boxplot()+
    ggtitle("Plot of expression value \n MAO vs. MNO")+
    xlab("Sample") + ylab("Expression value")+
    stat_summary(fun.y="mean",geom="point",shape=24,size=2,fill="orange", alpha =0.5)+
    theme_set(theme_set(theme_bw(base_size=14)))+
    theme(text=element_text(face='plain'),axis.text.x=element_text(angle=60,hjust=1))+
    theme(axis.title.x = element_text(color="black", size=14, face="bold"))+
    theme(axis.title.y = element_text(color="black", size=14, face="bold"))+
    theme(plot.title = element_text(hjust = 0.5))
}
#density
ggplot(m_efilt,aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 4)
#hclust
if(T){
  library(factoextra)
  colnames(efilt) <- paste(group_list,1:36,sep='')
  dd <- dist(scale(t(efilt)), method = "euclidean")
  hc <- hclust(dd, method = "ward.D2")
  fviz_dend(hc, k = 4, # Cut in four groups
            cex = 1, # label size
            k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
            color_labels_by_k = TRUE, # color labels by groups
            rect = TRUE, # Add rectangle around groups
            rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), 
            rect_fill = TRUE,
            horiz = TRUE)
}
#PCA
if(T){
  library(ggfortify)
  df <- as.data.frame(t(efilt))
  df$group <- group_list
  autoplot(prcomp( df[,1:(ncol(df)-1)]), data=df, colour = 'group') 
}
