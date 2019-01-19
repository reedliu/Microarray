#################################
####### 芯片数据分析 ###########
# made by Reedliu 18.9.10#####
#################################
rm(list=ls())
# prepare packages
if(F){
  source("https://bioconductor.org/biocLite.R")
  options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
  biocLite(c("GEOquery","limma","DESeq2", "clusterProfiler")) #软件包
  install.packages(c("reshape2","tidyverse")) #工具包
  install.packages("ggplot2") #绘图包
}
#################################
# 1.下载数据
#################################
#方法一：自己下载【自己定义函数】#自己再修改下
if(F){
  'get_GSE_links(x)' <- 
  function(studyID = x, down = F, destdir = "./") {
    ## studyID destdir ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62832/matrix/
    matrix_link = paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/", substr(studyID, 1, nchar(studyID) - 3), "nnn/", studyID, "/matrix/", studyID, "_series_matrix.txt.gz")
    print(paste0("The URL for matrix   is : ", matrix_link))
  }
get_GSE_links('GSE62832')
}
#方法二：用GEOquery下载
if(F){
  library(GEOquery)
  eSet <- getGEO('GSE62832', destdir = '.', getGPL = F, AnnotGPL = F)
  save(eSet, file = 'GSE62832.eSet.Rdata')
}

#################################
# 2.将下载的表达矩阵读进来
#################################
#方法一：针对自己下载好的矩阵文件
if(F){
e1 <- read.csv('GSE62832_series_matrix.txt.gz',
               comment.char = '!', fill = T, sep = '\t')
#再将e1的行名转为ID号，也就是把目前的第一列复制行名的位置，然后再把第一列去掉
rownames(e1) <- e1[,1]
e <- e1[,-1]
}
#方法二：针对GEOquery下载的eSet
load('GSE62832.eSet.Rdata')
library(GEOquery)
class(eSet) #先查看数据种类【列表还是数据框】
str(eSet) #再看数据结构【可以看到第二行涉及到了表达矩阵】
#既然是列表，那么从列表中提取信息就是 eSet[[1]]
e <- exprs(eSet[[1]])

#################################
# 附加.一键下载+得到表达矩阵
#################################
if(T){
  'prepGSE(x)' <- function(GEO = x, dir = ".") {
    library(GEOquery)
    eSet <- getGEO(GEO, destdir = dir, getGPL = F)
    e <- exprs(eSet[[1]])
    pdata <- pData(eSet[[1]])
    write.csv(e, paste0(GEO, "_expr.csv"))
    write.csv(pdata, paste0(GEO, "_pdata.csv"))
    save(eSet, file = paste0(GEO, ".eSet.Rdata"))
    return(eSet)
  }
  `prepGSE(x)`('GSE62832')
}

#当有多个GSE要处理时,先造一个GEO_lst, 再用函数prepGEO
if(F){
  for (i in 1:length(GEO_lst)){
  id <- GEO_lst[i]
  `prepGSE(x)`('id')
}
}


### 写成函数
GSE_expr <- function(GSE){
  library(GEOquery)
  if(!file.exists(GSE)){
    geo <<- getGEO(GSE, destdir = '.', getGPL = F, AnnotGPL = F)
    gdata <- paste0(GSE,'.eSet.Rdata')
    save(geo, file = gdata)
  }
  load(gdata)
  expr <<- exprs(geo[[1]])
}
# e.g.
GSE_expr("GSE17215")

