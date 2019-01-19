### 读取Affymetrix芯片的CEL-- 利用Affy包
#【支持平台有限：一般是hgu 95系列和133系列】
# 目的：利用CEL生成表达矩阵
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1428/suppl/
dir_cels <- "/Users/reedliu1/Downloads/tmp/GSE1428_RAW"
library(affy) #可以用于hgu 95系列和133系列

# 标准化方法一：mas5
affy_data = ReadAffy(celfile.path=dir_cels)
eset.mas5 = mas5(affy_data)
exprSet.nologs = exprs(eset.mas5)
exprSet = log(exprSet.nologs, 2)  #transform to Log_2 if needed

# 标准化二：rma
library(affy)

# 读取
rma_data<- ReadAffy(celfile.path=dir_cels) 
eset.rma <- rma(rma_data)
# 统计在所有样本中表达都为0的基因
calls <- mas5calls(rma_data) # get PMA calls
calls <- exprs(calls)
absent <- rowSums(calls == 'A') # how may samples are each gene 'absent' in all samples
absent <- which (absent == ncol(calls)) # which genes are 'absent' in all samples

# 过滤
rmaFiltered <- eset.rma[-absent,] # filters out the genes 'absent' in all samplesc

# 输出
write.exprs(rmaFiltered,file="data.txt")

###################################################################
### 读取Affymetrix 其他芯片，如 Human Gene 1.1 ST Array, 用oligo包
library(oligo)
celFiles <- list.celfiles("~/Downloads/tmp/GSE48452_RAW/", listGzipped = T, full.names=T)
# 读取cel文件 
affy_data <- read.celfiles(celFiles) 
# 标准化【比较慢】
geneCore <- rma(affy_data, target = "core")
genePs <- rma(affy_data, target = "probeset")
# 探针ID赋值给表达矩阵
featureData(genePs) <- getNetAffx(genePs, "probeset")
featureData(geneCore) <- getNetAffx(geneCore, "transcript")

# 探针ID转基因ID
# 获取GSE对应平台GPL信息
# 然后下载、加载对应注释包
if(T){
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("hugene11sttranscriptcluster.db", version = "3.8")
    library( hugene11sttranscriptcluster.db) 
    s <- toTable(hugene11sttranscriptclusterSYMBOL)
    
    #过滤【只保留和注释文件探针id相同的探针】
    e <- exprs(geneCore)
    efilt <- e[rownames(e)%in%s$probe_id,]
    maxp = by(efilt,s$symbol,function(x) rownames(x)[which.max(rowMeans(x))]) 
    uniprobes = as.character(maxp)
    efilt=efilt[rownames(efilt)%in%uniprobes,]
    rownames(efilt)=s[match(rownames(efilt),s$probe_id),2]
}

###################################################################
# 读取更新的芯片数据， 如HST2.0（多方面优于转录组）、Illumina HumanHT-12 V4.0 expression beadchip
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("lumi", version = "3.8")
library(lumi)
celFile <- '~/Downloads/tmp/GSE30669_HEK_Sample_Probe_Profile.txt.gz'
lumi_data <- lumiR.batch(celFile)
pData(phenoData(lumi_data))
lumi.N.Q <- lumiExpresso(lumi_data)
dataMatrix <- exprs(lumi.N.Q)
# 以下代码相同效果
rm(list=ls())
library(GEOquery)
library(limma)
GSE30669 <- getGEO('GSE30669', destdir=".",getGPL = F)
exprSet=exprs(GSE30669[[1]])
pdata=pData(GSE30669[[1]])
exprSet=exprs(GSE30669[[1]])


