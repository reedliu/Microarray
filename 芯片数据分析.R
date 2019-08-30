#################################
####### 芯片数据分析 ###########
#################################

#################################
# 1.下载数据
#################################
#方法一：自己下载【自己定义函数】#自己再修改下
'get_GSE_links(x)' <- 
  function(studyID = x, down = F, destdir = "./") {
    ## studyID destdir ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62832/matrix/
    matrix_link = paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/", substr(studyID, 1, nchar(studyID) - 3), "nnn/", studyID, "/matrix/", studyID, "_series_matrix.txt.gz")
    print(paste0("The URL for matrix   is : ", matrix_link))
  }
get_GSE_links('GSE62832')

#方法二：用GEOquery下载
if(T){
source("http://bioconductor.org/biocLite.R")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
biocLite("GEOquery")

library(GEOquery)
eSet <- getGEO('GSE62832', destdir = '.', getGPL = F, AnnotGPL = F)
save(eSet, file = 'GSE62832.eSet.Rdata')
}
#当有多个GSE要处理时
if(F){for (i in 1:length(GEO_lst)){
  id <- GEO_lst[i]
  e <- getGEO(id, destdir='.')
}
}

#################################
# 2.将下载的表达矩阵读进来
#################################
#方法一：针对自己下载好的矩阵文件
e1 <- read.csv('GSE62832_series_matrix.txt.gz',
              comment.char = '!', fill = T, sep = '\t')
#再将e1的行名转为ID号，也就是把目前的第一列复制行名的位置，然后再把第一列去掉
rownames(e1) <- e1[,1]
e <- e1[,-1]

#方法二：针对GEOquery下载的eSet
load('GSE62832.eSet.Rdata')
class(eSet) #先查看数据种类【列表还是数据框】
str(eSet) #再看数据结构【可以看到第二行涉及到了表达矩阵】
#既然是列表，那么从列表中提取信息就是 eSet[[1]]
e <- exprs(eSet[[1]])

#################################
# 2附加.将表达矩阵、表型分组一键输出
#################################
'downGSE(x)' <- function(GEO = x, dir = ".") {
  library(GEOquery)
  e <- getGEO(GEO, destdir = dir, getGPL = F)
  eSet = exprs(e[[1]])
  pdata = pData(e[[1]])
  write.csv(eSet, paste0(GEO, "_eSet.csv"))
  write.csv(pdata, paste0(GEO, "_metadata.csv"))
  save(eSet, file = paste0(GEO, ".eSet.Rdata"))
  return(e)
}
`downGSE(x)`('GSE62832')

####################################
# 3.将表达矩阵的探针ID转换为gene ID
####################################
#3.1首先需要知道GSE62832对应平台是GPL6244 
eSet

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
#整合1【目的：保证一个基因对应一个探针；如果基因和探针一一对应很好说，但如果一个基因对应多个探针：每个探针取一行的均值-》对应同一基因的探针取表达量最大的探针-》按照基因名给他们建索引，因为是按照基因来过滤探针（不用s$probe_id构建索引的原因是，看清楚我们的目的是让注释包的一个基因对应我们自己表达矩阵的一个探针。如果用s$probe_id那么结果就成了让注释包的一个探针对应我们自己表达矩阵的一个探针，当然这样也运行不成功，因为自己表达矩阵的探针过滤后的数量和注释包的探针数量不相等，这样没法一一对应。但基因名数量是不变的，什么是索引？以不变应万变的就是索引）】
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

####################################
# 5.对表达矩阵进行差异分析
####################################
#5.1 只需要提供表达矩阵efilt、分组信息group_list，就能使用limma进行分析
suppressMessages(library(limma))
#limma需要三个矩阵：表达矩阵（efilt）、分组矩阵（design）、比较矩阵（contrast）
#先做一个分组矩阵～design，说明MAO是哪几个样本，MNO又是哪几个，其中1代表“是”
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(efilt)
design
#再做一个比较矩阵【一般是case比control】
contrast<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast

#5.2 准备就绪，就可以开始差异分析
DEG <- function(efilt,design,contrast){
  ##step1
  fit <- lmFit(efilt,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast) 
  fit2 <- eBayes(fit2)  
  ##step3
  mtx = topTable(fit2, coef=1, n=Inf)
  deg_mtx = na.omit(mtx) 
  return(deg_mtx)
}
DEG_mtx <- DEG(efilt,design,contrast) #得到全部的差异基因矩阵
save(DEG_mtx, efilt, file = "GSE62832.DEG.Rdata")
#5.3 对一个小表达矩阵，如30、50个基因，可以用热图
top30_gene=head(rownames(DEG_mtx),30)
top30_matrix=efilt[top30_gene,] #得到top30的表达量矩阵
top30_matrix=t(scale(t(top30_matrix)))
#这里做个top30 gene heatmap
if(T){
  library(ComplexHeatmap)
  #install.packages("dendextend")
  library(dendextend)
  library(dendsort)
  e_mean = tail(sort(apply(efilt,1,mean)),30)
  top30_cluster_cols = hclust(dist(t(top30_matrix)))
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...))) #聚类的过滤函数
  top30_cluster_cols = sort_hclust(hclust(dist(t(top30_matrix))))
  top30_cluster_rows <- sort_hclust(hclust(dist(top30_matrix)))
  base_mean = e_mean
  ef_lable = HeatmapAnnotation(text = anno_text(colnames(efilt), rot = 45,offset = unit(23, "mm"), gp = gpar(fontsize	= 8.5)))
  Heatmap(top30_matrix, name = "expresssion", 
          column_title = "Samples", 
          column_title_gp = gpar(fontsize = 14, fontface = "bold"),
          row_title = "Genes",
          row_title_gp = gpar(fontsize = 14, fontface = "bold"),
          row_names_gp = gpar(fontsize = 8.5),
          row_names_side = "left",
          cluster_rows = color_branches(top30_cluster_rows, k = 4),
          cluster_columns = color_branches(top30_cluster_cols, k = 3),
          show_column_names = F,
          bottom_annotation = ef_lable, bottom_annotation_height = unit(0.5, "cm"))+
    Heatmap(base_mean, name = "base_mean", show_row_names = FALSE, width = unit(5, "mm"))
}

#5.4 对一个大的表达矩阵，如全部的差异基因，可以用火山图
#火山图实际上就是根据两列进行作图：logFC、pvalue
plot(DEG_mtx$logFC, -log10(DEG_mtx$P.Value)) #最简单的图（能说明原理）
#一个好看点的图
DEG=DEG_mtx
if(T){logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) ) 
DEG$result = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3), #round保留小数位数
                    '\nThe number of up gene is ',nrow(DEG[DEG$result =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$result =='DOWN',])
)
library(ggplot2)
ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=result)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red')) ## 这里要注意和之前设置的result三个因子相对应，DOWN就设为blue，NOT就设为black
}

#使用DESeq2+EnhancedVolcano【待完善！】
if(F){
library("DESeq2")
dds <- DESeqDataSet(airway, design = ~cell + dex)
dds <- DESeq(dds, betaPrior = FALSE)
res2 <- results(dds, contrast = c("cell", "N061011", "N61311"))
res2 <- lfcShrink(dds, contrast = c("cell", "N061011", "N61311"), res = res2)
install.packages("devtools")
library(devtools)
devtools::install_github("kevinblighe/EnhancedVolcano")
library(EnhancedVolcano)
EnhancedVolcano(DEG,
                lab = rownames(DEG),
                x = "log2FoldChange", y = "padj",
                selectLab = top30_gene,
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = 0.0001,
                FCcutoff = 2.0,
                xlim = c(-6,6),
                transcriptPointSize = 1.8,
                transcriptLabSize = 5.0,
                colAlpha = 1,
                legend=c("NS","Log2 FC","Adjusted p-value",
                         "Adjusted p-value & Log2 FC"),
                legendPosition = "bottom",
                legendLabSize = 10,
                legendIconSize = 3.0)
}                
                
####################################
# 6.对结果进行注释～富集分析（基于超几何分布检验）
####################################      
rm(list=ls())
load("GSE62832.DEG.Rdata")
#6.1 得到上调、下调基因
DEG=DEG_mtx
if(T){
  logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) ) 
  DEG$result = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
}
gene <- rownames(DEG[DEG$result != 'NOT', ])
#要再具体细分的话也可以
up_gene <- rownames(DEG[DEG$result =='UP',])
down_gene <- rownames(DEG[DEG$result =='DOWN',])

#6.2 ID转换【这个包接受的gene是Entrez ID，因此要先转换SYMBOL=》Entrez】
#几种主流的ID
#########
#Ensemble id：由欧洲生物信息数据库提供，一般以ENSG开头，后边跟11位数字。如TP53基因：ENSG00000141510
#Entrez id：由美国NCBI提供，通常为纯数字。如TP53基因：7157
#Symbol id：为我们常在文献中报道的基因名称。如TP53基因的symbol id为TP53
#Refseq id：NCBI提供的参考序列数据库：可以是NG、NM、NP开头，代表基因，转录本和蛋白质。如TP53基因的某个转录本信息可为NM_000546
#########
library(clusterProfiler) 
library(org.Hs.eg.db)
if(T){
  gene_tr <- bitr(gene, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  up_gene_tr <- bitr(up_gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
  down_gene_tr <- bitr(down_gene, fromType = "SYMBOL",
                   toType = c("ENSEMBL", "ENTREZID"),
                   OrgDb = org.Hs.eg.db)
  gene <- gene_tr$ENTREZID
  up_gene <- up_gene_tr$ENTREZID
  down_gene <- down_gene_tr$ENTREZID
} #结果发现：总体有3.07%没有对应，上调有2.13%没有对应，下调3.76%没有对应
head(gene_tr)
#转换完可以去检验一下，比如检查Entrez ID: 23336, https://www.ncbi.nlm.nih.gov/gene/23336。确实得到的是SYNM这个基因

##############################################################################################
#进行注释GO terms/Pathway（KEGG、BIOCARTA、Reactome）/MSigDB【最常见的就是GO/KEGG数据库注释】
#每个富集分析数据库对应一个R包

#6.3 KEGG
#关于KEGG的解释
#######
#使用超几何分布检验：
#比如背景基因这里有18837个，选出来的差异基因（上调235个+下调319个）554个，约占总体3%；
#目前KEGG 包含了530个pathway，其中一个04110是Cell Cycle通路，其中包含124个基因【怎么统计？打开链接-》https://www.kegg.jp/dbget-bin/www_bget?hsa04110，找到gene栏-》linux cat+wc】
#如果是随机抽样的话，cell cycle在我们这554个基因中，应该能抽到124*3%=3个左右。也就是说，在得到的554个差异基因中，应该只有3个是cell cycle的基因；但结果找到了30个cell cycle的基因，那么就说明差异基因的产生都有cell cycle这个因素。如果是药物处理组的结果发现这样，就说明确实可能是药物起了作用，当然后续还需要结合其他分析
#并非所有的通路都需要检测：有的通路比较小，只有十几个基因，我们这500多个差异基因中能抽到的理论上最多一个，这样一般就不考虑；还有的通路很大，有上千个基因（有个GO pathway有5k+基因），我们的背景基因才不到2w，那我们这586个差异基因中在这个pathway中的理论上就得有一百多个，这样和实际的变化差异就不明显
######
#6.3.1 enrichKEGG进行候选基因通路分析
if(T){
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
kk_gene <- (kk)[,1:6]
}
dotplot(kk,showCategory = 14,color="pvalue",font.size=14)  
#按照包的说明做一个geneList（标准：Entrez ID对应logFC，然后logFC降序排列）
if(T){
geneList <- DEG_mtx$logFC #先获得logFC【没有重复的18837个数值】
names(geneList) <- rownames(DEG_mtx) #获得基因名，但这里是SYMBOL格式【也是没有重复的18837个名称】
##问题就出在（下一步）SYMBOL转Entrez名字上，会多出来许多的重复：如果只是SYMBOL转Entrez，那么会得到18841个名称，原因就是一个symbol对应多个entrez ID，用count(geneList_tr,geneList_tr$SYMBOL) %>% arrange(desc(n))这样来查看到底哪些基因，结果查到有4个基因对应了2个entrez ID，其中一个是“HBD”，然后使用filter(geneList_tr, geneList_tr$SYMBOL == "HBD")来看HBD对应了哪个ID，果然发现HBD对应了两个ID；如果将SYMBOL转ENSEMBL+Entrez，那么会得到20198个名称，因为同一个SYNBOL，同一个Entrez，会有不同的Ensembl ID
#另外还有别人的友情提示：ENTREZID(具有唯一性）：基因编号来自ANNOVAR的注释结果，建议别用SYMBOL，因为这种名称特异性较差，在转成ENTREZID时可能出现不唯一的现象。symbol与entrezID并不是绝对的一一对应的
geneList_tr <- bitr(names(geneList), 
                    fromType = "SYMBOL",
                    toType = c("ENSEMBL","ENTREZID"),
                    OrgDb = org.Hs.eg.db) 
#下一步目的：将Entrez ID与logFC联系起来
#因为转换后的名称已经不是和原来的geneList行名一一对应了，所以不能直接将原来SYMBOL行名变为Entrez ID。但是可以用SYMBOL名作为桥梁，一边连接logFC（做一个数据框），另一边连接Entrez ID（已经有的geneList_tr数据框），再将两个数据库按照SYMBOL联系起来
new_list <- data.frame(SYMBOL=names(geneList), logFC = as.numeric(geneList)) 
new_list <- merge(new_list, geneList_tr, by  = "SYMBOL")

geneList <- new_list$logFC #重新生成geneList之logFC
names(geneList) <- geneList_tr$ENTREZID #重新生成geneList之Entrez ID

geneList <- sort(geneList,decreasing = T) #降序排列
}
#6.3.2 gseKEGG GSEA分析
if(T){
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.2,
               verbose      = FALSE)
}
head(kk2)[,1:6]
gseaplot(kk2, geneSetID = "hsa04310")
#6.3.3 查看特定通路图
library(pathview)
hsa01212 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa01212", #上述结果中的hsa01212通路
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))

#6.4 GO
if(T){
ego_CC <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
}
ego_gene <- head(ego_CC)[,1:6]
##可视化--点图
dotplot(ego_CC,title="EnrichmentGO_CC_dot")#点图，按富集的数从大到小的
##可视化--条形图
barplot(ego_CC, showCategory=20,title="EnrichmentGO_CC")#条状图，按p从小到大排，绘制前20个Term
#通路图
library(Rgraphviz)
library(topGO)
plotGOgraph(ego_CC)



















