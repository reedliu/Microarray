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
  #top30_cluster_cols = hclust(dist(t(top30_matrix)))
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...))) #聚类的过滤函数
  top30_cluster_cols = sort_hclust(hclust(dist(t(top30_matrix))))
  top30_cluster_rows <- sort_hclust(hclust(dist(top30_matrix)))
  base_mean = e_mean
  ef_lable = HeatmapAnnotation(text = anno_text(colnames(efilt), rot = 45,offset = unit(23, "mm"), gp = gpar(fontsize	= 8.5)))
  p <- Heatmap(top30_matrix, name = "expresssion", 
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
  print(p)
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
p <- ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=result)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red')) ## 这里要注意和之前设置的result三个因子相对应，DOWN就设为blue，NOT就设为black
print(p)
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
