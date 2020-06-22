rm(list = ls())
options(stringsAsFactors = F)

library(stringr)
library(ggplot2)

## 读取数据
exp =rio::import("test.xlsx")
exp = exp[,!stringr::str_detect(colnames(exp),"Average")]
exp = exp[,!stringr::str_detect(colnames(exp),"Ratio")]
rownames(exp) = exp$TargetID
exp = exp[,-1]
exp = log2(exp+1)
exp = exp[apply(exp, 1, sum) > 0 ,]

## limma需要三个矩阵：表达矩阵（exp）、分组矩阵（design）、比较矩阵（contrast）
suppressMessages(library(limma))
# 设置分组信息(一个control两个treat组，需要分别比较 treat1-control 和 treat2-control)
group_list = rep(c("control","treat1","treat2"),each = 2)
group_list=factor(group_list,levels = c("control","treat1","treat2"))
#先做一个分组矩阵～design，其中1代表“是”
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- as.character(group_list)
design
#再做一个多组比较矩阵【一般是case比control】
unique(group_list)
contrast <- makeContrasts(contrasts=c("treat1-control", "treat2-control"),levels = design) 
contrast


##开始差异分析
DEG <- function(efilt,design,contrast,coef=1){
  ##step1
  fit <- lmFit(efilt,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast) 
  fit2 <- eBayes(fit2)  
  ##step3
  mtx = topTable(fit2, coef=coef, n=Inf)
  deg_mtx = na.omit(mtx) 
  return(deg_mtx)
}

trt1_mtx <- DEG(exp,design,contrast,coef = 1) 
trt2_mtx <- DEG(exp,design,contrast,coef = 2) 

## 画图
data <- trt1_mtx
library(ggrepel)
library(ggplot2)
data$gene <- rownames(data)
data$change = as.factor(ifelse(data$P.Value < 0.05 & abs(data$logFC) > logFC_cutoff,
                               ifelse(data$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
logFC_cutoff=0.6
# 图片title
title <- paste0('CBX6.6 log2FoldChange cutoff: ',round(logFC_cutoff,3),
                '\nUp-regulated genes: ',nrow(data[data$change =='UP',]) ,
                '\nDown-regulated genes: ',nrow(data[data$change =='DOWN',]))

vol2 = ggplot(data=data, aes(x=logFC, y =-log10(P.Value),color =change)) +
  #设置点的透明度、大小等
  geom_point(alpha=0.4, size=1.75) +
  #设置点的颜色
  scale_color_manual(values =c("blue","black","red"))+
  #设置横竖阈值线（横：p为0.05，竖：FC的阈值）
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=1)+
  geom_vline(xintercept = c(logFC_cutoff,-logFC_cutoff),lty=4,lwd=0.6,alpha=1)+
  #更换背景主题
  theme_bw()+
  #更换背景格子
  # theme(panel.border = element_blank(),
  #       panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank(),   
  #       axis.line = element_line(colour = "black"))+
  #加标题(见上面title)
  # title：
  labs(title=title, x="log2 (Fold Change)",y="-log10 (P-value)")+
  #标题居中
  theme(plot.title = element_text(hjust = 0.5))
  # 设置x、y轴坐标范围
  # xlim(-10,10)

write.csv(data,file = "treat1-control.csv")
ggsave(vol2,filename = "treat1-control.png")




