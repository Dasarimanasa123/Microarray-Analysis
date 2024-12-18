
library(ggplot2)
library(ggrepel)
library(dplyr)

de_genes <- read.csv("deseq_results_filtered.csv",sep ="," )

head(de_genes)
colnames(de_genes)

de_genes$deffexpressed <- "NO"

dim(de_genes)
de_genes$deffexpressed[de_genes$log2FoldChange  > 0.58 & de_genes$padj < 0.01] <- "up"

de_genes$deffexpressed[de_genes$log2FoldChange  < 0.58 & de_genes$padj < 0.01] <- "down"

head(de_genes)

de_genes$delable <- NA

ggplot(data = de_genes,aes(x =log2FoldChange, y = -log10(pvalue),col=deffexpressed ,label = delable))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values = c('blue','black','red'))+
  theme(text=element_text(size=30))



ggplot(data = de_genes,aes(x =log2FoldChange , y = -log10(pvalue),col=deffexpressed ,label = delable))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values = c('blue','black','red'))+
  geom_vline(xintercept = c(-0.8,0.8),col="red")+
  geom_hline(yintercept = -log10(0.001),col="red")+
  theme(text=element_text(size=20))


head(de_genes)

arrange(de_genes,pvalue)
head(arrange(de_genes,pvalue),10)

head(arrange(de_genes,pvalue),10)$pvalue
thres = head(arrange(de_genes,pvalue),10)$pvalue[10]
head(arrange(de_genes,pvalue),10)$X

de_genes$delable[de_genes$pvalue <= thres] <- de_genes$X[de_genes$pvalue <= thres]

head(arrange(de_genes,pvalue),10)

#ggplot(data = de_genes,aes(x =log2FoldChange , y = -log10(pvalue),col=deffexpressed ,label = delable))+
ggplot(data = de_genes,aes(x =log2FoldChange, y = -log10(pvalue),col=deffexpressed ,label = delable))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values = c('blue','black','red'))+
  theme(text=element_text(size=30))

write.csv(de_genes , "de_genes.csv")

