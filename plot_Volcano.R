##### import R packages
rm(list=ls())
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))


##### read data
dataset <- fread("hs_HCMvsNCM.csv")
dataset <- fread("hs_HFvsNCM.csv")
dataset <- fread("mm_TAC_HTvsNHCM.csv")
dataset <- fread("mm_TAC_HTvsNHT.csv")

##### set pvalue and logFC
cut_off_pvalue = 0.05
cut_off_logFC = 1

##### save‘up’, ‘Down’ and ‘Stable’ to change column;change column for color setting
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC, 
                        ifelse(dataset$log2FoldChange> cut_off_logFC ,'Up','Down'),
                        'Stable')
##### plot
p <- ggplot(
  dataset, 
  aes(x = log2FoldChange, 
      y = -log10(pvalue), 
      colour=change)) +
  geom_point(alpha=0.8, size=1) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  
  # line
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  
  # axis
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  
  # theme
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )

##### mark genes
dataset$label = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= 5, as.character(dataset$gene),"")
morbid_index <- which(dataset$gene == "Morrbid")
if (length(morbid_index) > 0) {
  dataset$label[morbid_index] <- "Morrbid"
}

p + geom_text_repel(data = dataset, aes(x = dataset$log2FoldChange, 
                                        y = -log10(dataset$pvalue), 
                                        label = label),
                    size = 3,box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"), 
                    segment.color = "black", 
                    max.overlaps = 100, 
                    show.legend = FALSE)

