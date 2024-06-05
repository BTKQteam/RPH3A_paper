library(TCGAWorkflowData)
# Load packages
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library("DESeq2")
library(dplyr)
library(pheatmap)
library(colorRamp2)
library(ggrepel)
library(IOBR)
library(GEOquery)
library(dplyr)
library(scrabble)
library(edgebundleR)

##### input data
lolipop<-read.csv('G:/Shared drives/BTKQ/Quy/Quy-AD-2024/Data/gsea_out/Lolipop/lolipop_SFG_new.csv',check.names = F)
lolipop$Gene<-factor(lolipop$Gene,levels = c('KCNQ2','RPH3A','SV2A'))
df <- lolipop %>%
  dplyr::filter(Gene %in% c( "KCNQ2","RPH3A","SV2A")) %>%
  dplyr::group_by(pathway, Gene) %>%
  dplyr::summarise(counts = n())

###### count
ggdotchart(lolipop, x = 'Gene', y = "pathway",
           color = "Gene",                                # Color by groups
           palette = c("#00AFBB", "#E7B800", "#FC4E07","#44AD56"), # Custom color palette
           #sorting = "descending",                       # Sort value in descending order
           #add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           group = "Gene",                                # Order by groups
           dot.size = 7,     
           rotate = F,# Large dot size
          # label = round(lolipop$NES,1),                        # Add mpg values as dot labels
          # font.label = list(color = "white", size = 12, vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                      # ggplot2 theme
)+ theme_cleveland()   +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  theme(
    legend.text = element_text(size = 15),      # Adjust the size of legend text
    axis.text = element_text(size = 20),        # Adjust the size of axis labels
    axis.title.y = element_blank(), 
    axis.title.x = element_text(size = 15,face = 'bold'),# Hide x-axis tick labels
    # Adjust the size of axis titles
    strip.text = element_text(size = 15),       # Adjust the size of strip (facet) labels
    plot.caption = element_text(size = 15),     # Adjust the size of the plot caption
    axis.text.x = element_text(size = 15),      # Adjust the size of x-axis tick labels
    axis.text.y = element_text(size = 12),      # Adjust the size of y-axis tick labels
    strip.text.x = element_text(size = 15)      # Adjust the size of x-axis facet labels
  )

