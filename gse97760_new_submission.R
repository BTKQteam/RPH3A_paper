library(affy)
library(limma)
library(sva)
library(GEOquery)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(dplyr)
library(ggrepel)
library(IOBR)
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
library(survival)
library(dplyr)
library(tidyr)
# load series and platform data from GEO
matrix_dir<-"G:/My Drive/R code for all project/Quy/matrixdata"

setwd("G:/My Drive/R code for all project/")

gset <- getGEO("GSE97760", GSEMatrix =TRUE, getGPL=T,destdir = matrix_dir)
if (length(gset) > 1) idx <- grep("GPL16699", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]


## calculate median expression level
cutoff <- median(exprs(gset))

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(gset) > cutoff

## Identify genes expressed in more than 2 samples

keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
gset <- gset[keep,]
# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

##clear data matrix
matrix<-as.data.frame(exprs(gset))

#create matrix with gene GENE_SYMBOL and expression value
features_new <- fData(gset)
features_new <- select(features_new,GENE_SYMBOL)
matrix <- cbind(features_new, matrix) %>%
  filter(GENE_SYMBOL != "" & !duplicated(gsub(" ///.*", "", GENE_SYMBOL))) %>%
  mutate(GENE_SYMBOL = gsub(" ///.*", "", GENE_SYMBOL)) %>%
  as.data.frame()

rownames(matrix) <- matrix$GENE_SYMBOL
matrix$GENE_SYMBOL <- NULL


# log2 transformation
ex <- matrix
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
matrix <- log2(ex) }


##Extract sample Info
sampleInfo <- pData(gset)

##Extract sample Info
sampleInfo <- pData(gset)
sampleInfo<-select(sampleInfo,geo_accession,characteristics_ch1.2)
##rename sampleInfo
sampleInfo<-sampleInfo%>%dplyr::rename(ID=geo_accession,phenotype=characteristics_ch1.2)
sampleInfo<-sampleInfo%>%mutate(phenotype=gsub("disease: healthy","control",phenotype,fixed = T))
sampleInfo<-sampleInfo%>%mutate(phenotype=gsub("disease: advanced Alzheimer's disease","AD",phenotype,fixed = T))

#BEFORE EXCLUDE OUTLIER, RUN DEG FIRST:

# Convert the phenotype data to a design matrix
design <- model.matrix(~0+ phenotype,data=sampleInfo)
# Create a linear model and fit it to the gene expression data
fit <- lmFit(matrix, design)
# Create contrasts for the groups of interest
contrast.matrix <- makeContrasts(AD_vs_control = phenotypeAD - phenotypecontrol, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)

# Compute moderated t-statistics for each contrast
fit2<- eBayes(fit2,0.01)
full_results <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
full_results$Gene<-rownames(full_results)
full_results$group <- "non-significant"
full_results$group[full_results$logFC > 0.5 & full_results$P.Value < 0.05] <- "Up (eoCRC)"
full_results$group[full_results$logFC < -0.5 & full_results$P.Value < 0.05] <- "Down (control)"


logFC <- full_results$logFC
Pval <- full_results$P.Value
adj <-full_results$adj.P.Val
ggplot(data = data.frame(full_results), aes(x = logFC, y = -log10(Pval),col = group)) +
  geom_point(size = 0.3, alpha = 0.8)+
  guides(color = guide_legend(override.aes = list(size = 3)))+ theme_bw()+
  scale_color_manual(values = c("darkseagreen", "grey80", "lightcoral"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + 
  labs(x = "Log2 Fold Change (AD vs Control)", y = "-log10(P-value)",
       title = "GSE97760") +theme_bw()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed")+ 
  theme(panel.border = element_blank(),
        axis.text = element_text(size = 15),  # Adjust the size of axis labels
        axis.title = element_text(size = 20),  # Adjust the size of axis titles
        plot.title = element_text(size = 18,face = 'bold',hjust = 0.5, vjust = 1))

full_results_new<-full_results[full_results$adj.P.Val<0.05 & abs(full_results$logFC)>0.5,]
full_results_new$group<-ifelse(full_results_new$logFC>0.5,"up","down")


topN <- 40
##
ids_of_interest <- mutate(full_results_new, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(Gene)
## Get the rows corresponding to ids_of_interest and all columns
gene_matrix <- as.matrix(matrix[ids_of_interest,])
annot<-data.frame(group=sampleInfo$phenotype,row.names = rownames(sampleInfo))
ComplexHeatmap::pheatmap(mat = gene_matrix,
         labels_row = rownames(gene_matrix),
         scale="row",annotation_col = annot,
         annotation_colors = list(group=c(control='darkseagreen',AD='lightcoral')),
         show_colnames = FALSE)


## PCA analysis
#MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(matrix))

## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=phenotype,label=paste("", ID))) + geom_point(size=0.1) +
  geom_text_repel(aes(label = ID), size = 2.5, nudge_x = 1, nudge_y = 1)+
  theme_bw()+stat_ellipse()+geom_hline(yintercept = 0, linetype = "dashed",lwd=0.3)+
  geom_vline(xintercept = 0, linetype = "dashed",lwd=0.3)+ theme(panel.border = element_blank())+
  scale_color_manual(values = c("#F39C12","#3498DB"))
#writexl::write_xlsx(full_results_new,path = '97760_adj.xlsx')


#### Check expression of stem cell pathway
stem_list<-c('WNT4','RIF1','IL6ST','BMPR1A','PIK3CA',
             'REST','SKIL','AKT3','KRAS','SMARCAD1',
             'BMPR2','ESRRB','BMI1','FZD6','APC',
             'NRAS','AXIN2','MEIS1','JAK2','IGF1',
             'SMAD2','FGFR4','WNT9B','PIK3R1')

# Filter the expression matrix for these genes
expression_data <- matrix[stem_list, ]

# Ensure 'expression_data' is in a suitable format for merging
expression_data <- as.data.frame(t(expression_data))
expression_data$ID <- rownames(expression_data)
merged_data <- merge(expression_data, sampleInfo, by = "ID")

# Melt the data for ggplot2
library(reshape2)
melted_data <- melt(merged_data, id.vars = c("ID", "phenotype"), variable.name = "gene", value.name = "expression")
# Create a new variable to distinguish between gene and group combinations
melted_data$gene_group <- interaction(melted_data$gene, melted_data$phenotype)

comparisons <- list(c("AD", "control"))

p1 <- ggplot(melted_data, aes(x = phenotype, y = expression, fill = phenotype)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("AD", "control")), method = "t.test",
                                       label = "p.format")+
  facet_wrap(~gene,ncol = 6) +ylim(0,18)+
  labs(title = "Signaling pathway regulating pluripotency of stem cells",
       x = "Group",
       y = "Expression Level") + theme_pubr()+
 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        strip.text = element_text(size = 10,face = 'bold'), # Enlarges the gene names in the facet labels
        legend.title = element_text(size = 12), # Enlarges the legend title
        legend.text = element_text(size = 15)) +
  scale_fill_manual(values = c("AD" = "lightcoral", "control" = "darkseagreen"))
# Add statistical comparison

p1
### combine in 1 coordinate
# p2  <- ggboxplot(melted_data, x = "gene", y = "expression", color = "phenotype", 
#                      add = "jitter", palette = "jco", 
#                      title = "Gene Expression Comparison between AD and Control",
#                      xlab = "Gene", ylab = "Expression Level", 
#                      legend.title = "Group") +
#   stat_kruskal_test(aes(group = phenotype), method = "t.test")
# 
# p2
############
############
## GSEA   ##
############
library(msigdbr)
library(fgsea)
kegg<-msigdbr(species = "Homo sapiens", category = "C5",subcategory = 'GO:BP')
kegg_list<-split(x = kegg$gene_symbol,f=kegg$gs_name)
a<-msigdbr_collections()
## define input
deg<-full_results[order(full_results$logFC,decreasing = T),]
genelist <- deg$logFC
names(genelist) <- rownames(deg)

gsea.re2 <- fgsea(pathways = kegg_list,
                  stats = genelist,
                  nperm=1000,
                  minSize=1,
                  maxSize=10000)

colnames(gsea.re2)
g2 <- gsea.re2[gsea.re2$pval<0.05&abs(gsea.re2$NES)>1.5,]
g2 <- g2[order(g2$NES,decreasing = T),]
kegg_new<-kegg[kegg$gs_name%in%g2$pathway,]
###get back to list
geneSetList <- kegg_new %>% split(x = .$gene_symbol, f = .$gs_name)
our_target <- c("SV2A", "KCNQ2", "RPH3A")

# Filter the geneSetList to include only pathways with at least one of the genes of interest
filtered_geneSetList <- geneSetList[sapply(geneSetList, function(genes) {
  any(genes %in% our_target)
})]

### Make the input for lolipop analysis
lolipop_input<- kegg_new[kegg_new$gs_name%in%names(filtered_geneSetList)&kegg_new$gene_symbol%in%our_target,]
lolipop_input<-select(lolipop_input,gs_name,gene_symbol)
#### Add NES score into lolipop input
lolipop_input <- lolipop_input %>%
  left_join(g2, by = c("gs_name" = "pathway"))%>%select(gs_name,gene_symbol,NES)



#writexl::write_xlsx(g2,path = 'G:/My Drive/R code for all project/Quy/excel_output_new/gsea_out/97760_GO.xlsx')
#writexl::write_xlsx(kegg_new,path = 'G:/My Drive/R code for all project/Quy/excel_output_new/gsea_out/97760_GO_genelist.xlsx')
#writexl::write_xlsx(lolipop_input,path = 'G:/My Drive/R code for all project/Quy/excel_output_new/gsea_out/97760_lolipop.xlsx')
###taking leading edge gene list



####visualization (gggsea)
library(gggsea)
library(ggsci)
names(kegg_list)
se_hall<-c(head(g2$pathway,4),tail(g2$pathway,4))


sig1<-gsea.re2[gsea.re2$pathway%in%se_hall,]

kegg.se<-kegg_list[sig1$pathway]
df.new <- gseaCurve(genelist , kegg.se, sig1)
df.new$set<-factor(df.new$set,levels = sig1$pathway)


pal_line<-pal_lancet()(9)

ggplot() + 
  geom_gsea(df.new,
            prettyGSEA=T,
            tickcolor='grey30',
            ticksize=0.2,
            
            linecolor='#5DADE2',
            linesize=1,lty=1,
            
            nrow = 4,
            ncol = 2
  ) + 
  theme_bw()+
  theme(strip.text = element_text(size = 8,face = 'bold'),
        strip.background = element_rect(fill = 'white'))+
  xlab(bquote(italic('Rank')))+ylab(bquote(italic('Enrichment Score')))+
  theme(axis.text.x = element_text(size=12,angle = 0,face = 'plain',hjust = 0.5),
        axis.text.y = element_text(size=12,angle = 0,face = 'plain',vjust = 0.5))

