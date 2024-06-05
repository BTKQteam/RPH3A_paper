library(Matrix)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(SCENIC)
library(celldex)
library(scRNAseq)
library(ggpubr)
library(hdf5r)
library(RColorBrewer)
library(devtools)
library(AnnoProbe)
#library(rjags)
#library(infercnv)
library(presto)
library(SCP)
library(reticulate)
library(forcats)
library(copykat)
library(harmony)
library(ggrepel)

#get dir
h5_files <- list.files(pattern = "*.h5")
h5_read <- lapply(h5_files, Read10X_h5)
h5_seurat <- lapply(h5_read, CreateSeuratObject)
h5_seurat  <- merge(h5_seurat[[1]], y = h5_seurat[2:length(h5_seurat)], 
                    add.cell.ids = c("EC2", "EC1", "EC3","EC4","EC6","EC7","EC5","EC9","EC8","EC10"), project = "project")
view(h5_seurat@meta.data)
# create a sample column
h5_seurat$sample <- rownames(h5_seurat@meta.data)
# split sample column
h5_seurat@meta.data <- separate(h5_seurat@meta.data, col = 'sample', into = c('Patient', 'Barcode'), 
                                    sep = '_')

# calculate mitochondrial percentage
h5_seurat$mitoPercent <- PercentageFeatureSet(h5_seurat, pattern='^MT-')

saveRDS(h5_seurat,file = 'gse147528_raw.rds')
# filtering
merged_seurat_filtered_EC <- subset(h5_seurat, subset = nCount_RNA > 500 &
                                   nFeature_RNA > 500 &
                                   mitoPercent < 20)

# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered_EC <- NormalizeData(object = merged_seurat_filtered_EC)
merged_seurat_filtered_EC <- FindVariableFeatures(object = merged_seurat_filtered_EC,nfeatures=2000)
merged_seurat_filtered_EC <- ScaleData(object = merged_seurat_filtered_EC)
merged_seurat_filtered_EC <- RunPCA(object = merged_seurat_filtered_EC)
ElbowPlot(merged_seurat_filtered_EC)
merged_seurat_filtered_EC <- FindNeighbors(object = merged_seurat_filtered_EC, dims = 1:20)
merged_seurat_filtered_EC <- FindClusters(object = merged_seurat_filtered_EC,resolution=1)
merged_seurat_filtered_EC <- RunUMAP(object = merged_seurat_filtered_EC, dims = 1:20)
merged_seurat_filtered_EC <- RunTSNE(object = merged_seurat_filtered_EC, dims = 1:20)

# Create a new column named 'group' in the metadata
merged_seurat_filtered_EC@meta.data$group <- ifelse(merged_seurat_filtered_EC@meta.data$Patient %in% c('EC1', 'EC2', 'EC3'), 'braak 0',
                                                 ifelse(merged_seurat_filtered_EC@meta.data$Patient %in% c('EC4', 'EC5', 'EC6', 'EC7'), 'braak 2',
                                                        'braak 6'))
saveRDS(merged_seurat_filtered_EC,file = 'D:/KHOA_recover/Alzheimer paper-Quy/GSE147528_entorhinal/gse147528_entorhinal_filtered.rds')



#Read RDS object
merged_seurat_filtered_EC<-readRDS('D:/KHOA_recover/Alzheimer paper-Quy/GSE147528_entorhinal/gse147528_entorhinal_filtered.rds')


####  merge 2 types of neuron into 1 type , names "NEURON"
merged_seurat_filtered_EC@meta.data <- merged_seurat_filtered_EC@meta.data %>%
  mutate(cluster_new = if_else(grepl("Neuron", cluster), "Neuron", cluster))
merged_seurat_filtered_EC$cluster_new <- factor(merged_seurat_filtered_EC$cluster_new, levels = c('OPC','Astrocyte','Microglia','Oligo',
                                                                                  'Endothelial','Neuron')) 


list_check<-c('ADD2','AP2M1','ATP6V0C','ATP6V1G2','C12orf10','CACNG2','CELF4',
              'FAM174B','GLT1D1','INA','KCNQ2','KLHL35','MCAT','RPH3A',
              'SEZ6L2','SV2A','SYT3','TMEM59L','NOTCH2NL','OSBPL11')

# plot
p1 <- DimPlot(merged_seurat_filtered_EC, reduction = 'umap', group.by = 'Patient')
DimPlot(merged_seurat_filtered_EC, reduction = 'umap', group.by = 'group',
              cols = c('#13d18c','#ffcc33','#df6f75'))

DimPlot(merged_seurat_filtered_EC, reduction = 'umap', group.by = 'seurat_clusters',label = TRUE,
        label.size = 4, label.color = "black")

FeaturePlot(merged_seurat_filtered_EC, 
            reduction = 'umap', features = c('ANKRD36','BAZ1A','CGN','CHD1','FAM84B','FGD4','PIP4K2A',
                                                                    
                                            'GPAM','ITPR2','KDM5A','KIF5B','MALAT1','NOTCH2NL',
                                            'PLGLB1','SLCO4A1','SNAP23','TOB1','ZFP36L1'))
FeaturePlot(merged_seurat_filtered_EC, 
            reduction = 'umap', features = c('ALDOA','AP2M1','ATP2B3','ATP6V0C','ATP6V1G2',
                                             'C12orf10','CALY','CHRM1','CPNE6','CX3CL1','DDAH1',
                                             'DHRS7B','DOCK3','DPCD','ELAVL4','GPI','HMP19','ICA1',
                                             'INA','KCNQ2','MAP7D2','MAST3','MCAT','NAPA','NPTXR',
                                             'NRGN','PDE2A','PRSS3','REEP1','SNCA','TUBA4A','WBSCR17'))

FeatureDimPlot(merged_seurat_filtered_EC, features =c('CELF4','KCNQ2','RPH3A','SV2A','TMEM59L'),
               reduction = 'umap',theme_use = theme_blank)


b<-DimPlot(merged_seurat_filtered_EC, reduction = 'umap', group.by = 'cluster')
b

OPC = c(13,29)
Astrocyte = c(4,6,15,30,31)
Microglia = c(3,20)
Oligo=c(0,2)

Endothelial = c(23)

Excit = c(1,8,9,14,16,17,18,19,21,22,24,25,27,28)
Inhibit = c(5,7,10,11,12,26)


cl.mat <- cbind(c( rep("OPC", length(OPC)), rep("Astrocyte", length(Astrocyte)), rep("Microglia", length(Microglia)),rep("Oligo", length(Oligo)),
                   rep("Endothelial", length(Endothelial)), rep("Excitatory Neuron", length(Excit)), rep("Inhibitory Neuron", length(Inhibit))),
                c(OPC, Astrocyte, Microglia,Oligo, Endothelial, Excit, Inhibit))

cl.vec <- merged_seurat_filtered_EC$RNA_snn_res.1
ct.vec <- rep(NA, length(cl.vec))
for(x in unique(cl.mat[,1])){
  cl.x <- cl.mat[cl.mat[,1]==x,2]
  ct.vec[which(cl.vec%in%cl.x)] <- x
}
merged_seurat_filtered_EC$cluster <- ct.vec

merged_seurat_filtered_EC$cluster <- factor(merged_seurat_filtered_EC$cluster, levels = c('OPC', 'Astrocyte','Microglia','Oligo', 'Endothelial', 
                                                                                    'Excitatory Neuron', 'Inhibitory Neuron'))
DimPlot(merged_seurat_filtered_EC, reduction = "umap", label = TRUE, pt.size = 0.5,group.by ='cluster',cols='Set2',label.size = 3,label.box = T,repel = T) 

comparisons <- list(c('braak 0','braak 2'),c('braak 2','braak 6'),c('braak 0','braak 6'))
Idents(object = merged_seurat_filtered_EC)<-merged_seurat_filtered_EC$cluster_new

### Violin plot 
VlnPlot(merged_seurat_filtered_EC, features = c("KCNQ2"),pt.size = 0,
        idents = c('Neuron'),group.by = 'group', split.by = "group",cols = c('#13d18c','#ffcc33','#df6f75')) +
  stat_compare_means(comparisons = comparisons,label = "p.signif",size=5)+ geom_boxplot(width=.1)+
  ylim(0,7)+color_palette('Set1')+
  theme(
    legend.text = element_text(size = 20),      # Adjust the size of legend text
    axis.text = element_text(size = 20),        # Adjust the size of axis labels
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 20,face = 'bold'),# Hide x-axis tick labels
    # Adjust the size of axis titles
    strip.text = element_text(size = 15),       # Adjust the size of strip (facet) labels
    plot.caption = element_text(size = 20),     # Adjust the size of the plot caption
    axis.text.x = element_text(size = 20),      # Adjust the size of x-axis tick labels
    axis.text.y = element_text(size = 15),      # Adjust the size of y-axis tick labels
    strip.text.x = element_text(size = 15),
    plot.title = element_text(size = 25,face = 'bold')      # Adjust the size of x-axis facet labels
  )



####### check list new --------

listnew<-c('ATP2B2','ATP6V0C','ATP6V1G2','FKBP1B','HAS1','INA','KCNQ2','KLHL35','PRKAR1B','RPH3A','SEZ6L2','ZFP36L1')

FeatureDimPlot(merged_seurat_filtered_EC,reduction = 'umap',features = list_4,label.size = 10,theme_use = theme_blank)
FeaturePlot(merged_seurat_filtered_EC, features = listnew)




list0.7<-c('ADD2','AP2M1','ATP2B2','ATP6V0C','ATP6V1G2','C12orf10','CACNG2','CCDC157',
           'CELF4','CGREF1','FAM174B','GLT1D1','INA','KCNQ2','KLHL35','MCAT','PRKAR1B',
           'RPH3A','SEZ6L2','SV2A','SYT3','TMEM59L','NOTCH2NL','OSBPL11','RHOBTB3','ZFP36L1','ZNF217')

list0.75<-c('ADD2','AP2M1','ATP6V0C','ATP6V1G2','C12orf10','CACNG2','FAM174B','GLT1D1',
            'INA','KCNQ2','KLHL35','MCAT','RPH3A','SEZ6L2','SV2A','SYT3','TMEM59L','NOTCH2NL','OSBPL11')
list_4<-c('KCNQ2','RPH3A','SV2A','TMEM59L')
########### loop to plot violin plot
for (gene in list_4){
  plot<-VlnPlot(merged_seurat_filtered_EC, features = gene,pt.size = 0,
                idents = c('Neuron'),group.by = 'group', split.by = "group",cols = c('#13d18c','#ffcc33','#df6f75')) +
    stat_compare_means(comparisons = comparisons,label = "p.signif",size=5)+ geom_boxplot(width=.1)+
    ylim(0,7)+color_palette('Set1')+
    
    theme(
      legend.text = element_text(size = 15),      # Adjust the size of legend text
      axis.text = element_text(size = 20),        # Adjust the size of axis labels
      axis.title.x = element_blank(), 
      axis.title.y = element_text(size = 20,face = 'bold'),# Hide x-axis tick labels
      # Adjust the size of axis titles
      strip.text = element_text(size = 15),       # Adjust the size of strip (facet) labels
      plot.caption = element_text(size = 15),     # Adjust the size of the plot caption
      axis.text.x = element_text(size = 20),      # Adjust the size of x-axis tick labels
      axis.text.y = element_text(size = 15),      # Adjust the size of y-axis tick labels
      strip.text.x = element_text(size = 15),     # Adjust the size of x-axis facet labels
      plot.title = element_text(size = 25,face = 'bold') 
    )
  ggsave(filename = paste0("G:/Shared drives/BTKQ/Quy-AD-2024/AD-Figs-Ppt/SC/EC_singlecell/neuron_new/list_4gene/", gene, ".png"), plot = plot)
}



#### Dotplot -------------
library(ggmin)
DotPlot(merged_seurat_filtered_EC,group.by = "cluster",dot.scale = 5 ,features = c("PDGFRA","GFAP","DOCK8","MBP",
                                                                                   "CLDN5","SLC17A7","GAD1"),
        cols = c("red","blue"),scale = T)+scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu")))+
  ggmin::theme_min()+ RotatedAxis()+theme(axis.title = element_blank(),
                                          axis.text.x = element_text(size = 15),      # Adjust the size of x-axis tick labels
                                          axis.text.y = element_text(size = 20))


#### Enrichment ==============
FeatureStatPlot(merged_seurat_filtered_EC, stat.by = list_4, group.by = "cluster_new",
                add_box = TRUE, stack = TRUE,palette = 'Set2',legend.position = "top", legend.direction = "horizontal",
                #theme_use = theme_blank,
                theme_args=list(axis.text.x = element_text(size = 20),
                                axis.title.x = element_blank(),
                                axis.title.y = element_blank(),
                                legend.key.size = unit(0.5, 'cm'), #change legend key size
                                legend.key.height = unit(0.5, 'cm'), #change legend key height
                                legend.key.width = unit(0.5, 'cm'), #change legend key width
                                legend.title = element_blank(), #change legend title font size
                                legend.text = element_text(size=15))) #change legend text font size))









### split neuron only ================
neuron <-subset(x = merged_seurat_filtered_EC, subset = (cluster_new == "Neuron"))

DimPlot(neuron, reduction = "umap", label = F, pt.size = 0.5,group.by = 'group',
        cols= c('#13d18c','#ffcc33','#df6f75'),label.size = 3) 




### DEG ----------

Idents(object = merged_seurat_filtered_EC)<-merged_seurat_filtered_EC$cluster_new
markers.cluster = FindAllMarkers(merged_seurat_filtered_EC,logfc.threshold = 0.5,test.use = "wilcox",only.pos = F)
neuron<-markers.cluster[markers.cluster$cluster=='Neuron'&markers.cluster$p_val_adj<0.05,]
neuron <- neuron[order(neuron$avg_log2FC,decreasing=TRUE),]




#### GSEA ============
library(msigdbr)
library(fgsea)
hallmark<-msigdbr(species = "Homo sapiens", category = "C5",subcategory = 'GO:BP')
hallmark_list<-split(x = hallmark$gene_symbol,f=hallmark$gs_name)

## define input
neuron <- neuron[order(neuron$avg_log2FC,decreasing=TRUE),]
genelist <- neuron$avg_log2FC
names(genelist) <- neuron$gene

gsea.re2 <- fgsea(pathways = hallmark_list,
                  stats = genelist,
                  nperm=1000,
                  minSize=1,
                  maxSize=10000)

colnames(gsea.re2)
g2 <- gsea.re2[gsea.re2$pval<0.05,]
g2 <- g2[order(g2$NES,decreasing = T),]
hallmark_new<-hallmark[hallmark$gs_name%in%g2$pathway,]
#writexl::write_xlsx(g2,path = 'G:/Shared drives/BTKQ/Quy-AD-2024/AD-Figs-Ppt/SC/DEG_neuron/EC_neuron.xlsx')



g3<-readxl::read_xlsx("G:/Shared drives/BTKQ/Quy-AD-2024/AD-Figs-Ppt/SC/DEG_neuron/EC_neuron_GOBP.xlsx")
ggbarplot(g3, x = "pathway", y = "NES",
          fill = "group",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "aaas",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          #x.text.angle = 90, 
          #y.text.angle = 90,# Rotate vertically x axis texts
          ylab = "NES-score",
          legend.title = "DEG in Neuron",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )+
                                # ggplot2 theme+
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

