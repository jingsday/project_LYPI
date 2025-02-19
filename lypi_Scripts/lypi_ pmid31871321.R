#Murine models https://www.nature.com/articles/s41556-019-0439-6#Sec10
#SM data 
#Previous loading
# data <- load(paste0(wkdir,'NicheData10x.rda'))
# seurat_object <- get(data)
# seurat_object <- UpdateSeuratObject(object = seurat_object)
# str(seurat_object)
#DE genes either negative or pos
#seurat.markers <- FindAllMarkers(seurat_object)

library(Seurat)
library(ggplot2)
wkdir <- '/home/jing/Downloads/lypi_DATA/41556_2019_439_MOESM4_ESM/'

#Previous ourputs
#saveRDS(seurat_object,paste0(wkdir,'lypi_1321.rds'))
#write.csv(seurat.markers, paste0(wkdir,'markers.csv'))

seurat_object <- readRDS(paste0(wkdir,'lypi_1321.rds'))
str(seurat_object)
DimPlot(seurat_object, reduction = "tsne")

#Violin plot
VlnPlot(seurat_object, features=c('Vdr','Cyp27a1','Cyp27b1','Cyp2r1'))
#'Vdr','Cyp27a1','Cyp2r1'  not present 'Cyp27b1'

#tsne 
FeaturePlot(seurat_object, features = c('Cyp27a1'))

for (i in c('Vdr','Cyp27a1','Cyp27b1','Cyp2r1')){

    if (nrow(seurat.markers[seurat.markers$gene == i,]) >0){
        print(paste0(i, ' differentially expressed in'))
        print(seurat.markers[seurat.markers$gene == i,])
    } else {
        print(paste0(i, ' not differentially expressed in any clusters.'))
    }

}

#tsne gene expression 'Vdr','Cyp27a1','Cyp27b1','Cyp2r1'
#cell counts per type


sorted_counts <- as.data.frame(sort(table(seurat_object@active.ident), decreasing = TRUE))

# Rename columns for clarity
colnames(sorted_counts) <- c("identity", "count")
sorted_counts <- sorted_counts[order(toupper(sorted_counts$identity)), ]
# Plot as a horizontal bar chart

ggplot(sorted_counts, aes(x = reorder(identity, count), y = count)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(x = "Cell Type", y = "Count", title = "Cell Type Distribution") +
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12))+



#Results
# Cyp27a1 differentially expressed in"
#                  p_val avg_log2FC pct.1 pct.2    p_val_adj           cluster  gene
# Cyp27a1   4.521488e-04  -1.423644 0.004 0.030 1.000000e+00         Adipo-CAR  Cyp27a1
# Cyp27a1.1 1.350501e-17   2.454020 0.098 0.025 2.255471e-13 Arteriolar fibro.  Cyp27a1
# Cyp27a1.2 1.206314e-06   1.974297 0.131 0.027 2.014665e-02   Dendritic cells  Cyp27a1
# Cyp27a1.3 4.714634e-07   1.783411 0.067 0.026 7.873910e-03  Endosteal fibro.  Cyp27a1
# Cyp27a1.4 2.398414e-03  -4.290459 0.000 0.029 1.000000e+00      Ery/Mk prog.  Cyp27a1
# Cyp27a1.5 5.652243e-08  -4.959557 0.003 0.032 9.439811e-04     Erythroblasts  Cyp27a1
# Cyp27a1.6 3.456504e-03  -1.832881 0.010 0.030 1.000000e+00             LMPPs  Cyp27a1
# Cyp27a1.7 9.902866e-42   2.164279 0.184 0.024 1.653878e-37         Monocytes  Cyp27a1
# Cyp27a1.8 2.797270e-85   3.674008 0.371 0.024 4.671720e-81     Schwann cells  Cyp27a1
# Cyp27a1.9 3.767580e-05   2.632817 0.092 0.027 6.292235e-01    Stromal fibro.  Cyp27a1
# 
# 
# Cyp27b1 not differentially expressed in any clusters.
# Cyp2r1 differentially expressed in 
#               p_val avg_log2FC pct.1 pct.2  p_val_adj         cluster  
# Cyp2r1 2.148645e-06   1.201577 0.031 0.009 0.03588453 Gran/Mono prog. 
# 
# Vdr differentially expressed in
#               p_val avg_log2FC pct.1 pct.2     p_val_adj        cluster gene
# Vdr    6.459633e-04  -3.847567 0.008 0.024  1.000000e+00  Erythroblasts  Vdr
# Vdr.1  1.430832e-03  -4.628682 0.000 0.023  1.000000e+00 Myofibroblasts  Vdr
# Vdr.2 1.432354e-229   5.575646 0.435 0.014 2.392175e-225      Ng2+ MSCs  Vdr
# Vdr.3  3.976234e-13   3.607292 0.158 0.020  6.640708e-09      Osteo-CAR  Vdr
# Vdr.4  1.201194e-95   3.936367 0.434 0.018  2.006114e-91    Osteoblasts  Vdr


#Output as a format readable in python

library(SeuratDisk)


# Convert Seurat object to AnnData format and save
SaveH5Seurat(seurat_object, filename = paste0(wkdir,'lypi_1321_anndata.h5Seurat'))
Convert(paste0(wkdir,'lypi_1321_anndata.h5Seurat'), dest = "h5ad")

