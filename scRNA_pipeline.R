setwd("T:/scRNA_brain_dataset")

library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)

install.packages("hdf5r")
library(hdf5r)

#human brain cancer tissue
# Load the dataset
hbct.m <- Read10X_h5(filename = 'CytAssist_11mm_FFPE_Human_Glioblastoma_raw_feature_bc_matrix.h5')
str(hbct.m)

#Create seurat object
hbct <- CreateSeuratObject(
  counts = hbct.m,
  project = "Glioblastoma_scRNAseq",
  min.cells = 3,
  min.features = 200)
str(hbct)
view(hbct)      #18993 features across 12949 samples
view(hbct@meta.data)

# 1. Quality control
#Add mitochondrial percentage
hbct[["percent.mt"]] <- PercentageFeatureSet(hbct,pattern = "^MT-")

plots <- VlnPlot(hbct,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),combine = FALSE)

wrap_plots(plots, ncol = 3)
view(wrap_plots)

# nFeature_RNA: Most cells show a moderate number of genes (good quality); very low or high values indicate poor cells or doublets
# nCount_RNA: Most cells have normal UMI counts; a few high-count cells likely represent doublets
# percent.mt: Most cells have low mitochondrial content (healthy); higher values indicate stressed or dying cells


# 2. Filtering
hbct <- subset(hbct,subset = nFeature_RNA > 500 &
    nFeature_RNA < 8000 &
    nCount_RNA < 50000 &
    percent.mt < 15)

# 3. Normalization
hbct <- NormalizeData(hbct,
  normalization.method = "LogNormalize",scale.factor = 10000)

# 4. Identify HVGs
hbct <- FindVariableFeatures(hbct,
  selection.method = "vst",nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(hbct), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hbct)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 5. Scaling
hbct <- ScaleData(hbct,features = VariableFeatures(hbct))

# 6. Dimentionality Reduction
hbct <- RunPCA(hbct,features = VariableFeatures(object = hbct))
ElbowPlot(hbct, ndims = 50)

# 7. Clustering cells
hbct <- FindNeighbors(hbct, dims = 1:12)
hbct <- FindClusters(hbct, resolution = 0.5)
table(Idents(hbct))


#UMAP
hbct <- RunUMAP(hbct, dims = 1:12)
library(ggplot2)
DimPlot(hbct, reduction = "umap", label = TRUE, pt.size = 0.5)

# Find positive marker genes for all clusters
markers <- FindAllMarkers(hbct,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)

# View top 10 markers per cluster
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

head(top_markers)
View(top_markers)

# Export the top 10 markers per cluster to a CSV file
write.csv(top_markers, file = "top_markers_per_cluster.csv",row.names = FALSE)

# Select the top 3 marker genes per cluster based on average log2 fold change
top3_genes <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 3)

# Export the top 3 markers per cluster to a CSV file
write.csv(top3_genes, file = "top3_markers_cluster.csv",row.names = FALSE)

# Save the Seurat object
saveRDS(hbct, file = "Glioblastoma_scRNAseq_seurat.rds")



#----------Rename clusters based on marker interpretation
hbct <- RenameIdents(hbct,`0`  = "Plasma_B",`1`  = "Fibroblast_Stromal",
                     `2`  = "Inflammatory_Immune",`3`  = "Activated_T_cells",
                     `4`  = "Neuronal",`5`  = "Erythroid_Stress",
                     `6`  = "Neuronal_like",`7`  = "Hypoxic_Tumor",
                     `8`  = "Mature_Neurons",`9`  = "Oligodendrocytes",
                     `10` = "Endothelial",`11` = "Vascular_SMC",
                     `12` = "Tumor_Associated",`13` = "Interneurons")


#------------annotated UMAP
DimPlot(hbct, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Cell populations in Glioblastoma")

# View all the cluster identities
levels(Idents(hbct))

# Find marker genes for all clusters and export them to a CSV file
write.csv(FindAllMarkers(hbct, only.pos = TRUE),"Final_Annotated_Markers.csv",row.names = FALSE)




