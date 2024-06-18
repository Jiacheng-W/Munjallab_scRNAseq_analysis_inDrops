install.packages("Seurat")
install.packages("anndata")
install.packages("tidyverse")

setwd('/Users/jiacheng/Documents/scRNAseq_analysis_indrops')
setwd('E:/Jiacheng/scRNAseq_analysis_indrops')
library(Matrix)
library(anndata)
library(Seurat)
library(tidyverse)
ad <- read_h5ad("adata_ear_indexed.h5ad")

# Get the expression data matrix
data_matrix <- ad$X
data_matrix <- t(data_matrix)
# # Test mitocondrial genes
# a <- colnames(data_matrix)
# 
# mitogenes <- grepl("^mt-", a)
# print(a[mitogenes])
# 
# mito_data <- data_matrix[, mitogenes]
# nonzero_cells <- rowSums(mito_data>0)>0
# 
# numbercells <- sum(nonzero_cells)


# Get the cell metadata
cell_metadata <- ad$obs

# # Get the gene metadata
# gene_metadata <- ad$var #This is empty.

# Create a Seurat object from the data matrix
seurat_object <- CreateSeuratObject(counts = data_matrix, project = "InDrops", 
                                    min.cell = 3, min.feature = 200)
for (meta in colnames(cell_metadata)) {
  seurat_object <- AddMetaData(seurat_object, metadata = cell_metadata[[meta]], col.name = meta)
}

head(rownames(seurat_object))

#Cells that have higher than 15% mitochondrial genes are filtered out according to the manuscript.
#Double check the mito percentage by plotting here.
seurat_object[['percent.mt']] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## Subsetting the object to get only 40 hpf cells.
wt_40hpf <- subset(seurat_object, subset = nFeature_RNA >200 & nFeature_RNA < 4000)
# test_ob <- subset(seurat_object, subset = nFeature_RNA >200 & nFeature_RNA < 5000)
wt_40hpf <- subset(wt_40hpf, subset = timepoint == "40hpf" & genotype == "wildtype")
# test_ob <- subset(test_ob, subset = timepoint == "40hpf")


# Check subsetting result
table(wt_40hpf$timepoint)
#Normalization has been done according to the manuscript, but rescaling to 10000 counts per cell here.
wt_40hpf <- NormalizeData(wt_40hpf, normalization.method = "LogNormalize", scale.factor = 10000)

#Identify variable features
wt_40hpf <- FindVariableFeatures(wt_40hpf, selection.method = "vst", nfeatures = 3000)
#Check top10 variable genes.
top10 <- head(VariableFeatures(wt_40hpf), 10)
top10
#Plot variable features (optional)
plot1 <- VariableFeaturePlot(wt_40hpf)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#Scale data. (Different than the normalization step, this is scaling distributions for PCA)
all.genes <- rownames(wt_40hpf)
wt_40hpf <- ScaleData(wt_40hpf, features = all.genes)

#Reduction
#Default setting of RunPCA has 50 PCs
wt_40hpf <- RunPCA(wt_40hpf, features = VariableFeatures(object = wt_40hpf))
DimPlot(wt_40hpf, reduction = "pca")

ElbowPlot(wt_40hpf, ndims = 50)

wt_40hpf <- RunUMAP(wt_40hpf, dims = 1:20) #Choose dims based on elbowplot.
FeaturePlot(wt_40hpf, features = c('stm'), reduction = 'umap')

#Clustering
wt_40hpf <- FindNeighbors(wt_40hpf, dims = 1:20) #Louvain cluster, graph based clustering
wt_40hpf <- FindClusters(wt_40hpf, resolution = 0.5)
DimPlot(wt_40hpf, reduction = 'umap', group.by = 'seurat_clusters', label = T)

#Get DEGs for each cluster.
wt_40hpf.markers <- FindAllMarkers(wt_40hpf, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)
marker_summary <- wt_40hpf.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

print(marker_summary, n=50)





#Get subset based on stm expression (cluster 13)
earcells <- subset(wt_40hpf, idents = 13)

table(earcells$seurat_clusters)

#Identify variable features
earcells <- FindVariableFeatures(earcells, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(earcells), 10)
top10
#Plot variable features
plot1 <- VariableFeaturePlot(earcells)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#Scale data
all.genes <- rownames(earcells)
earcells <- ScaleData(earcells, features = all.genes)

#Reduction
#Default setting of RunPCA has 50 PCs
earcells <- RunPCA(earcells, features = VariableFeatures(object = earcells))
DimPlot(earcells, reduction = "pca")

ElbowPlot(earcells, ndims = 50)

earcells <- RunUMAP(earcells, dims = 1:20) #Choose dims based on elbowplot.
FeaturePlot(earcells, features = c('cyr61l1', 'tmprss5', 'lmx1ba', 'lmx1bb', 'vcana', 'vcanb'), reduction = 'umap')
FeaturePlot(earcells, features = c('rcn3', 'cnmd', 'fkbp11', 'fkbp7', 'creb3l1', 'tmem45a'), reduction = 'umap')
FeaturePlot(earcells, features = c('cyr61l1', 'thbs1b', 'cnmd', 'wnt2', 'col2a1a', 'col11a2'), reduction = 'umap')
FeaturePlot(earcells, features = c('cyr61l1', 'kcnq1', 'cnmd', 'wnt2', 'col2a1a', 'col11a2'), reduction = 'umap')

#Clustering
earcells <- FindNeighbors(earcells, dims = 1:20) #Louvain cluster, graph based clustering
earcells <- FindClusters(earcells, resolution = 0.7)
DimPlot(earcells, reduction = 'umap', group.by = 'seurat_clusters', label = T, pt.size = 3)

#Get DEGs for each cluster.
earcells.markers <- FindAllMarkers(earcells, only.pos = FALSE, min.pct = 0.25)
min(earcells.markers$avg_log2FC)

marker_summary <- earcells.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

print(marker_summary, n=16)

# Arrange the dataframe by log2 fold change within each cluster
sortedDEG_earcells <- earcells.markers %>%
  arrange(cluster, desc(avg_log2FC))

# Export DEGs
write.csv(sortedDEG_earcells, "wt_40hpf_earcells_DEG.csv", row.names = FALSE)

#Generate volcano plot.
DEG2 <- sortedDEG_earcells[sortedDEG_earcells$cluster==7,]
DEG <- FindMarkers(earcells, ident.1 = 7, min.pct = 0.25)










