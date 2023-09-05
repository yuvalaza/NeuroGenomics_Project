library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratObject)
library(Matrix)
library(ggplot2)
library(cowplot)


exp_data<-read.csv("~/NeuroGenom/Final Project/Data/BreastCancerExpressionInCells.csv",row.names=1)
exp_data<-t(as.matrix(exp_data))
brst<- CreateSeuratObject(counts = exp_data, project = "breast", min.cells =1, min.features =20)
brst[["percent.mt"]] <- PercentageFeatureSet(brst, pattern = "^MT-")
sum(brst[["percent.mt"]])#Zero!

VlnPlot(brst, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(brst, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#Filtering Cells By Unique Features
brst <- subset(brst, subset = nFeature_RNA > 25 & nFeature_RNA < 225)

brst <- NormalizeData(brst)

brst <- FindVariableFeatures(brst, selection.method = "vst", nfeatures = 200)
top5 <- head(VariableFeatures(brst), 5)
plot1 <- VariableFeaturePlot(brst)
plot2 <- LabelPoints(plot = plot1, points = top5, repel = TRUE)
plot2

genes_names <- rownames(brst)
brst <- ScaleData(brst, features = genes_names)

brst <- RunPCA(brst, features = VariableFeatures(object = brst), npcs = 100)
print(brst[["pca"]], dims = 1:5, nfeatures = 5)

DimHeatmap(brst, dims = 1:20, cells = 500, balanced = TRUE)

brst <- JackStraw(brst, num.replicate = 100)
brst <- ScoreJackStraw(brst, dims = 1:20)
JackStrawPlot(brst, dims = 1:20,xmax = 0.1,ymax = 0.5)

#VizDimLoadings(brst, dims = 1:20, reduction = "pca")
#ElbowPlot(brst)



#examine different resolutions
plot_list <- list()

for(resolution in seq(0.3, 1.1, by = 0.1)) {
  brst <- FindNeighbors(brst, dims = 1:14,k.param=5)
  brst <- FindClusters(brst, resolution = resolution)
  plot <- DimPlot(brst, group.by = "seurat_clusters") +
    ggtitle(paste("Resolution:", round(resolution, 1)))
  plot_list[[length(plot_list) + 1]] <- plot
}

final_plot <- plot_grid(plotlist = plot_list, ncol = 3) 
print(final_plot)


brst <- RunUMAP(brst, dims = 1:14)
brst <- FindClusters(brst, resolution = 0.42)
dim_plot<-DimPlot(brst, group.by = "seurat_clusters")
dim_plot
#order by avg_log2FC
brst.markers <- FindAllMarkers(brst, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
brst.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

top_genes <- brst.markers %>%
  group_by(cluster) %>%
  slice_max(n =20, order_by = avg_log2FC)

for(i in 0:(length(levels(brst))-1)){
  print(paste("Cluster number: ", i))
  print(top_genes[top_genes$cluster==i,],n=20)
}
#order by p_val_adj 

top_genes <- brst.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = -p_val_adj)

print(top_genes, n = 60)

#cluster0.markers <- FindMarkers(brst, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

feature_plot1 <- FeaturePlot(brst, features = c("CD3G", "FOXP3", "CD8A", "CD3D", "CD3E"))
feature_plot2 <- FeaturePlot(brst, features = c("HLA.DRA","CD68", "CD4"))
feature_plot3 <- FeaturePlot(brst, features = c("IGHG1", "IGHG4", "IGKC", "IGHM"))
feature_plot4.1 <- FeaturePlot(brst, features = c("EGFR", "GRB7", "ERBB2", "PGR", "CD44","CD24"))
feature_plot4.2 <- FeaturePlot(brst, features = c( "ALDH1A3", "EPCAM", "KRT19", "KRT18","CDH1"))
feature_plot5 <- FeaturePlot(brst, features = c("HSPG2","SULF1"))

#see each gene on the dim plot
combined_plot <-feature_plot1+dim_plot
combined_plot
combined_plot <-feature_plot2+dim_plot
combined_plot
combined_plot <-feature_plot3+dim_plot
combined_plot
combined_plot <-feature_plot4.1+dim_plot
combined_plot
combined_plot <-feature_plot4.2+dim_plot
combined_plot
combined_plot <-feature_plot5+dim_plot
combined_plot




# VlnPlot(brst, features = c("CD3G", "FOXP3", "CD8A", "CD3D", "CD3E"))
# VlnPlot(brst, features = c("HLA.DRA","CD68", "CD4"))
# VlnPlot(brst, features =  c("IGHG1", "IGHG4", "IGKC", "IGHM"))
# VlnPlot(brst, features =c("EGFR", "GRB7", "ERBB2", "PGR", "CD44","CD24", "ALDH1A3", "EPCAM", "KRT19", "KRT18","CDH1"))
# VlnPlot(brst, features =  c("HSPG2","SULF1"))



cluster_names <- c("Fibroblas", "Tumor", "B_cells", "T_cells", "T_cells", "Macrophage")
names(cluster_names) <- levels(brst)
brst <- RenameIdents(brst, cluster_names)
DimPlot(brst,label = TRUE, pt.size = 0.5) + NoLegend()

#Q1
immune_clusters <- c("B_cells", "T_cells", "Macrophage")
immune_cells <- sum(Idents(brst) %in% immune_clusters)
percentage_immune <- (immune_cells / ncol(brst)) * 100
print(paste("The sample contains", percentage_immune, "% immune cells."))


#clusters proportion check
# Define marker gene groups
marker_genes <- list(
  T_cells = c("CD3G", "FOXP3", "CD8A", "CD3D", "CD3E"),
  Macrophage = c("HLA.DRA","CD68", "CD4"),
  B_cells =c("IGHG1", "IGHG4", "IGKC", "IGHM"),
  Tumor =c("EGFR", "GRB7", "ERBB2", "PGR", "CD44","CD24", "ALDH1A3", "EPCAM", "KRT19", "KRT18","CDH1"),
  Fibroblast = c("HSPG2","SULF1")
)

proportion_df <- data.frame(Group = character(),
                            Cluster = character(),
                            Proportion = numeric(),
                            Count = integer(),
                            Total = integer(),
                            stringsAsFactors = FALSE)

# Extracting the expression matrix
expr_matrix <- brst@assays$RNA@counts

for (group in names(marker_genes)) {
  group_genes <- marker_genes[[group]]
  
  for (cluster in clusters) {
    # Find the indices of the cells in this cluster
    cluster_indices <- which(brst$seurat_clusters == cluster)
    cluster_data <- expr_matrix[group_genes, cluster_indices]
    count <- sum(rowSums(cluster_data > 0) > 0)
    total <- length(cluster_indices)
    proportion <- count / total
    
    proportion_df <- rbind(proportion_df, data.frame(Group = group,
                                                     Cluster = cluster,
                                                     Proportion = proportion,
                                                     Count = count,
                                                     Total = total,
                                                     stringsAsFactors = FALSE))
  }
}
print(proportion_df)


#Q2
spatial_data<-read.csv("~/NeuroGenom/Final Project/Data/BreastCancerLocationOfCells.csv")
#in order to filter the cells Seurat filtered from the spatial data
retained_cells <- as.numeric(names(Idents(brst)))
all_cells <- 1:2748
missing_cells <- setdiff(all_cells, retained_cells)
spatial_data_filtered <- spatial_data[-missing_cells,]
spatial_data_filtered$clusters<-Idents(brst)

colnames(spatial_data_filtered)[2] <- "X_Coordinate"
colnames(spatial_data_filtered)[3] <- "Y_Coordinate"


ggplot(spatial_data_filtered, aes(x = X_Coordinate, y = Y_Coordinate, color = as.factor(clusters))) +
  geom_point() +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) + # You can adjust the colors
  labs(
    title = "Spatial Distribution of Cells",
    x = "X Coordinate (pixels)",
    y = "Y Coordinate (pixels)",
    color = "Cell Type"
  )
ggplot(spatial_data_filtered, aes(x = X_Coordinate, y = Y_Coordinate, color = as.factor(clusters))) +
  geom_point() +
  scale_color_manual(values = c("red", "blue", "green", "green", "green")) + # You can adjust the colors
  labs(
    title = "Spatial Distribution of Cells",
    x = "X Coordinate (pixels)",
    y = "Y Coordinate (pixels)",
    color = "Cell Type"
  )

#Q3-PD-L1
FeaturePlot(brst, features = c("CD274"))
VlnPlot(brst, features =  c("CD274"))
#Find PD-L1 expression in each cell
pdl1_expression <- brst[["RNA"]]@counts["CD274",]
count_pdl1_cells <- sum(pdl1_expression!=0)
percentage_pdl1 <- (count_pdl1_cells / length(pdl1_expression)) * 100

print(paste(percentage_pdl1, "% of the cells in the biopsy express the gene for PD-L1"))

