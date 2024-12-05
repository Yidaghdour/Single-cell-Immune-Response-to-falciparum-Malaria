options(future.globals.maxSize = 540 * 1024^3) # Adjusted memory limit

library(Seurat) 
library(SeuratObject)
library(SeuratDisk)
library(sctransform)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(ggplot2) 
library(patchwork)
library(harmony)
library(future)

# Load Seurat object
seurat_object <- readRDS("/scratch/mv83/sc-eQTLgen-consortium/WG1-pipeline-QC_imputation/WG1-pipeline-QC_wgpipeline/RESULTS2/QC_figures/seurat_object_all_pools_singlet_barcodes_final_assignments.rds")

# Load metadata
metadata <- read.csv("subset_metadata.txt", sep = "\t")
samples_to_include <- metadata$Sample_ID  
batch_to_include <- metadata$scRNA_Batch 
pool_to_include <- metadata$scRNA_Pool        


seurat_object_subset <- subset(seurat_object, subset = (Assignment %in% samples_to_include))

# Pre-QC Plots
pag.combined <- PercentageFeatureSet(seurat_object_subset, pattern = "^MT-", col.name = 'percent.mt', assay = "RNA")
DT <- table(pag.combined$Assignment)
write.table(DT, file="RawCount_Table.txt", sep="\t")

# Violin and Feature Plots
png("Raw_ViolinPlot.png", width = 2500, height = 2500, res = 300)
VlnPlot(pag.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

png("Raw_FeaturePlot.png", width = 2500, height = 2500, res = 300)
FeatureScatter(pag.combined, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(pag.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()


# Step 1: Split Seurat object by Batch and Pool
combined_list <- SplitObject(seurat_object_subset, split.by = "Pool")



# Add mitochondrial percentage and filter
combined_list <- lapply(combined_list, function(x) {
  x <- PercentageFeatureSet(x, pattern = "^MT-", col.name = 'percent.mt')
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 & nCount_RNA > 250)
  return(x)
})


# Step 2: Normalize with SCTransform and run PCA
combined_list <- lapply(combined_list, function(x) {
  x <- SCTransform(x, verbose = FALSE, vars.to.regress = c("percent.mt", "nFeature_RNA"), return.only.var.genes = FALSE)  # Normalize all genes
  x <- RunPCA(x, verbose = FALSE)  # PCA
  return(x)
})

# Step 2.5: Select integration features that exist across all objects
features <- SelectIntegrationFeatures(object.list = combined_list, nfeatures = 2000)

# Verify features are in scale.data for each object in combined_list
validated_features <- features
for (i in seq_along(combined_list)) {
  available_features <- rownames(GetAssayData(combined_list[[i]], slot = "scale.data"))
  validated_features <- intersect(validated_features, available_features)
}

# Step 3: Feature Selection & Integration
combined_list <- PrepSCTIntegration(object.list = combined_list, anchor.features = validated_features)


plan()

plan("multicore", workers = 20)

plan()


features <- SelectIntegrationFeatures(object.list = combined_list, nfeatures = 2000)

# Ensure selected features are consistent across objects
validated_features <- features
for (i in seq_along(combined_list)) {
  available_features <- rownames(GetAssayData(combined_list[[i]], assay = "SCT", slot = "scale.data"))
  validated_features <- intersect(validated_features, available_features)
}


# Merge Seurat objects without integration to prepare for Harmony
seurat_combined <- merge(combined_list[[1]], y = combined_list[-1])

VariableFeatures(seurat_combined) <- validated_features

# Set Default Assay and Run PCA
DefaultAssay(seurat_combined) <- "SCT"
seurat_combined <- RunPCA(seurat_combined, npcs = 50, verbose = FALSE)

# Save PCA output (optional)
saveRDS(seurat_combined, file = "RunPCA.Rds")

# Run Harmony for batch correction
seurat_combined <- RunHarmony(
    object = seurat_combined, 
    group.by.vars = "Pool",
    assay.use = "SCT",
    reduction.use = "pca",
    dims = 1:40,
    verbose = FALSE,
    theta = 2
)

# UMAP and Clustering using Harmony embeddings
seurat_combined <- RunUMAP(seurat_combined, reduction = "harmony", dims = 1:40, min.dist = 0.5, n.neighbors = 50L)
seurat_combined <- FindNeighbors(seurat_combined, reduction = "harmony", dims = 1:40)

saveRDS (seurat_combined, file="FindNeighbours.Rds")

# Optional parallel processing
plan("multicore", workers = 20)

# Clustering with the resolution of 0.6
seurat_combined <- FindClusters(seurat_combined, resolution = 0.8)


saveRDS (seurat_combined, file="Harmony_0.8_40Dims.Rds")


#### Annotating with Reference ######

reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")

DefaultAssay(seurat_combined) <- 'SCT'

anchors.Harmony <- FindTransferAnchors(
  reference = reference,
  query = seurat_combined,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50,
  recompute.residuals = FALSE
)


seurat_combined <- MapQuery(
  anchorset = anchors.Harmony,
  query = seurat_combined,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca"
)


saveRDS (seurat_combined, file="Harmony_0.8_40Dims_Annotated.Rds")


## Adding metadata 
existing_metadata <- seurat_combined@meta.data
metadata <- read.csv("metadata_toadd_Tala.csv", header = TRUE, stringsAsFactors = FALSE)


merged_metadata <- merge(
  existing_metadata, 
  metadata, 
  by = c("Assignment"), 
  all.x = TRUE, 
  sort = FALSE
)
head(existing_metadata)
seurat_combined@meta.data <- merged_metadata
row.names(seurat_combined@meta.data) <- seurat_combined@meta.data$Barcode


Idents(seurat_combined) <- "seurat_clusters"

# UMAP colored by PoolID
png("UmapPlot_Clusture.png", width = 2500, height = 2500,res=300)
DimPlot(seurat_combined, reduction = "umap", group.by = "seurat_clusters",label = TRUE, pt.size = 0.5)
dev.off()

Idents(seurat_combined) <- "predicted.celltype.l1"

# UMAP colored by cluster
png("UmapPlot_predicted.celltype.l1.png", width = 2500, height = 2500,res=300)
DimPlot(seurat_combined, reduction = "umap", group.by = "predicted.celltype.l1",label = TRUE, pt.size = 0.5)
dev.off()

png("UmapPlot_predicted.celltype.l2.png", width = 2500, height = 2500,res=300)
DimPlot(seurat_combined, reduction = "umap", group.by = "predicted.celltype.l2",label = TRUE, pt.size = 0.5)
dev.off()


png("UmapPlot_Ethnicity.png", width = 2500, height = 2500,res=300)
DimPlot(seurat_combined, reduction = "umap", group.by = "Ethnicity")
dev.off()


png("Filtered_ViolinPlot.png", width = 2500, height = 2500,res=300)
VlnPlot(seurat_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

png("Filtered_FeaturePlot.png", width = 2500, height = 2500,res=300)
FeatureScatter(seurat_combined, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter(seurat_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

DT <-table(seurat_combined$Assignment)
write.table(DT,file="FilteredCount_Table_PerSample.txt",sep="\t")


DT <-table(seurat_combined$Ethnicity)
write.table(DT,file="FilteredCount_Table_Ethnicity.txt",sep="\t")


cluster_markers <- FindAllMarkers(seurat_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, file = "all_cluster_markers.csv", row.names = FALSE)

## Highy differentaited Markers
refined_markers <- subset(cluster_markers, avg_log2FC > 0.5 & p_val_adj < 0.05)
write.csv(refined_markers, file = "refined_cluster_markers.csv", row.names = FALSE)

png("Top_markeres_HeatmapPlot.png", width = 2500, height = 2500,res=300)
DoHeatmap(seurat_combined, features = refined_markers$gene)
dev.off()

png("Top_markeres_DotPlot.png", width = 2500, height = 2500,res=300)
DotPlot(seurat_combined, features = refined_markers$gene) + RotatedAxis()
dev.off()

png("Top_markeres_VlnPlot.png", width = 2500, height = 2500,res=300)
VlnPlot(seurat_combined, features = unique(refined_markers$gene), ncol = 5)
dev.off()

png("Top_markeres_FeaturePlot.png", width = 2500, height = 2500,res=300)
FeaturePlot(seurat_combined, features = unique(refined_markers$gene), ncol = 5)
dev.off()








