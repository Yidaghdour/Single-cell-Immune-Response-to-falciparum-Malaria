library(Seurat)
library(SeuratDisk)
library(DESeq2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(muscat)
library(scater)
library(UpSetR)
library(CATALYST)
library(purrr)
library(limma)
library(tidyverse)
library(ggrepel)
library(readxl)
library(dendextend)
library(circlize)
library(GraftingScripts)
library(ComplexHeatmap)
library(stats)
library(writexl)
library(ggpubr)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(car)

#### Create pseudobulk matrices and perform differential state analysis ####

# Read the RDS seurat object file
setwd("~/Documents/file_location/")
GData<-readRDS("Seurat_filtered_126.RDS")

#Assign cell subtype identity classes
Idents(GData) <- GData$predicted.celltype.l2
Idents(GData)

#Subset cells that belong to infected "Inf" or non-infected "NI" samples 
expr <- FetchData(GData, vars = 'Infection')
GData <- GData[, which(expr == "Inf")]
#GData <- GData[, which(expr == "NI")]

#Set up single cell experiment object for pseudobulk analysis
GData.SE<-as.SingleCellExperiment(GData)
GData.SE$cluster_id <- GData.SE$predicted.celltype.l2
GData.SE$sample_id <- GData.SE$Assignment
GData.SE$group_id <- GData.SE$Ethnicity
GData.SE$id <-  GData.SE$Assignment
(GData.SE <- prepSCE(GData.SE, 
                     kid = "cluster_id", # subpopulation assignments
                     gid = "group_id",  # group IDs (Mossi/Fulani)
                     sid = "sample_id",   # sample IDs 
                     drop = TRUE))

nk <- length(kids <- levels(GData.SE$cluster_id))
ns <- length(sids <- levels(GData.SE$sample_id))
ng <- length(gids <- levels(GData.SE$group_id))
names(kids) <- kids; names(sids) <- sids; names(gids) <- gids
assayNames(GData.SE)

#Get and save number of cells per cluster and per sample
CellInfo<-t(table(GData.SE$cluster_id, GData.SE$sample_id))
write.table(CellInfo,file="Cells_PerCelltype_PerSample.csv",sep = ",")
CellInfo

#Aggregate single cell data (sum of counts) by cluster and by sample_ID into pseudobulk data
pb <-NULL
pb <- aggregateData(GData.SE, assay = "counts", fun = "sum",  by = c("cluster_id", "sample_id"))

# write one sheet per cell subtype/cluster
CellTypes<-assayNames(pb)
for(i in 1:length(CellTypes)) {# Head of for-loop
  Fname<-NULL
  ROutput<-NULL
  Fname<-paste(CellTypes[i],"_Count_Sum.csv",sep=",")
  ROutput<-t((assay(pb[],i)))
  write.table(ROutput,file=Fname,sep = ",",row.names = T,col.names=T)
}

# construct design & contrast matrix by adding metadata (age and sex)
Metadata <- metadata(GData.SE)$experiment_info 
meta <- read.csv("metadata_for_mm.csv")
meta_subset <- meta[,c('Assignment','Age','Sex')]
#Merge the Metadata
df3 <- left_join(Metadata, meta_subset, by=c('sample_id'='Assignment'))
GData.SE@metadata$experiment_info<-df3
ei <- metadata(GData.SE)$experiment_info
#create model matrix to compare ethnic groups and account for age and sex
mm <- model.matrix(~ 0 + ei$group_id + ei$Sex + ei$Age)
dimnames(mm) <- list(ei$sample_id, c("Fulani", "Mossi", "Male", "Age"))
contrast <- makeContrasts("Fulani-Mossi", levels = mm)

# run differential state analysis
GData.SE$group_id <- factor(GData.SE$group_id,levels=c('Fulani','Mossi'))
res <- pbDS(pb, method="DESeq2", design = mm, contrast = contrast, verbose = FALSE)

# one data.frame per cluster
tbl <- res$table[[1]]
tbl
CellType<-names(tbl)
CellType
for(i in 1:length(CellType)) {# Head of for-loop
  Fname<-NULL
  ROutput<-NULL
  Fname<-paste(CellType[i],"_Sample-level_Deseq2.csv",sep="")
  ROutput<-tbl[[i]]
  write.table(ROutput,file=Fname,sep = ",",row.names = F)
  print(Fname)                    
}
# filter FDR < 5%, abs(logFC) > 0.263 & sort by adj. p-value
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 0.263)
  dplyr::arrange(u, p_adj.loc)
})

for(i in 1:length(CellType)) {# Head of for-loop
  Fname<-NULL
  ROutput<-NULL
  Fname<-paste(CellType[i],"_Sample-level_Deseq2_Filterd_pval_0.05_logFC_0.263.csv",sep="")
  ROutput<-tbl_fil[[i]]
  write.table(ROutput,file=Fname,sep = ",",row.names = F)
  print(Fname)                    
}


#### Run DeSeq2 on count matrices per cell subtype to get normalized VST counts for heatmaps and violin plots ####

# Read the sum count matrices generated from the aggregation of cell for each cell subtype from line 68 above for the infected and non-infected counterparts, and combine the two matrices 
# The following code will be applied on the CD14 monocyte cell subtype as an example
CD14_counts_INF <- as.matrix(t(read.csv("CD14 Mono_INF_Count_Sum.csv")))
CD14_counts_NI <- as.matrix(t(read.csv("CD14 Mono_NI_Count_Sum.csv")))
CD14_counts <- cbind(CD14_counts_INF, CD14_counts_NI)

# Read the metadata file and match the sorting of the sample IDs
ei_dds <- read.csv("metadata.csv")
rownames(ei_dds) <- ei_dds[,1]
ei_dds <- ei_dds[,-1]
idx <- match(rownames(ei_dds), colnames(CD14_counts))
CD14_counts <- CD14_counts[,idx]

#Split matrix into 4 groups per ethnicity and infection status and filter to keep genes with more than 5 reads in at least half of the samples in that group
# For the Infected Fulani group of samples (n = 35)
CD14_Counts_FuINF <- CD14_counts[,1:35]
# Calculate the threshold for the number of elements, half the number of samples in that group
threshold <- ncol(CD14_Counts_FuINF) * 0.5
# Function to check if a gene has greater than 5 counts (x) in at least half of the number of samples (threshold)
more_than_5 <- function(x) {
  sum(x > 5) >= threshold
}
# Filter rows with genes containing counts greater than 5 in at east half of the number of samples
filtered_rows <- CD14_Counts_FuINF[apply(CD14_Counts_FuINF, 1, more_than_5), ]
# Subset the filtered rows from the total count matrix of that group of samples
keep <- rownames(filtered_rows)
CD14_Counts_FuINF <- CD14_Counts_FuINF[keep, ]

# For the Infected Mossi group of samples (n = 52)
CD14_Counts_MoINF <- CD14_counts[,36:87]
threshold <- ncol(CD14_Counts_MoINF) * 0.5
filtered_rows <- CD14_Counts_MoINF[apply(CD14_Counts_MoINF, 1, more_than_5), ]
keep <- rownames(filtered_rows)
CD14_Counts_MoINF <- CD14_Counts_MoINF[keep, ]

# For the Non-infected Fulani group of samples (n = 23)
CD14_Counts_FuNI <- CD14_counts[,88:110]
threshold <- ncol(CD14_Counts_FuNI) * 0.5
filtered_rows <- CD14_Counts_FuNI[apply(CD14_Counts_FuNI, 1, more_than_5), ]
keep <- rownames(filtered_rows)
CD14_Counts_FuNI <- CD14_Counts_FuNI[keep, ]

# For the Non-infected Mossi group of samples (n = 16)
CD14_Counts_MoNI <- CD14_counts[,111:126]
threshold <- ncol(CD14_Counts_MoNI) * 0.5
filtered_rows <- CD14_Counts_MoNI[apply(CD14_Counts_MoNI, 1, more_than_5), ]
keep <- rownames(filtered_rows)
CD14_Counts_MoNI <- CD14_Counts_MoNI[keep, ]

# Get all gene names kept from all 4 groups of samples 
rownames_CD14_Counts_FuINF <- rownames(CD14_Counts_FuINF)
rownames_CD14_Counts_MoINF <- rownames(CD14_Counts_MoINF)
rownames_CD14_Counts_FuNI <- rownames(CD14_Counts_FuNI)
rownames_CD14_Counts_MoNI <- rownames(CD14_Counts_MoNI)
union <- Reduce(union, list(rownames_CD14_Counts_FuINF, rownames_CD14_Counts_MoINF, rownames_CD14_Counts_FuNI, rownames_CD14_Counts_MoNI))

# Filter each of the infected and non-infected count matrices for genes we decided to keep, and generate the corresponding metadata dataframes for each infection group
CD14_counts_INF_filtered <- CD14_counts_INF[union, ]
CD14_counts_NI_filtered <- CD14_counts_NI[union, ]
ei_dds_INF <- ei_dds[1:87, ]
ei_dds_NI <- ei_dds[88:126, ]
idx1 <- match(rownames(ei_dds_INF), colnames(CD14_counts_INF_filtered))
CD14_counts_INF_filtered <- CD14_counts_INF_filtered[,idx1]
idx2 <- match(rownames(ei_dds_NI), colnames(CD14_counts_NI_filtered))
CD14_counts_NI_filtered <- CD14_counts_NI_filtered[,idx2]

# Add 1 to each count in the matrix to remove zeros
CD14_counts_INF_filtered_1 <- as.matrix(CD14_counts_INF_filtered)+1
CD14_counts_NI_filtered_1 <- as.matrix(CD14_counts_NI_filtered)+1

# Perform DeSeq2 to compare ethnic groups per infection state, accounting for age and sex, here is the example for the infected group
dds <- DESeqDataSetFromMatrix(countData = CD14_counts_INF_filtered_1, colData = ei_dds_INF, design = ~ Sex + Age + group_id)
dds <- estimateSizeFactors(dds)
dds
norm_counts <- counts(dds, normalize = TRUE)
vsd <- vst(dds, blind = TRUE)

#Generate PCA plot based on all genes included 
PCA <- DESeq2::plotPCA(vsd, intgroup = "group_id", ntop= 4781)
pdf("./PCAplot_Inf.pdf")
print(PCA)
dev.off()

# Save normalized variant stabilizing transformation dataframe
vst1 <- vsd@assays@data@listData[[1]]
vst1 <- as.data.frame(vst1)
write.csv(vst1, file = "VST_INF.csv")

# Run DeSeq2 comparing Fulani to Mossi 
dds <- DESeq(dds)
result <- results(dds, contrast = c("group_id", "Fulani", "Mossi"), alpha = 0.05)
plotMA(result, ylim = c(-8,8))
results_all <- data.frame(result)
write.csv(results_all, file = "manual_deseq_result_CD14_INF.csv")

# Filter for deseq2 results with padj < 0.05 and |Log2FC| > 0.263 and save the results
sigs <- na.omit(results_all)
sigs <- sigs[sigs$padj < 0.05,]
df <- as.data.frame(sigs)
df.top <- df[ (abs(df$log2FoldChange) > 0.263),]
df.top <- df.top[order(df.top$log2FoldChange, decreasing = TRUE),]
write.csv(df.top, file = "Man_CD14_INF_Padj_0.1_logFC_0.263.csv")


# Perform the same analysis for the non-infected group

dds <- DESeqDataSetFromMatrix(countData = CD14_counts_NI_filtered_1, colData = ei_dds_NI, design = ~ Sex + Age + group_id)
dds <- estimateSizeFactors(dds)
dds
norm_counts <- counts(dds, normalize = TRUE)
vsd <- vst(dds, blind = TRUE)

PCA <- DESeq2::plotPCA(vsd, intgroup = "group_id", ntop= 4781)
pdf("./PCAplot_NI.pdf")
print(PCA)
dev.off()

vst2 <- vsd@assays@data@listData[[1]]
vst2 <- as.data.frame(vst2)
write.csv(vst1, file = "VST_NI.csv")

dds <- DESeq(dds)
result <- results(dds, contrast = c("group_id", "Fulani", "Mossi"), alpha = 0.05)
plotMA(result, ylim = c(-8,8))
results_all <- data.frame(result)
write.csv(results_all, file = "manual_deseq_result_CD14_NI.csv")

sigs <- na.omit(results_all)
sigs <- sigs[sigs$padj < 0.05,]
df <- as.data.frame(sigs)
df.top <- df[ (abs(df$log2FoldChange) > 0.263),]
df.top <- df.top[order(df.top$log2FoldChange, decreasing = TRUE),]
write.csv(df.top, file = "Man_CD14_NI_Padj_0.1_logFC_0.263.csv")

# Combined variant stabilized trasnformation matrix
vst <- cbind(vst1, vst2)
write.csv(vst, file = "vst_all.csv")


# Make heatmaps with top 25 upregulated and top 25 downregulated genes in Fulani vs Mossi

# Create gene list for DEGs in each of the comparisons (Infected and Non-infected groups)
gene_list_INF <- read.csv("Man_CD14_INF_Padj_0.05_logFC_0.263.csv")
gl_INF <- gene_list_INF[,1]
samples_INF <- colnames(vst1)
gene_list_NI <- read.csv("Man_CD14_NI_Padj_0.05_logFC_0.263.csv")
gl_NI <- gene_list_NI[,1]
samples_NI <- colnames(vst2)

# Below is the code only for CD14 infected group comparison of Fulani vs Mossi samples as an example
#Subset from the normalized vst matrix the genes from the gene list
CD14mono_mat_INF <- vst1[gl_INF, samples_INF]
CD14mono_mat_INF <- na.omit(CD14mono_mat_INF)
base_mean <- rowMeans(CD14mono_mat_INF)
#center and scale each column (Z-score) then transpose
CD14_mat_scaled_INF <- t(apply(CD14mono_mat_INF, 1, scale)) 
colnames(CD14_mat_scaled_INF)<-colnames(CD14mono_mat_INF)
rownames <- gsub("\\..*", "", rownames(CD14_mat_scaled_INF))
rownames(CD14_mat_scaled_INF) <- rownames

# Top 25 
num_keep <- 25
rows_keep_INF <- c(seq(1:num_keep), seq((nrow(CD14_mat_scaled_INF)-num_keep), nrow(CD14_mat_scaled_INF)) )

# Read in the metadata to layer onto the heatmap for infected samples
ei <- read.csv("metadata.csv")
ei_INF <- ei[1:87, ]
d <- as.data.frame(ei_INF[,c("sample_id", "group_id", "Infection", "EthInf")])
d <- d[order(d$group_id, decreasing = F),]
d <- d[order(d$Infection, decreasing = F),]
new_order <- d$sample_id
CD14_mat_scaled_INF <- CD14_mat_scaled_INF[, new_order]
EthInf_info <- data.frame(EthInf = d$EthInf)
EthInf_info
EthInf <- EthInf_info[,1]
EthInf.colors <- c(rep(c("red","green", "blue", "yellow")))
EthInf.colors
names(EthInf.colors) <- paste(c("FuInf", "MoInf", "FuNI","MoNI"),sep="")
EthInf.colors

ha <- HeatmapAnnotation(df = EthInf_info, col = list(EthInf = EthInf.colors))
dend1 = cluster_within_group(CD14_mat_scaled_INF, EthInf)

# Read in the log2 fold changes, subsetting the genes of interest from the gene list, to layer onto the heatmap
LFC_mat_INF <- read.csv("Man_CD14_INF_Padj_0.05_logFC_0.263.csv")
LFC_mat_INF <- LFC_mat_INF[order(LFC_mat_INF$log2FoldChange, decreasing = TRUE),]
rownames(LFC_mat_INF) <- LFC_mat_INF[,1]
LFC_mat_INF <- LFC_mat_INF[gl_INF,]
LFC_mat_INF <- LFC_mat_INF[1:2374,]
l2_val <- as.matrix(LFC_mat_INF[rows_keep_INF,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val)<-"logFC"
#maps values between b/w/r for min and max l2 values
col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red")) 

# Create heatmap h1 with top 25 up and down-regulated genes
h1 <- Heatmap(CD14_mat_scaled_INF[rows_keep_INF,], cluster_rows = F, top_annotation = ha,
              name="Z-score",
              row_names_gp = gpar(fontsize = 8),
              column_names_side = "top",
              show_column_names = F,
              cluster_columns = dend1,
              column_split = 8)
pdf("CD14_INF_heatmap_cluster_padj-0.05_VST.pdf", width = 14, height = 16)
print(h1)
dev.off()

# Create column to add to heatmap with the log2 fold changes of the top 25 genes
h1_lfc <- Heatmap(l2_val, row_labels = rownames(CD14_mat_scaled_INF[rows_keep_INF,]),
                  row_names_gp = gpar(fontsize = 4),
                  cluster_rows = F, name="logFC", col = col_logFC,
                  heatmap_width = unit(2, "cm"),
                  cell_fun = function(j, i, x, y, w, h, col) {
                    grid.text(round(l2_val[i, j],2), x, y, gp = gpar(fontsize = 10))
                  })

# Combine both to form the final heatmap
h_INF <- h1 + h1_lfc
pdf("CD14_INF_LFCheatmap_cluster_padj-0.05_VST.pdf", width = 12, height = 16)
print(h_INF)
dev.off()


# Create violin plots for comparison of normalized expression values for each gene between the 4 groups and perform t test
# Read in the variant stabilized transformed matrix from the ethnic comparison for both infected and non-infected samples
vst <- read.csv("vst_all.csv")
rownames(vst) <- vst[,1]
vst <- vst[, -1]
vst_t <- t(vst)
# Merge normalized expression matrix with metadata 
df1 <- merge(vst_t, ei_dds, by = "row.names")
rownames(df1) <- df1[,1]
df1 <- df1[, -1]
# Violin plot for IL6 gene for each of the 4 groups
gene.plot <- ggplot(df1, aes(x = factor(EthInf, levels = c("FuInf", "FuNI", "MoInf", "MoNI")), y = IL6, fill = EthInf)) + 
  geom_violin() +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black") +
  geom_jitter(shape=16, position=position_jitter(0.2), size = 1.5)
pdf("CD14_IL6.pdf", width = 12, height = 16)
print(gene.plot)
dev.off()

# Form dubsets from the dataframe for each ethnic and infection group for t tests
df1_Inf <- df1 %>%
  filter(Infection == "Inf")
df1_NI <- df1 %>%
  filter(Infection == "NI")
df1_F <- df1 %>%
  filter(group_id == "Fulani")
df1_M <- df1 %>%
  filter(group_id == "Mossi")

leveneTest(df1_Inf$IL6, df1_Inf$group_id)
t.test(df1_Inf$IL6 ~ df1_Inf$group_id, var.eq = F, Paired = F)

leveneTest(df1_NI$IL6, df1_NI$group_id)
t.test(df1_NI$IL6 ~ df1_NI$group_id, var.eq = F, Paired = F)

leveneTest(df1_F$IL6, df1_F$Infection)
t.test(df1_F$IL6 ~ df1_F$Infection, var.eq = F, Paired = F)

leveneTest(df1_M$IL6, df1_M$Infection)
t.test(df1_M$IL6 ~ df1_M$Infection, var.eq = F, Paired = F)


#### Perform Gene Set Enrichment Analysis ####

library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(fgsea)
library(data.table)
library(ggplot2)

# Read the filtered DeSeq2 results, we are showing it for CD14 monocytes infected samples here
setwd("~/Documents/file_location/")
CD14_Mono <- read.csv("CD14 Mono_Sample-level_Deseq2_Filterd_padj_0.05_logFC_0.263.csv")

# Get "Hallmark" gene set database
H <- as.data.frame(msigdbr(species = "Homo sapiens", 
                           category = "H"))
H.genes.ls <- H %>% 
  #Keep gene ID that match expression data gene ID
  select(gs_name, gene_symbol) %>% 
  #Collapse all genes in each gene set into 1 row each
  group_by(gs_name) %>%
  summarise(all.genes = list(unique(gene_symbol))) %>%
  #Convert to list
  deframe()

# Subset gene names and log2 fold change columns into vector
FC_CD14_Mono <- CD14_Mono[,1:4]
FC_CD14_Mono <- FC_CD14_Mono[,-2]
FC_CD14_Mono <- FC_CD14_Mono[,-2]
FC_CD14_Mono.vec <- FC_CD14_Mono$logFC
names(FC_CD14_Mono.vec) <- FC_CD14_Mono$gene

# Set Score type to standard
scoreType <- "std"

# Run GSEA
gsea.H.CD14_Mono <- as.data.frame(fgseaSimple(pathways = H.genes.ls, 
                                              stats = FC_CD14_Mono.vec,
                                              nperm = 1000,
                                              scoreType = scoreType))

# Visualize significant GSEA
# Bar plot `NES` (normalized enrichment score) for the subset of significant GSEA at FDR < 0.1
gsea.H.CD14_Mono %>% 
  filter(padj <= 0.1) %>% 
  mutate(pathway = gsub("HALLMARK_","", pathway),
         pathway = gsub("_"," ", pathway)) %>% 
  
  ggplot(aes(x=reorder(pathway, NES), #Reorder gene sets by NES values
             y=NES)) +
  geom_col() +
  theme_classic() +
  lims(y=c(-3.2,3.2)) +
  coord_flip() +
  labs(y="Normalized enrichment score (NES)",
       x="Gene set",
       title = "Hallmark GSEA (FDR < 0.1)\nDown in Fulani <--         --> Up in Fulani")
write_csv(gsea.H.CD14_Mono, file = "GSEA-H_CD14.csv")


