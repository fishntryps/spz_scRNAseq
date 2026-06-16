
#### This script generates an Seurat RDS object from multiple 10X single-cell RNA-seq samples
#### it includes QC filtering, normalization, dimensional reduction, and RPCA-based integration

# Let's begin with a clean environment
rm(list=ls())

# Load required packages
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library("netdiffuseR")

# Set working directory
setwd("C:/Users/philipp/OneDrive - Harvard University/Desktop/single_cell_rnaseq")

# Read in raw 10X count matrices for each sample
R1_d11_WT_1BF_counts <- Read10X(data.dir = "mats_for_seurat/R1-d11-WT-1BF-SPZ_raw_feature_bc_matrix")
R1_d11_dsGFP_2BF_counts <- Read10X(data.dir = "mats_for_seurat/R1-d11-dsGFP-2BF-SPZ_raw_feature_bc_matrix")
R2_d11_WT_1BF_counts <- Read10X(data.dir = "mats_for_seurat/R2-d11-WT-1BF_raw_feature_bc_matrix")
R2_d11_dsGFP_2BF_counts <- Read10X(data.dir = "mats_for_seurat/R2-d11-dsGFP-2BF_raw_feature_bc_matrix")
R3_d11_WT_1BF_counts <- Read10X(data.dir = "mats_for_seurat/R3-d11-WT-1BF_raw_feature_bc_matrix")
R3_d11_dsGFP_2BF_counts <- Read10X(data.dir = "mats_for_seurat/R3-d11-dsGFP-2BF_raw_feature_bc_matrix")
R4_d11_WT_1BF_counts <- Read10X(data.dir = "mats_for_seurat/R4-d11-WT-1BF_raw_feature_bc_matrix")
R4_d11_WT_2BF_counts <- Read10X(data.dir = "mats_for_seurat/R4-d11-WT-2BF_raw_feature_bc_matrix")
R5_d11_WT_1BF_counts <- Read10X(data.dir = "mats_for_seurat/R5-d11-WT-1BF_raw_feature_bc_matrix")
R5_d11_WT_2BF_counts <- Read10X(data.dir = "mats_for_seurat/R5-d11-WT-2BF_raw_feature_bc_matrix")

# note that I refer to bioreps 5-6 as 4-5 in some analysis stages

# Use Seurat v5 assay structure
options(Seurat.object.assay.version = "v5")

# Create Seurat objects for each sample
R1_d11_WT_1BF <- CreateSeuratObject(counts = R1_d11_WT_1BF_counts,
                                    min.features = 100)

R1_d11_dsGFP_2BF <- CreateSeuratObject(counts = R1_d11_dsGFP_2BF_counts,
                                       min.features = 100)

R2_d11_WT_1BF <- CreateSeuratObject(counts = R2_d11_WT_1BF_counts,
                                    min.features = 100)

R2_d11_dsGFP_2BF <- CreateSeuratObject(counts = R2_d11_dsGFP_2BF_counts,
                                       min.features = 100)

R3_d11_WT_1BF <- CreateSeuratObject(counts = R3_d11_WT_1BF_counts,
                                    min.features = 100)

R3_d11_dsGFP_2BF <- CreateSeuratObject(counts = R3_d11_dsGFP_2BF_counts,
                                       min.features = 100)

R4_d11_WT_1BF <- CreateSeuratObject(counts = R4_d11_WT_1BF_counts,
                                    min.features = 100)

R4_d11_WT_2BF <- CreateSeuratObject(counts = R4_d11_WT_2BF_counts,
                                    min.features = 100)

R5_d11_WT_1BF <- CreateSeuratObject(counts = R5_d11_WT_1BF_counts,
                                    min.features = 100)

R5_d11_WT_2BF <- CreateSeuratObject(counts = R5_d11_WT_2BF_counts,
                                    min.features = 100)

# Inspect metadata
head(R1_d11_WT_1BF@meta.data)
head(R1_d11_dsGFP_2BF@meta.data)

# Merge all samples into one Seurat object
merged_seurat <- merge(
  x = R1_d11_WT_1BF,
  y = c(
    R1_d11_dsGFP_2BF,
    R2_d11_WT_1BF,
    R2_d11_dsGFP_2BF,
    R3_d11_WT_1BF,
    R3_d11_dsGFP_2BF,
    R4_d11_WT_1BF,
    R4_d11_WT_2BF,
    R5_d11_WT_1BF,
    R5_d11_WT_2BF),
  add.cell.id = c(
    "R1_d11_WT_1BF",
    "R1_d11_dsGFP_2BF",
    "R2_d11_WT_1BF",
    "R2_d11_dsGFP_2BF",
    "R3_d11_WT_1BF",
    "R3_d11_dsGFP_2BF",
    "R4_d11_WT_1BF",
    "R4_d11_WT_2BF",
    "R5_d11_WT_1BF",
    "R5_d11_WT_2BF"))

# Inspect merged metadata
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

# Load external metadata table
xx <- read.table("C:/Users/philipp/Downloads/mcameta.txt", header=T)

# Store cell names in metadata
merged_seurat@meta.data$cells <- rownames(merged_seurat@meta.data)

# Add developmental stage metadata
merged_seurat@meta.data$stage <- xx$stage[
  match(merged_seurat@meta.data$cells,
        paste("mca_", xx$cell, sep=""))]

# Label current experiment cells
merged_seurat@meta.data$stage[
  grepl("^R", merged_seurat@meta.data$cells)] <- "sgSpz"

# Add day metadata
merged_seurat@meta.data$day <- xx$day[
  match(merged_seurat@meta.data$cells,
        paste("mca_", xx$cell, sep=""))]

# Override day labels for d11 samples
merged_seurat@meta.data$day[
  grepl("_d11_", merged_seurat@meta.data$cells)] <- "day11"

# Combine day and stage into one label
merged_seurat@meta.data$day_stage <- paste(
  merged_seurat@meta.data$day,
  merged_seurat@meta.data$stage,
  sep="_")

################################################
### QUALITY CONTROL ############################
################################################

# Calculate genes per UMI metric
merged_seurat$log10GenesPerUMI <-
  log10(merged_seurat$nFeature_RNA) /
  log10(merged_seurat$nCount_RNA)

# Calculate mitochondrial transcript percentage
merged_seurat$mitoRatio <- PercentageFeatureSet(
  object = merged_seurat,
  pattern = "^PF3D7-MIT"
)

merged_seurat$mitoRatio <-
  merged_seurat@meta.data$mitoRatio / 100

# Calculate apicoplast transcript percentage
merged_seurat$apicoRatio <- PercentageFeatureSet(
  object = merged_seurat,
  pattern = "^PF3D7-API")

merged_seurat$apicoRatio <-
  merged_seurat@meta.data$apicoRatio / 100

# Display object summary
merged_seurat

# Copy metadata into dataframe
metadata <- merged_seurat@meta.data

# Add cell names
metadata$cells <- rownames(metadata)

# Create sample labels
metadata$sample <- NA

metadata$sample[
  which(str_detect(metadata$cells, "^R1_d11_WT_1BF"))] <- "R1_d11_WT_1BF"

metadata$sample[
  which(str_detect(metadata$cells, "^R1_d11_dsGFP_2BF"))] <- "R1_d11_dsGFP_2BF"

metadata$sample[
  which(str_detect(metadata$cells, "^R2_d11_WT_1BF"))] <- "R2_d11_WT_1BF"

metadata$sample[
  which(str_detect(metadata$cells, "^R2_d11_dsGFP_2BF"))] <- "R2_d11_dsGFP_2BF"

metadata$sample[
  which(str_detect(metadata$cells, "^R3_d11_WT_1BF"))] <- "R3_d11_WT_1BF"

metadata$sample[
  which(str_detect(metadata$cells, "^R3_d11_dsGFP_2BF"))] <- "R3_d11_dsGFP_2BF"

metadata$sample[
  which(str_detect(metadata$cells, "^R4_d11_WT_1BF"))] <- "R4_d11_WT_1BF"

metadata$sample[
  which(str_detect(metadata$cells, "^R4_d11_WT_2BF"))] <- "R4_d11_WT_2BF"

metadata$sample[
  which(str_detect(metadata$cells, "^R5_d11_WT_1BF"))] <- "R5_d11_WT_1BF"

metadata$sample[
  which(str_detect(metadata$cells, "^R5_d11_WT_2BF"))] <- "R5_d11_WT_2BF"

# Create feeding-condition labels
metadata$num_feeds <- NA
metadata$num_feeds[
  which(str_detect(metadata$cells, "1BF"))] <- "1BF"

metadata$num_feeds[
  which(str_detect(metadata$cells, "2BF"))] <- "2BF"

# Add experiment metadata
metadata$dpbm <- NA
metadata$dpbm[
  which(str_detect(metadata$cells, "d11"))] <- "d11"

metadata$inject <- NA
metadata$inject[
  which(str_detect(metadata$cells, "dsGFP"))] <- "dsGFP"

metadata$inject[
  which(str_detect(metadata$cells, "WT"))] <- "WT"

# Rename metadata columns
metadata <- metadata %>%
  dplyr::rename(
    seq_folder = orig.ident,
    nUMI = nCount_RNA,
    nGene = nFeature_RNA)

# Store metadata back into Seurat object
merged_seurat@meta.data <- metadata

# Plot number of cells per sample
metadata %>%
  ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +
  theme_classic()

# Plot UMI distribution
metadata %>%
  ggplot(aes(color=sample, x=nUMI, fill= sample)) +
  geom_density(alpha = 0.2)

# Plot gene-count distribution
metadata %>%
  ggplot(aes(color=sample, x=nGene, fill= sample)) +
  geom_density(alpha = 0.2)

# Plot genes-per-UMI complexity metric
metadata %>%
  ggplot(aes(x=log10GenesPerUMI,
             color = sample,
             fill=sample)) +
  geom_density(alpha = 0.2)

# Plot mitochondrial content distribution
metadata %>%
  ggplot(aes(color=sample,
             x=mitoRatio,
             fill=sample)) +
  geom_density(alpha = 0.2)

# Plot gene vs UMI relationship
metadata %>%
  ggplot(aes(x=nUMI,
             y=nGene,
             color=mitoRatio)) +
  geom_point()

# Filter low-quality cells
filtered_seurat_downsampled <- subset(
  x = merged_seurat,
  subset =
    (nUMI >= 400) &
    (nGene >= 125) &
    (log10GenesPerUMI > 0.80) &
    (mitoRatio < 0.20))

# Join assay layers
obj_new = JoinLayers(filtered_seurat_downsampled)

# Extract count matrix
counts <- GetAssayData(
  object = obj_new,
  layer = "counts")

# Identify genes expressed in at least 10 cells
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Subset count matrix
filtered_counts <- counts[keep_genes, ]

# Check dimensions before and after filtering
dim(counts)
dim(filtered_counts)

# Rebuild Seurat object with filtered genes
filtered_seurat_downsampled <- CreateSeuratObject(
  filtered_counts,
  meta.data = filtered_seurat_downsampled@meta.data)

# Save filtered metadata
metadata_clean <- filtered_seurat_downsampled@meta.data

# Determine minimum cell count per replicate pair
R1_min <- min(
  length(metadata_clean$sample[
    metadata_clean$sample=="R1_d11_WT_1BF"
  ]),
  length(metadata_clean$sample[
    metadata_clean$sample=="R1_d11_dsGFP_2BF"
  ]))

# Set cell identities
filtered_seurat_downsampled <- SetIdent(
  filtered_seurat_downsampled,
  value = "sample")

# Downsample cells evenly across conditions
R1.cell.list <- WhichCells(
  filtered_seurat_downsampled,
  idents = c("R1_d11_WT_1BF",
             "R1_d11_dsGFP_2BF"),
  downsample = R1_min)

# Store cleaned metadata
metadata_clean_downsampled <-
  filtered_seurat_downsampled@meta.data

# Load plotting helper package
library('gcookbook')

# Plot QC metrics across samples
metadata_clean_downsampled %>%
  ggplot(aes(x=sample2, fill=num_feeds)) +
  geom_bar()

################################################
### NORMALIZATION ##############################
################################################

# Load packages again
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Normalize expression values
seurat_phase <- NormalizeData(filtered_seurat_downsampled)

# Inspect normalized RNA assay
seurat_phase[["RNA"]]$counts

# Identify variable genes
seurat_phase <- FindVariableFeatures(
  seurat_phase,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE)

# Scale expression values
seurat_phase <- ScaleData(seurat_phase)

# Run PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot PCA split by sample
DimPlot(seurat_phase, split.by = "orig.ident")

# Plot PCA variance explained
ElbowPlot(seurat_phase)

# Run UMAP
early_umap <- RunUMAP(
  seurat_phase,
  dims = 1:5,
  reduction = "pca"
)

# Plot UMAP colored by feeding condition
DimPlot(
  early_umap,
  split.by = 'orig.ident',
  group.by = 'num_feeds',
  shape.by = 'num_feeds',
  pt.size = 1.6,
  reduction = "umap"
)

# Plot QC features on UMAP
FeaturePlot(early_umap, features = "nGene")
FeaturePlot(early_umap, features = "mitoRatio")
FeaturePlot(early_umap, features = "PF3D7-1032700")

# Summarize mitochondrial ratio distribution
summary(seurat_phase@meta.data$mitoRatio)

# Categorize cells by mitochondrial content
seurat_phase@meta.data$mitoFr <- cut(
  seurat_phase@meta.data$mitoRatio,
  breaks = c(
    -Inf,
    summary(seurat_phase@meta.data$mitoRatio)[3],
    summary(seurat_phase@meta.data$mitoRatio)[5],
    Inf
  ),
  labels = c("Low", "Elevated", "High"))

# Split object by sample
split_seurat <- SplitObject(
  seurat_phase,
  split.by = "sample")

# Keep selected samples
split_seurat <- split_seurat[c(
  "R1_d11_WT_1BF",
  "R1_d11_dsGFP_2BF",
  "R2_d11_WT_1BF",
  "R2_d11_dsGFP_2BF",
  "R3_d11_WT_1BF",
  "R3_d11_dsGFP_2BF",
  "R4_d11_WT_1BF",
  "R4_d11_WT_2BF",
  "R5_d11_WT_1BF",
  "R5_d11_WT_2BF")]

# # Increase memory limit for SCTransform
# options(future.globals.maxSize = 4000 * 1024^2)
# 
# # Run SCTransform on each sample
# for (i in 1:length(split_seurat)) {
#   split_seurat[[i]] <- SCTransform(split_seurat[[i]])
# }
# 
# # Find variable genes in SCT assay
# for (i in 1:length(split_seurat)) {
#   split_seurat[[i]] <- FindVariableFeatures(
#     split_seurat[[i]],
#     assay="SCT"
#   )
# }

# Merge samples
obj <- merge(
  x=split_seurat[[1]],
  y=c(
    split_seurat[[2]],
    split_seurat[[3]],
    split_seurat[[4]],
    split_seurat[[5]],
    split_seurat[[6]],
    split_seurat[[7]],
    split_seurat[[8]],
    split_seurat[[9]],
    split_seurat[[10]]
  )
)

# Switch back to RNA assay
DefaultAssay(obj) <- "RNA"

# Normalize merged object
obj <- NormalizeData(obj)

# Find variable features
obj <- FindVariableFeatures(obj)

# Scale data
obj <- ScaleData(obj)

# Run PCA
obj <- RunPCA(obj)

# Run RPCA integration
obj <- IntegrateLayers(
  object = obj,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = TRUE
)

obj <- RunUMAP(obj, dims = 1:9, reduction = "integrated.rpca", reduction.name="integrated.rpca.umap")

# Save integrated object
saveRDS(obj, "results/5reps_RPCA_9dimUMAP_400umi_125gene_only10x.rds")