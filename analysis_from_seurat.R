#### This script operates on an RDS object saved from a Seurat pipeline
#### It also uses 2 annotation files (annotations.txt and go_terms.txt)
#### And it calls on some functions saved in 10x_functions_github.R

#### you might want to clear your workspace before you start

rm(list=ls())

### Set your working directory to your directory of choice and make sure these files are present 

setwd("C:/Users/philipp/OneDrive - Harvard University/Desktop/single_cell_rnaseq")

#### load softwares
#### you will first have needed to install these to your computer
#### some can be installed from CRAN via install.packages("packagename")
#### otherwise BiocManager::install("packagename") is super useful

library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)
library(scater)
library(scran)
library(batchelor)
library(edgeR)
library(bluster)
library(grid)
library(ggplot2)
library(slingshot)
library(monocle3)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(magrittr)
library(corrplot)
library(gplots)
library(rstatix)
library(ggforce)
library(data.table)
library(slingshot)

#### here are some notes for some of the trickier installations above, e.g. monocle3
#install.packages("devtools")
#library(devtools)
#install.packages('rlang')
#devtools::install_github('cole-trapnell-lab/monocle3')
#devtools::install_github('satijalab/seurat-data')
#remotes::install_github('satijalab/seurat-wrappers') ### had to set API for this

### first let's load R functions used in this script

source("10x_functions_github.R")

### here also some colors we use later on

yan_colors <- c(
'#a6cee3', ## Ld2
'#1f78b4', ## pseudo midpoint
'#b2df8a', ## Ld3
'#33a02c', ## pseudo highpoint
'#fb9a99', ## 1BF
'#e31a1c',## heat highpoint
'#fdbf6f', ## Ld1
'#ff7f00',
'#cab2d6', ## 2BF
'#ffff99', ## heat lowpoint and pseudo lowpoint
'#6a3d9a')

colLd1 <- '#fdbf6f'
colLd2 <- '#a6cee3'
colLd3 <- '#b2df8a'

col1BF <- '#fb9a99'
col2BF <- '#cab2d6'

colHeatLow <- '#ffff99'
colHeatHigh <- '#e31a1c'

colSling <- c('#ffff99', '#1f78b4', '#33a02c')  
colSlingAlt1 <- c('#ff7f00','#e31a1c','#6a3d9a')
colSlingOld <- c('turquoise', 'blue', 'hotpink')

### here we read in the Seurat object containing all the consolidated scRNAseq info from 5 bioreps

seurat_integrated <- readRDS("results/5reps_RPCA_9dimUMAP_400umi_125gene_only10x_with_clusters_and_pseutodimes.rds")

### and gene annotations...

annotations <- read.table("annotations.txt", header=T)

seurat_integrated@meta.data$repset <- NA
seurat_integrated@meta.data$repset[seurat_integrated@meta.data$orig.ident %in% c("R1","R2","R3")] <- "R123"
seurat_integrated@meta.data$repset[seurat_integrated@meta.data$orig.ident %in% c("R4","R5")] <- "R45" ### note that I refer to bioreps 5-6 as 4-5


seurat_integrated <- JoinLayers(seurat_integrated)
seurat_integrated[["integrated.rpca.umap"]]@cell.embeddings[,1] <- -1*seurat_integrated[["integrated.rpca.umap"]]@cell.embeddings[,1]
seurat_integrated[["pca"]]@cell.embeddings[,1] <- -1*seurat_integrated[["pca"]]@cell.embeddings[,1]
seurat_integrated[["umap"]] <- seurat_integrated[["integrated.rpca.umap"]]

cds <- as.cell_data_set(seurat_integrated)
cds <- cluster_cells(cds, resolution=0.00015)
plot_cells(cds, color_cells_by = "cluster", cell_size = 1.1, show_trajectory_graph = FALSE, group_label_size=10, reduction_method = "UMAP")
rowData(cds)$gene_short_name <- row.names(rowData(cds))
plot_cells(cds,genes =  c("PF3D7-0629300","PF3D7-1147800"),cell_size = 2)
plot_cells(cds, scale_to_range=T, genes =  c(
  "PF3D7-1133700", "PF3D7-0718300", "PF3D7-0404800",
  "PF3D7-1431500", "PF3D7-0812300", "PF3D7-0818600",
  "PF3D7-0202100", "PF3D7-0206900", "PF3D7-0902900",
  "PF3D7-0922200", "PF3D7-0829500"),
  cell_size = 2.1, show_trajectory_graph = F, label_cell_groups = F) + 
            xlim(-9.4,6.6) + ylim(-8,8) +
             #theme(legend.key.height=unit(2, "cm"), legend.text = element_text(size=12)) +
             coord_fixed() + theme(aspect.ratio=1) +
  scale_color_gradientn(colors=colorRampPalette(c("grey","cornsilk", "mistyrose", "lightpink", "maroon", "maroon4", "#330033", "black"))(100)) 

seurat_integrated$leiden_res00015_clusters <- as.character(cds@clusters$UMAP$clusters)
seurat_integrated$leiden_res00015_clusters[seurat_integrated$leiden_res00015_clusters %in% "3"] <-"Ld1" 
seurat_integrated$leiden_res00015_clusters[seurat_integrated$leiden_res00015_clusters %in% "1"] <-"Ld2" 
seurat_integrated$leiden_res00015_clusters[seurat_integrated$leiden_res00015_clusters %in% "2"] <-"Ld3" 
seurat_integrated$leiden_res00015_clusters <- as.factor(seurat_integrated$leiden_res00015_clusters)

seurat_integrated$leiden_res00015_mergedclusters <- NA
seurat_integrated$leiden_res00015_mergedclusters[seurat_integrated$leiden_res00015_clusters %in% c("Ld2", "Ld3")] <-"Ld2+3" 
seurat_integrated$leiden_res00015_mergedclusters[seurat_integrated$leiden_res00015_clusters %in% "Ld1"] <-"Ld1"
seurat_integrated$leiden_res00015_mergedclusters <- as.factor(seurat_integrated$leiden_res00015_mergedclusters)

# Plot cell clusters
DimPlot(seurat_integrated, group.by="leiden_res00015_clusters",
        reduction = "umap",
        label = TRUE,
        label.size = 10,
        pt.size=1.5) + scale_color_manual(values=c("Ld1"=colLd1, "Ld2"=colLd2, "Ld3"=colLd3))

# Plot cell clusters
DimPlot(seurat_integrated, group.by="leiden_res00015_clusters",
        reduction = "integrated.jointpca.umap",
        label = TRUE,
        label.size = 10,
        pt.size=1.5)

#### Find markers for each cluster number using FindConservedMarkers function
#### These would be cluster marker signals conserved across all samples
#### One can also use more lenient criteria, e.g. changing groupvar to "orig.ident" such that conservation must occur across bioreps instead of samples


# Assign identity of clusters
Idents(object = seurat_integrated) <- "leiden_res00015_clusters"

clus1_cons_mark <- get_conserved(clus_num="Ld1", groupvar="sample")
clus2_cons_mark <- get_conserved(clus_num="Ld2", groupvar="sample")
clus3_cons_mark <- get_conserved(clus_num="Ld3", groupvar="sample")

all_cons_mark <- rbind(clus1_cons_mark,clus2_cons_mark, clus3_cons_mark)

### add column representing average log2FC for each marker gene  

all_cons_mark$avg_log2FC_btw_samp <- (all_cons_mark$R1_d11_WT_1BF_avg_log2FC + all_cons_mark$R1_d11_dsGFP_2BF_avg_log2FC + all_cons_mark$R2_d11_WT_1BF_avg_log2FC + all_cons_mark$R2_d11_dsGFP_2BF_avg_log2FC + all_cons_mark$R3_d11_WT_1BF_avg_log2FC + all_cons_mark$R3_d11_dsGFP_2BF_avg_log2FC + all_cons_mark$R4_d11_WT_1BF_avg_log2FC + all_cons_mark$R4_d11_WT_2BF_avg_log2FC + all_cons_mark$R5_d11_WT_1BF_avg_log2FC + all_cons_mark$R5_d11_WT_2BF_avg_log2FC)/10

### subselect markers with...
max_pval_to_use=0.05
log2FC_to_use=0.584962501 ### log2FC = log2(1.5) = 0.584962501

### Show results for columns of interest
### colnames(all_cons_mark)

filtered_cons_mark <- all_cons_mark %>% filter(max_pval<max_pval_to_use & avg_log2FC_btw_samp>log2FC_to_use) %>% 
  select("cluster_id", "gene", "description", "max_pval", "avg_log2FC_btw_samp")

#### Have a look at the marker results

filtered_cons_mark[filtered_cons_mark$cluster_id=="Ld1",] %>% arrange(-avg_log2FC_btw_samp)  %>% filter(max_pval<0.00001)

### And plot marker expression of interest in UMAP
### You can type feature(s) in manually (make sure to use hyphens in place of underscores) or use a subsection of the filtered_cons_mark table

# cluster_id          gene                                                description     max_pval avg_log2FC_btw_samp
# 1         Ld1 PF3D7-1351800           conserved_Plasmodium_protein%2C_unknown_function 2.539931e-26           4.9985844
# 2         Ld1 PF3D7-0315900           conserved_Plasmodium_protein%2C_unknown_function 8.523644e-08           4.3668274
# 3         Ld1 PF3D7-0404800                sporozoite-specific_protein_S10%2C_putative 1.727134e-19           2.7752727
# 4         Ld1 PF3D7-0820300           conserved_Plasmodium_protein%2C_unknown_function 2.841974e-06           2.6082287
# 5         Ld1 PF3D7-0620200           conserved_Plasmodium_protein%2C_unknown_function 6.221067e-09           2.5636488
# 6         Ld1 PF3D7-0315800                            zinc_finger_protein%2C_putative 1.012367e-18           2.5343904
# 7         Ld1 PF3D7-0822700          TRP1_thrombospondin-related_protein_1%2C_putative 8.698480e-08           2.3336369
# 8         Ld1 PF3D7-1236500           conserved_Plasmodium_protein%2C_unknown_function 1.081041e-06           2.2592755
# 9         Ld1 PF3D7-0718300                    CRMP2_cysteine_repeat_modular_protein_2 8.996846e-06           2.2459309
# 10        Ld1 PF3D7-0720400                  AIF_apoptosis-inducing_factor%2C_putative 3.107140e-06           2.2297328
# 11        Ld1 PF3D7-1442600          TREP_sporozoite-specific_transmembrane_protein_S6 1.028505e-07           2.2182937
# 12        Ld1 PF3D7-0323800           conserved_Plasmodium_protein%2C_unknown_function 5.762631e-06           2.2037280
# 13        Ld1 PF3D7-0802000                   GDH3_glutamate_dehydrogenase%2C_putative 4.719360e-06           1.9526926
# 14        Ld1 PF3D7-1147800 MAEBL_membrane_associated_erythrocyte_binding-like_protein 1.350129e-09           1.8915299
# 15        Ld1 PF3D7-1008600                zinc_finger_(CCCH_type)_protein%2C_putative 6.561558e-08           1.8692827
# 16        Ld1 PF3D7-1026500           conserved_Plasmodium_protein%2C_unknown_function 3.047842e-06           0.9335927


filtered_cons_mark[filtered_cons_mark$cluster_id=="Ld2",] %>% arrange(-avg_log2FC_btw_samp)  %>% filter(max_pval<0.00001)

# cluster_id          gene                                                              description     max_pval avg_log2FC_btw_samp
# 1         Ld2 PF3D7-0913800                         conserved_Plasmodium_protein%2C_unknown_function 3.809682e-21           1.3255403
# 2         Ld2 PF3D7-1431500                                 MAPK1_mitogen-activated_protein_kinase_1 2.252258e-06           1.2582663
# 3         Ld2 PF3D7-0919500 MFS3_major_facilitator_superfamily_domain-containing_protein%2C_putative 3.438402e-11           1.2437179
# 4         Ld2 PF3D7-1016000                         conserved_Plasmodium_protein%2C_unknown_function 1.806142e-07           1.1792336
# 5         Ld2 PF3D7-0931700                               PIH1_domain-containing_protein%2C_putative 6.282863e-06           1.1354795
# 6         Ld2 PF3D7-0812300                            SSP3_sporozoite_surface_protein_3%2C_putative 1.019750e-22           1.1336293
# 7         Ld2 PF3D7-0814600                                                SIMP_concavin%2C_putative 1.347489e-06           1.0890687
# 8         Ld2 PF3D7-0931600                                    conserved_protein%2C_unknown_function 2.002794e-08           1.0868281
# 9         Ld2 PF3D7-0818600                                      PBLP_BEM46-like_protein%2C_putative 9.887292e-08           1.0476806
# 10        Ld2 PF3D7-1207400                         conserved_Plasmodium_protein%2C_unknown_function 9.508154e-19           0.9182901
# 11        Ld2 PF3D7-1335900                                        TRAP_sporozoite_surface_protein_2 1.627745e-09           0.5878101
# 12        Ld2 PF3D7-1011500                                    conserved_protein%2C_unknown_function 1.860444e-09           0.5873966

filtered_cons_mark[filtered_cons_mark$cluster_id=="Ld3",] %>% arrange(-avg_log2FC_btw_samp) %>% filter(max_pval<0.00001)
# cluster_id          gene                                                    description     max_pval avg_log2FC_btw_samp
# 1         Ld3 PF3D7-0202100 LSAP2_Plasmodium_exported_protein_(PHISTc)%2C_unknown_function 5.445905e-08           4.8898704
# 2         Ld3 PF3D7-1436300                           PTEX150_translocon_component_PTEX150 8.889446e-07           3.4859510
# 3         Ld3 PF3D7-1311600               conserved_Plasmodium_protein%2C_unknown_function 4.179444e-06           2.8062445
# 4         Ld3 PF3D7-0306100               conserved_Plasmodium_protein%2C_unknown_function 1.212233e-07           2.6846505
# 5         Ld3 PF3D7-1333100               conserved_Plasmodium_protein%2C_unknown_function 5.523517e-08           2.6505816
# 6         Ld3 PF3D7-1231200               conserved_Plasmodium_protein%2C_unknown_function 6.764627e-08           1.7600271
# 7         Ld3 PF3D7-1142900               conserved_Plasmodium_protein%2C_unknown_function 2.080469e-10           1.7379623
# 8         Ld3 PF3D7-0911100                    START_domain-containing_protein%2C_putative 2.143564e-06           1.6124914
# 9         Ld3 PF3D7-0616500                                          TLP_TRAP-like_protein 5.077124e-14           1.3832566
# 10        Ld3 PF3D7-1016900             ETRAMP10.3_early_transcribed_membrane_protein_10.3 1.589105e-06           1.2706308
# 11        Ld3 PF3D7-1446900                glutaminyl-peptide_cyclotransferase%2C_putative 6.118370e-06           0.8843082


markers.to.plot <- 
  
  (c(filtered_cons_mark[filtered_cons_mark$cluster_id=="Ld1",] %>% arrange(-avg_log2FC_btw_samp) %>% filter(max_pval<0.001 & !grepl("conserved.*protein%2C_unknown_function|^NA$", description)) %>% pull(gene) %>% head(6),
     filtered_cons_mark[filtered_cons_mark$cluster_id=="Ld2",] %>% arrange(-avg_log2FC_btw_samp) %>% filter(max_pval<0.001 & !grepl("conserved.*protein%2C_unknown_function|^NA$", description)) %>% pull(gene) %>% head(6),
     filtered_cons_mark[filtered_cons_mark$cluster_id=="Ld3",] %>% arrange(-avg_log2FC_btw_samp) %>% filter(max_pval<0.001 & !grepl("conserved.*protein%2C_unknown_function|^NA$", description)) %>% pull(gene) %>% head(6)))

tempseurat <- SetIdent(seurat_integrated, value = "leiden_res00015_clusters")

mylabs <- as.data.frame(markers.to.plot) %>% 
  left_join(annotations[, c("gene_name", "description")] %>% 
              mutate(mylab=paste(str_extract(description, 
                                             "^............."), "...", gene_name)), 
            by = c("markers.to.plot"="gene_name")) %>% 
  pull(mylab)

DotPlot(tempseurat, features = markers.to.plot, scale = F,
        cols=c("blue","red"), dot.scale = 15) + scale_x_discrete(labels=mylabs) +
  theme(axis.text.x = element_text(angle=34, hjust=1, vjust=1)) +
  theme(plot.margin = margin(2,2,2,4, "cm")) + labs(y="Cluster\n", x="Gene") +
  theme(axis.text.x=element_text(size=15, 
                                 color=(c(rep("#729ECE",6),rep("#FF9E4A",6),rep("#67BF5C",6))))) +
  theme(axis.text.y=element_text(size=15, 
                                 color=c("#729ECE","#FF9E4A","#67BF5C")))

markers.to.plot <- 
  
  rev(c(filtered_cons_mark[filtered_cons_mark$cluster_id=="Ld1",] %>% arrange(-avg_log2FC_btw_samp) %>% filter(max_pval<0.001 & !grepl("conserved.*protein%2C_unknown_function|^NA$", description)) %>% pull(gene) %>% head(6),
        filtered_cons_mark[filtered_cons_mark$cluster_id=="Ld2",] %>% arrange(-avg_log2FC_btw_samp) %>% filter(max_pval<0.001 & !grepl("conserved.*protein%2C_unknown_function|^NA$", description)) %>% pull(gene) %>% head(6),
        filtered_cons_mark[filtered_cons_mark$cluster_id=="Ld3",] %>% arrange(-avg_log2FC_btw_samp) %>% filter(max_pval<0.001 & !grepl("conserved.*protein%2C_unknown_function|^NA$", description)) %>% pull(gene) %>% head(6)))

tempseurat <- SetIdent(seurat_integrated, value = "leiden_res00015_clusters")

mylabs <- as.data.frame(markers.to.plot) %>% 
  left_join(annotations[, c("gene_name", "description")] %>% 
              mutate(mylab=paste(str_extract(description, 
                                             "^............."), "...", gene_name)), 
            by = c("markers.to.plot"="gene_name")) %>% 
  pull(mylab)

DotPlot(tempseurat, features = markers.to.plot, scale = F,
        cols=c("blue","red"), dot.scale = 8) + scale_x_discrete(labels=mylabs) + 
  coord_flip() +
  #theme(plot.margin = margin(1,2,1,1, "cm")) + 
  labs(y="\nCluster", x="Gene") +
  theme(axis.text.y=element_text(size=15, 
                                 color=rev(c(rep("#729ECE",6),rep("#FF9E4A",6),rep("#67BF5C",6))))) +
  theme(axis.text.x=element_text(size=15, 
                                 color=c("#729ECE","#FF9E4A","#67BF5C")))

# Plot gene expression on UMAP

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("PF3D7-1032700"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE, pt.size = 2.2)


FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = filtered_cons_mark$gene[filtered_cons_mark$cluster_id=="Ld1" & filtered_cons_mark$avg_log2FC_btw_samp >2.5][1:4], 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE, pt.size = 2.2)

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = filtered_cons_mark$gene[filtered_cons_mark$cluster_id=="Ld2" & filtered_cons_mark$avg_log2FC_btw_samp >1][1:4], 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE, pt.size = 2.2)

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = filtered_cons_mark$gene[filtered_cons_mark$cluster_id=="Ld3" & filtered_cons_mark$avg_log2FC_btw_samp >3][1:4], 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE, pt.size = 2.2)

annotations[annotations$gene_name %in% filtered_cons_mark$gene[filtered_cons_mark$cluster_id=="Ld3" & filtered_cons_mark$avg_log2FC_btw_samp >2],] 

day11_1BF_sgSpz_x <- Embeddings(object = seurat_integrated[["umap"]])[rownames(Embeddings(object = seurat_integrated[["umap"]])) %in% seurat_integrated@meta.data$cells[
  grepl("d11_.*1BF", seurat_integrated@meta.data$sample)],1]
day11_1BF_sgSpz_y <- Embeddings(object = seurat_integrated[["umap"]])[rownames(Embeddings(object = seurat_integrated[["umap"]])) %in% seurat_integrated@meta.data$cells[
  grepl("d11_.*1BF", seurat_integrated@meta.data$sample)],2]
day11_1BF_sgSpz <- as.data.frame(cbind(day11_1BF_sgSpz_x, day11_1BF_sgSpz_y))
colnames(day11_1BF_sgSpz) = c("x","y")

day11_2BF_sgSpz_x <- Embeddings(object = seurat_integrated[["umap"]])[rownames(Embeddings(object = seurat_integrated[["umap"]])) %in% seurat_integrated@meta.data$cells[
  grepl("d11_.*2BF", seurat_integrated@meta.data$sample)],1]
day11_2BF_sgSpz_y <- Embeddings(object = seurat_integrated[["umap"]])[rownames(Embeddings(object = seurat_integrated[["umap"]])) %in% seurat_integrated@meta.data$cells[
  grepl("d11_.*2BF", seurat_integrated@meta.data$sample)],2]
day11_2BF_sgSpz <- as.data.frame(cbind(day11_2BF_sgSpz_x, day11_2BF_sgSpz_y))
colnames(day11_2BF_sgSpz) = c("x","y")

r45_day11_1BF_sgSpz_x <- Embeddings(object = seurat_integrated[["umap"]])[rownames(Embeddings(object = seurat_integrated[["umap"]])) %in% seurat_integrated@meta.data$cells[
  grepl("R[4-5]_d11_.*1BF", seurat_integrated@meta.data$sample)],1]
r45_day11_1BF_sgSpz_y <- Embeddings(object = seurat_integrated[["umap"]])[rownames(Embeddings(object = seurat_integrated[["umap"]])) %in% seurat_integrated@meta.data$cells[
  grepl("R[4-5]_d11_.*1BF", seurat_integrated@meta.data$sample)],2]
r45_day11_1BF_sgSpz <- as.data.frame(cbind(r45_day11_1BF_sgSpz_x, r45_day11_1BF_sgSpz_y))
colnames(r45_day11_1BF_sgSpz) = c("x","y")

r45_day11_2BF_sgSpz_x <- Embeddings(object = seurat_integrated[["umap"]])[rownames(Embeddings(object = seurat_integrated[["umap"]])) %in% seurat_integrated@meta.data$cells[
  grepl("R[4-5]_d11_.*2BF", seurat_integrated@meta.data$sample)],1]
r45_day11_2BF_sgSpz_y <- Embeddings(object = seurat_integrated[["umap"]])[rownames(Embeddings(object = seurat_integrated[["umap"]])) %in% seurat_integrated@meta.data$cells[
  grepl("R[4-5]_d11_.*2BF", seurat_integrated@meta.data$sample)],2]
r45_day11_2BF_sgSpz <- as.data.frame(cbind(r45_day11_2BF_sgSpz_x, r45_day11_2BF_sgSpz_y))
colnames(r45_day11_2BF_sgSpz) = c("x","y")

sce.sling <- slingshot(as.SingleCellExperiment(seurat_integrated), reducedDim='UMAP')
head(sce.sling$slingPseudotime_1)
embedded <- embedCurves(sce.sling, "UMAP")
embedded <- slingCurves(embedded)[[1]] # only 1 path.
embedded <- data.frame(embedded$s[embedded$ord,])
sling_df <- cbind(as.data.frame(sce.sling$cells), as.data.frame(sce.sling$slingPseudotime_1))
seurat_integrated@meta.data$slingPseudotime_1 <- sling_df$`sce.sling$slingPseudotime_1`[match(seurat_integrated@meta.data$cells, sling_df$`sce.sling$cells`)]

p1 <- ggplotGrob(plotPCA(sce.sling, colour_by="leiden_res00015_clusters", point_size=1.25) +
                          theme(legend.position = "none") + 
                   #xlim(-9.4,6.6) + ylim(-8,8) + 
                   scale_color_manual(values=c("Ld1"=colLd1, "Ld2"=colLd2, "Ld3"=colLd3)) +
                          theme(text = element_text(size = 7.5),
                                              axis.title = element_text(size = 7.5),
                                              axis.text = element_text(size = 7.5)))
p2 <- ggplotGrob(plotUMAP(sce.sling, colour_by="leiden_res00015_clusters", point_size=1.25) +
                          theme(legend.position = "none") + xlim(-9.4,6.6) + ylim(-8,8) + coord_fixed() + theme(text = element_text(size = 7.5),
                                                                                  axis.title = element_text(size = 7.5),
                                                                                  axis.text = element_text(size = 7.5)) + theme(aspect.ratio=1) +
                          scale_color_manual(values=c("Ld1"=colLd1, "Ld2"=colLd2, "Ld3"=colLd3)))
p3 <- ggplotGrob(plotUMAP(sce.sling, colour_by="slingPseudotime_1", point_size=1.25) +
                          geom_path(data=embedded, aes(x=umap_1, y=umap_2), size=0.8) + 
                          geom_text(label="\u25BA", data=tail(head(embedded,2000),1),
                                    mapping=aes(x=umap_1, y=umap_2), size=4, angle=6) +
                          geom_text(label="\u25BA", data=tail(head(embedded,4060),1),
                                    mapping=aes(x=umap_1, y=umap_2), size=4, angle=-60) +
                   theme(text = element_text(size = 7.5),
                   axis.title = element_text(size = 7.5),
                   axis.text = element_text(size = 7.5)) +
                          theme(legend.position = "none") + 
                   xlim(-9.4,6.6) + ylim(-8,8) +
                          #theme(legend.key.height=unit(2, "cm"), legend.text = element_text(size=12)) +
                          coord_fixed() + theme(aspect.ratio=1) +
                          #scale_color_gradient2(low = "turquoise", mid = "blue", high = "hotpink", midpoint = 11,
                          scale_color_gradientn(colors=colorRampPalette(colSling)(10),
                                                breaks=c(seq(min(sce.sling$slingPseudotime_1),
                                                             max(sce.sling$slingPseudotime_1, length.out=10)))))

# Plot cell clusters
p4 <- ggplotGrob(DimPlot(subset(seurat_integrated, cells=seurat_integrated@meta.data$cells[seurat_integrated@meta.data$orig.ident %in% "R1"]), 
                         label.size =0.3,
                         pt.size=3, alpha=0) + 
                   stat_density2d(data=r45_day11_1BF_sgSpz, aes(x=x,y=y, alpha=factor(..level..)), fill='#fb9a99', lwd=0.07,
                                  color=c(rep(NA, 300), rep("red", 591)),
                                  geom="polygon", bins=6) +
                   stat_density2d(data=r45_day11_2BF_sgSpz, aes(x=x,y=y, alpha=factor(..level..)), fill='#cab2d6', lwd=0.07,
                                  color=c(rep(NA, 300), rep("purple", 593)),
                                  geom="polygon", bins=6) +
                   scale_alpha_manual(values=rep(c(0,0.5,0.5,0.5,0.5),6)) +
                   theme(legend.position = "none") + xlim(-9.4,6.6) + ylim(-8,8) + coord_fixed() + 
                   theme(text = element_text(size = 7.5),
                         axis.title = element_text(size = 7.5),
                         axis.text = element_text(size = 7.5)) + 
                   theme(aspect.ratio=1) +
                   labs(x= "UMAP1", y="UMAP2")) 

grid.draw(rbind(p2,p3,p4, size = "first"))


# Plot cell clusters
DimPlot(seurat_integrated, group_by="slingPseudotime_1", point_size=1.25) +
  geom_path(data=embedded, aes(x=umap_1, y=umap_2), size=0.8) + 
  geom_text(label="\u25BA", data=tail(head(embedded,2000),1),
            mapping=aes(x=umap_1, y=umap_2), size=4, angle=6) +
  geom_text(label="\u25BA", data=tail(head(embedded,4060),1),
            mapping=aes(x=umap_1, y=umap_2), size=4, angle=-50)

vln_df <- as.data.frame(seurat_integrated@assays$RNA$data[, grepl("^R4|R5", colnames(seurat_integrated@assays$RNA$data))]) %>% tibble::rownames_to_column('gene') %>% 
  pivot_longer(cols = starts_with("R"), values_to = "count") %>% 
  mutate(num_feeds =str_extract(name, ".BF"),
         biorep = str_extract(name, "^R."))
vln_df$cluster <- seurat_integrated@meta.data$leiden_res00015_mergedclusters[match(vln_df$name, seurat_integrated@meta.data$cells)]
vln_df$subcluster <- seurat_integrated@meta.data$leiden_res00015_clusters[match(vln_df$name, seurat_integrated@meta.data$cells)]

# Build replicate labels 
vln_df <- vln_df %>%
  mutate(
    replicate = case_when(
      grepl('^R4', name) ~ 'R4',
      grepl('^R5', name) ~ 'R5',
      TRUE ~ NA_character_))

# Cluster-level summaries
cluster_summary <- vln_df %>%
  group_by(cluster, num_feeds) %>%
  summarise(
    clusterMemSize = n_distinct(name),
    clusterCountSum = sum(count),
    .groups = 'drop')

vln_df <- left_join(
  vln_df,
  cluster_summary,
  by = c('cluster', 'num_feeds'))

# Subcluster-level summaries
subcluster_summary <- vln_df %>%
  group_by(subcluster, num_feeds) %>%
  summarise(
    subclusterMemSize = n_distinct(name),
    subclusterCountSum = sum(count),
    .groups = 'drop')

vln_df <- left_join(
  vln_df,
  subcluster_summary,
  by = c('subcluster', 'num_feeds'))

# Replicate + cluster summaries
rep_cluster_summary <- vln_df %>%
  group_by(cluster, replicate, num_feeds) %>%
  summarise(
    RepClusterCountSum = sum(count),
    RepClusterMemSize = n_distinct(name),
    .groups = 'drop')

vln_df <- left_join(
  vln_df,
  rep_cluster_summary,
  by = c('cluster', 'replicate', 'num_feeds'))

# Replicate + subcluster summaries
rep_subcluster_summary <- vln_df %>%
  group_by(subcluster, replicate, num_feeds) %>%
  summarise(
    RepSubclusterCountSum = sum(count),
    RepSubclusterMemSize = n_distinct(name),
    .groups = 'drop')

vln_df <- left_join(
  vln_df,
  rep_subcluster_summary,
  by = c('subcluster', 'replicate', 'num_feeds'))


vln_df.2 <- vln_df
noise <- rnorm(n = length(x = vln_df.2[, "count"])) / 100000
vln_df.2$count <- vln_df.2$count + noise
# violin plot without noise
#vln_df <- vln_df[vln_df$gene %in% mygenes$gene[c(3,4,5,9)],]

vln_sub <- vln_df[vln_df$gene %in% 
                   c(
                   #   "PF3D7-0408700", 
                   #   "PF3D7-1318800", 
                   #   "PF3D7-0708400", 
                   #   "PF3D7-1147000", ### slarp
                   # "PF3D7-0202100", 
                      "PF3D7-0629300", ### phospho
                   #   "PF3D7-1436300", 
                   #   "PF3D7-0930200", 
                   #   "PF3D7-1129100",
                    #  "PF3D7-1216600", ### celtos early
                    
                   #  "PF3D7-1016900", ### etramp10
                     # "PF3D7-1436300",
                     # "PF3D7-1336700",
                     #  "PF3D7-0829600",
                     # "PF3D7-1137800",
                     # 
                     # "PF3D7-1121700",
                   "PF3D7-0922200",

              #   "PF3D7-1461200",
              #   "PF3D7-1027100",
              # "PF3D7-1031400",
              #   "PF3D7-1351000",
              #   "PF3D7-1121700",
              #   "PF3D7-0202100",
              #   "PF3D7-1338900",
              #   "PF3D7-1428500",
              #   "PF3D7-0311300",
              #   "PF3D7-0726500",
             #   "PF3D7-0506300",
             #  "PF3D7-0818900",
             #  "PF3D7-0507300",
             #  "PF3D7-1148620",
             #  "PF3D7-0404900",
             # "PF3D7-1235100",
            #  "PF3D7-1148620", ### ribosomal
             # "PF3D7-0528400",
            #  "PF3D7-0818900",
            #   "PF3D7-0404800", # S10
             #  "PF3D7-1442600", # TREP
               "PF3D7-1147800" # MEABL
              ) & vln_df$count>=0,]


# mygenes <- c("PF3D7-1449600",
# "PF3D7-1338700",
# "PF3D7-0209500",
# "PF3D7-0622700",
# "PF3D7-0512400",
# "PF3D7-0803500",
# "PF3D7-1031900",
# "PF3D7-0209500",
# "PF3D7-1442700",
# "PF3D7-1351800",
# "PF3D7-0624800",
# "PF3D7-0522900",
# "PF3D7-1449600",
# "PF3D7-0718300",
# "PF3D7-1008600",
# "PF3D7-0404800",
# "PF3D7-0720400",
# "PF3D7-1442600",
# "PF3D7-1409100",
# "PF3D7-1236500",
# "PF3D7-0622800",
# "PF3D7-0503600",
# "PF3D7-0822700",
# "PF3D7-1456000",
# "PF3D7-1147800",
# "PF3D7-1462800",
# "PF3D7-0931700",
# "PF3D7-0931600",
# "PF3D7-0919300")

mygenes <- c(
  "PF3D7-0404800",
  "PF3D7-1442600",
  "PF3D7-1147800",
  "PF3D7-1351800",
  "PF3D7-1008600",
  "PF3D7-0209500",
  "PF3D7-0725600",
  "PF3D7-1236500",
  "PF3D7-0624800",
  "PF3D7-0931600",
  "PF3D7-1475400",
  "PF3D7-0315800",
  "PF3D7-1218100",
  "PF3D7-0622800",
  "PF3D7-0531600",
  "PF3D7-0726000",
  "PF3D7-0802000",
  "PF3D7-1456000",
  "PF3D7-0521600",
  "PF3D7-0718300",
  "PF3D7-1017700",
  "PF3D7-0911900",
  "PF3D7-0720400",
  "PF3D7-1148640",
  "PF3D7-1317200",
  "PF3D7-1462800",
  "PF3D7-1449600",
  "PF3D7-MIT04000",
  "PF3D7-1023000",
  "PF3D7-0532000",
  "PF3D7-0315900",
  "PF3D7-0408600",
  "PF3D7-1127500",
  "PF3D7-1023900",
  "PF3D7-0818600",
  "PF3D7-0925100",
  "PF3D7-0931700",
  "PF3D7-0822700",
  "PF3D7-1211000",
  "PF3D7-0704600",
  "PF3D7-0808200",
  "PF3D7-0930300",
  "PF3D7-1022300",
  "PF3D7-0814600",
  "PF3D7-1407700",
  "PF3D7-1371300",
  "PF3D7-0801900",
  "PF3D7-0620200",
  "PF3D7-0304600")
  



violin_g <- ggplot(vln_sub, aes(x = num_feeds, y = count, col=num_feeds, group=num_feeds)) + 
  geom_violin(aes(fill=num_feeds), col="black", alpha=0.45, lwd=0.45, adjust=1, trim=TRUE, scale="width")  + 
  geom_sina(scale="width", trim=TRUE, data=vln_sub[vln_sub$count>0.5,], aes(x = num_feeds, y = count, fill=num_feeds), alpha=0.3, col="black", shape=21, size=1.5) + 
  geom_pointrange(stat="summary",
                  fun.min = function(z) { quantile(z,0.25) },
                  fun.max = function(z) { quantile(z,0.5) },
                  fun = function(z) { quantile(z,0.5) }, aes(fill=num_feeds), lwd=1, col="black", shape=22, size=1.1) +
  geom_sina(scale="width", trim=TRUE, data=vln_sub[vln_sub$count<=0.5,], aes(x = num_feeds, y = count, fill=num_feeds), alpha=0.3, col="black", shape=21, size=1.5) +
  geom_pointrange(stat="summary",
                  fun.min = function(z) { quantile(z,0.5) },
                  fun.max = function(z) { quantile(z,0.75) },
                  fun = function(z) { quantile(z,0.5) }, aes(fill=num_feeds), lwd=1, col="black", shape=22, size=1.1) +
  facet_grid(paste(gene, str_extract(description, "^................"), sep="\n")~subcluster) + 
  #facet_grid(paste(gene, str_extract(description, "^................"), sep="\n")~.) + 
  scale_fill_manual(values=c("goldenrod1","seagreen3"))



vln_sub <- vln_df[vln_df$gene %in% 
                    c(
                      "PF3D7-1461200",
                      "PF3D7-1027100",
                      "PF3D7-1031400",
                      "PF3D7-1351000",
                      "PF3D7-1121700",
                      "PF3D7-0202100",
                      "PF3D7-1338900",
                      "PF3D7-1428500",
                      "PF3D7-0311300",
                      "PF3D7-0726500",
                      "PF3D7-0506300",
                      "PF3D7-0818900"
                      ),]

mygenes <- c("PF3D7-1148620",
"PF3D7-1217400",
"PF3D7-1216600",
"PF3D7-0629300",
"PF3D7-0616500",
"PF3D7-1310700",
"PF3D7-1312500",
"PF3D7-0922200",
"PF3D7-1016900",
"PF3D7-1142900",
"PF3D7-1244600",
"PF3D7-1444800",
"PF3D7-0932200",
"PF3D7-1322400",
"PF3D7-1327300",
"PF3D7-0408700",
"PF3D7-0206200",
"PF3D7-0316600",
"PF3D7-1212800",
"PF3D7-1446900",
"PF3D7-1116700",
"PF3D7-1218000",
"PF3D7-1318800",
"PF3D7-1114100",
"PF3D7-0519700",
 "PF3D7-0206100",
 "PF3D7-1329700",
 "PF3D7-1147000",
 "PF3D7-0511200",
 "PF3D7-1339300",
 "PF3D7-1231200",
 "PF3D7-1336700",
 "PF3D7-0723300",
 "PF3D7-1206700",
 "PF3D7-1331000",
 "PF3D7-0503400",
 "PF3D7-1208900",
 "PF3D7-1365700",
 "PF3D7-1105100",
 "PF3D7-1333100",
 "PF3D7-1324900",
 "PF3D7-0708400",
 "PF3D7-1022700",
 "PF3D7-1307500",
 "PF3D7-1133400",
 "PF3D7-0214600",
 "PF3D7-0528400",
 "PF3D7-1244700",
 "PF3D7-0507300",
 "PF3D7-1209300",
 "PF3D7-1112100",
 "PF3D7-1368700",
 "PF3D7-1246200",
 "PF3D7-1404100",
 "PF3D7-0407600")

sig_norm <- as.data.frame(seurat_integrated@assays$RNA$data[, grepl("^R4|R5", colnames(seurat_integrated@assays$RNA$data))]) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% mygenes[1:25]  & !gene %in% c("PF3D7-1216600","PF3D7-0408700","PF3D7-0616500"))
rownames(sig_norm) <- sig_norm$gene


mycellorder <- c(
  # colnames(sig_norm)[colnames(sig_norm) %in% seurat_integrated@meta.data$cells[
  #  seurat_integrated@meta.data$leiden_res00015_mergedclusters %in% "Ld1" & 
  #    grepl("^R[4,5].*1BF", seurat_integrated@meta.data$cells)]]
  #  ,colnames(sig_norm)[colnames(sig_norm) %in% seurat_integrated@meta.data$cells[
  #  seurat_integrated@meta.data$leiden_res00015_mergedclusters %in% "Ld1" & 
  #    grepl("^R[4,5].*2BF", seurat_integrated@meta.data$cells)]],
  colnames(sig_norm)[colnames(sig_norm) %in% seurat_integrated@meta.data$cells[
    seurat_integrated@meta.data$leiden_res00015_mergedclusters %in% "Ld2+3" &
      grepl("^R[4,5].*1BF", seurat_integrated@meta.data$cells)]]
  ,colnames(sig_norm)[colnames(sig_norm) %in% seurat_integrated@meta.data$cells[
    seurat_integrated@meta.data$leiden_res00015_mergedclusters %in% "Ld2+3" &
      grepl("^R[4,5].*2BF", seurat_integrated@meta.data$cells)]]
)
sig_norm <- sig_norm[, c("gene", mycellorder)]


top25 <- scran_fM_up1BF_repset45$gene[scran_fM_up1BF_repset45$TRUE.logFC.FALSE<(-log2(1.3))][1:26]
top25 <- scran_fM_up2BF_repset45$gene[scran_fM_up2BF_repset45$TRUE.logFC.FALSE>(log2(1.3))][1:33]
top25 <- annotations$gene_name[grepl("TRAMP", annotations$description) & annotations$gene_name %in% vln_df$gene]
top25 <- c("PF3D7-1218000",
           #"PF3D7-0206200",
           #"PF3D7-0316600",
           "PF3D7-1310700",
           "PF3D7-1016900",
      #     "PF3D7-1147000"
           
           "PF3D7-1148620",
      #    "PF3D7-1336700",
           #"PF3D7-0404800",
      #     "PF3D7-1147800",
      #     "PF3D7-0725600"
           #"PF3D7-0531600",
           #"PF3D7-0726000"
           )

annotations[annotations$gene_name %in% top25,]


vln_sub <- vln_df[vln_df$gene %in% top25,]

RepsubclusterMemSizeOnlyPos <- vln_sub %>% mutate(positive = case_when(count > 0 ~ 1, .default=0)) %>%
  dplyr::group_by(subcluster, biorep, gene, num_feeds) %>%
  dplyr::summarise(RepsubclusterMemSizeOnlyPos = sum(positive))

vln_sub$RepsubclusterMemSizeOnlyPos <- RepsubclusterMemSizeOnlyPos$RepsubclusterMemSizeOnlyPos[match(
  paste(vln_sub$gene,
        vln_sub$subcluster,
        vln_sub$biorep,
        vln_sub$num_feeds),
  paste(RepsubclusterMemSizeOnlyPos$gene,
        RepsubclusterMemSizeOnlyPos$subcluster,
        RepsubclusterMemSizeOnlyPos$biorep,
        RepsubclusterMemSizeOnlyPos$num_feeds))]


temp <- as.data.frame(rep(top25,12))
myplot <- ggplot(vln_sub, aes(fill = num_feeds, alpha=biorep,
                    #weight = count/(RepsubclusterMemSizeOnlyPos),
                    weight = count/(RepSubclusterMemSize),
                    x=subcluster)) + 
  geom_bar(position = position_dodge(preserve = "single")) + 
  
  facet_wrap(.~factor(gene, levels=top25)) + 
  scale_fill_manual(values=c("goldenrod1","seagreen3")) +
  scale_alpha_manual(values=c(0.6,1)) +
  labs(y="Mean UMI\n")
g <- ggplot_build(myplot)
g$data[[1]]$subcluster <- rep(c(rep("Ld1", 4),rep("Ld2", 4),rep("Ld3", 4)),length(top25))
g$data[[1]]$gene <- temp[ order(match(temp$`rep(top25, 12)`, top25)), ]
g$data[[1]]$num_feeds <- rep(c(rep("1BF", 2), rep("2BF", 2)),3*length(top25))
g$data[[1]]$count <- g$data[[1]]$y
g$data[[1]]$subclusterMemSize <- 1
g$data[[1]]$biorep <- rep(c("R4", "R5"),6*length(top25))

#ggsave("fM_1-25_bars_3subclus_weighted_by_numcells_MergeRepBars.pdf", width=15, height=25)
g2 <- ggplot() + 
  geom_bar(vln_sub, mapping=aes(fill = num_feeds,
                        weight = count/(subclusterMemSize),
                        x=subcluster), position = position_dodge(preserve = "single"), color="black") + 
   geom_segment(cbind(g$data[[1]][g$data[[1]]$biorep %in% "R4",], g$data[[1]]$y[g$data[[1]]$biorep %in% "R5"]),
               mapping=aes(y=y, yend=g$data[[1]]$y[g$data[[1]]$biorep %in% "R5"], group = num_feeds, x=subcluster),
               position = position_dodge(width=0.9), lwd=0.2, lty="dashed") +
  geom_point(g$data[[1]][g$data[[1]]$biorep %in% "R4",], mapping=aes(y=y, group = num_feeds, x=subcluster),
             position = position_dodge(width=0.9), size=3.2, fill="black", shape=21) +
  geom_point(g$data[[1]][g$data[[1]]$biorep %in% "R5",], mapping=aes(y=y, group = num_feeds, x=subcluster),
             position = position_dodge(width=0.9), size=3.2, fill="grey", shape=21) +
 facet_wrap(.~factor(gene, levels=top25)) + 
  scale_fill_manual(values=c("goldenrod1","seagreen3")) +
  scale_color_manual(values=c("goldenrod1","seagreen3")) +
  labs(y="Mean UMI\n", x="\nSubcluster") + theme(strip.text = element_text(size=20)) +
  theme(axis.text = element_text(size=20)) + theme(axis.title = element_text(size=20))
#dev.off()


g3 <- ggplot_build(g2)
g3[[1]][[2]]$y_old <- g3[[1]][[2]]$y
for (i in (1:length(g3[[1]][[2]]$yend))[1:length(g3[[1]][[2]]$yend) %% 2 == 1]){
g3[[1]][[2]]$yend[i] <- g3[[1]][[2]]$ymax[i+1]}

g3[[1]][[2]]$xend <- NA
for (i in (1:length(g3[[1]][[2]]$yend))[1:length(g3[[1]][[2]]$yend) %% 2 == 1]){
g3[[1]][[2]]$xend[i] <- g3[[1]][[2]]$x[i+1]}


 for (i in (1:length(g3[[1]][[2]]$yend))[1:length(g3[[1]][[2]]$yend) %% 2 == 0]){
  g3[[1]][[2]]$y[i] <- g3[[1]][[4]]$y[i-1]}

 for (i in (1:length(g3[[1]][[2]]$yend))[1:length(g3[[1]][[2]]$yend) %% 2 == 0]){
 g3[[1]][[2]]$xend[i] <- g3[[1]][[2]]$x[i]}

 for (i in (1:length(g3[[1]][[2]]$yend))[1:length(g3[[1]][[2]]$yend) %% 2 == 0]){
   g3[[1]][[2]]$x[i] <- g3[[1]][[2]]$x[i-1]}

grid::grid.draw(ggplot_gtable(g3))


gx <- ggplot(vln_sub %>% filter(gene %in% top25), aes(x = num_feeds, y = count, col=num_feeds, group=num_feeds)) + 
  # geom_pointrange(stat="summary",
  #                 fun.min = function(z) { quantile(z,0.9) },
  #                 fun.max = function(z) { quantile(z,0.1)  },
  #                 fun = function(z) { median(z) }, aes(x = num_feeds, y = count, fill=num_feeds, col=num_feeds), lwd=2, shape="-", size=0) +
  geom_violin(aes(fill=num_feeds), col="black", alpha=0.45, lwd=0.45, adjust=1, trim=TRUE, scale="width")  + 
  geom_sina(scale="width", trim=TRUE, data=vln_sub %>% filter(gene %in% top25 & count>0.5), aes(x = num_feeds, y = count, fill=num_feeds), alpha=0.3, col="black", shape=21, size=1.1, stroke=0.2) + 
  geom_sina(scale="width", trim=TRUE, data=vln_sub %>% filter(gene %in% top25 & count<=0.5), aes(x = num_feeds, y = count, fill=num_feeds), alpha=0.3, col="black", shape=21, size=1.1, stroke=0.2) +
  
  
  #facet_grid(.~subcluster)
   # facet_grid(factor(paste(gene, "\n", str_extract(description, "^......."), "...", sep=""),
   #                   levels=paste(top25, "\n", str_extract(annotations$description[match(top25, annotations$gene_name)], "^......."), "...", sep="")
   #                   )~subcluster) + 
  facet_grid(factor(gene, levels=top25)~subcluster) +
  scale_fill_manual(values=c(col1BF,col2BF)) + scale_color_manual(values=c(col1BF,col2BF))

for (mynum in 1:length(top25)){
  gx <- gx + 
    # need to streamline facettes 
    geom_line(data=as.data.frame(cbind(x=c(1,2), y=as.data.frame(g3[[1]][[3]])[(1:2)+6*(mynum-1),"y"])) %>% mutate(subcluster="Ld1", num_feeds=c("1BF", "2BF"), gene=top25[mynum]) %>% left_join(annotations[, c("gene_name", "description")], by=c("gene"="gene_name")), aes(x=num_feeds, y=y, group=subcluster), size=0.2, color="black") +
    geom_line(data=as.data.frame(cbind(x=c(1,2), y=as.data.frame(g3[[1]][[4]])[(1:2)+6*(mynum-1),"y"])) %>% mutate(subcluster="Ld1", num_feeds=c("1BF", "2BF"), gene=top25[mynum]) %>% left_join(annotations[, c("gene_name", "description")], by=c("gene"="gene_name")), aes(x=num_feeds, y=y, group=subcluster), size=0.2, color="black") +
    
    geom_line(data=as.data.frame(cbind(x=c(1,2), y=as.data.frame(g3[[1]][[3]])[(3:4)+6*(mynum-1),"y"])) %>% mutate(subcluster="Ld2", num_feeds=c("1BF", "2BF"), gene=top25[mynum]) %>% left_join(annotations[, c("gene_name", "description")], by=c("gene"="gene_name")), aes(x=num_feeds, y=y, group=subcluster), size=0.2, color="black") +
    geom_line(data=as.data.frame(cbind(x=c(1,2), y=as.data.frame(g3[[1]][[4]])[(3:4)+6*(mynum-1),"y"])) %>% mutate(subcluster="Ld2", num_feeds=c("1BF", "2BF"), gene=top25[mynum]) %>% left_join(annotations[, c("gene_name", "description")], by=c("gene"="gene_name")), aes(x=num_feeds, y=y, group=subcluster), size=0.2, color="black") +
    
    geom_line(data=as.data.frame(cbind(x=c(1,2), y=as.data.frame(g3[[1]][[3]])[(5:6)+6*(mynum-1),"y"])) %>% mutate(subcluster="Ld3", num_feeds=c("1BF", "2BF"), gene=top25[mynum]) %>% left_join(annotations[, c("gene_name", "description")], by=c("gene"="gene_name")), aes(x=num_feeds, y=y, group=subcluster), size=0.2, color="black") +
    geom_line(data=as.data.frame(cbind(x=c(1,2), y=as.data.frame(g3[[1]][[4]])[(5:6)+6*(mynum-1),"y"])) %>% mutate(subcluster="Ld3", num_feeds=c("1BF", "2BF"), gene=top25[mynum]) %>% left_join(annotations[, c("gene_name", "description")], by=c("gene"="gene_name")), aes(x=num_feeds, y=y, group=subcluster), size=0.2, color="black") +
    
    geom_point(data=as.data.frame(cbind(x=c(1,2), y=as.data.frame(g3[[1]][[3]])[(1:2)+6*(mynum-1),"y"])) %>% mutate(subcluster="Ld1", num_feeds=c("1BF", "2BF"), gene=top25[mynum]) %>% left_join(annotations[, c("gene_name", "description")], by=c("gene"="gene_name")), aes(x=num_feeds, y=y, group=subcluster, fill=num_feeds), shape=21, size=2, stroke=0.3, color="black")  +
    geom_point(data=as.data.frame(cbind(x=c(1,2), y=as.data.frame(g3[[1]][[4]])[(1:2)+6*(mynum-1),"y"])) %>% mutate(subcluster="Ld1", num_feeds=c("1BF", "2BF"), gene=top25[mynum]) %>% left_join(annotations[, c("gene_name", "description")], by=c("gene"="gene_name")), aes(x=num_feeds, y=y, group=subcluster, fill=num_feeds), shape=22, size=2, stroke=0.3, color="black")  +
    
    geom_point(data=as.data.frame(cbind(x=c(1,2), y=as.data.frame(g3[[1]][[3]])[(3:4)+6*(mynum-1),"y"])) %>% mutate(subcluster="Ld2", num_feeds=c("1BF", "2BF"), gene=top25[mynum]) %>% left_join(annotations[, c("gene_name", "description")], by=c("gene"="gene_name")), aes(x=num_feeds, y=y, group=subcluster, fill=num_feeds), shape=21, size=2, stroke=0.3, color="black")  +
    geom_point(data=as.data.frame(cbind(x=c(1,2), y=as.data.frame(g3[[1]][[4]])[(3:4)+6*(mynum-1),"y"])) %>% mutate(subcluster="Ld2", num_feeds=c("1BF", "2BF"), gene=top25[mynum]) %>% left_join(annotations[, c("gene_name", "description")], by=c("gene"="gene_name")), aes(x=num_feeds, y=y, group=subcluster, fill=num_feeds), shape=22, size=2, stroke=0.3, color="black")  +
    
    geom_point(data=as.data.frame(cbind(x=c(1,2), y=as.data.frame(g3[[1]][[3]])[(5:6)+6*(mynum-1),"y"])) %>% mutate(subcluster="Ld3", num_feeds=c("1BF", "2BF"), gene=top25[mynum]) %>% left_join(annotations[, c("gene_name", "description")], by=c("gene"="gene_name")), aes(x=num_feeds, y=y, group=subcluster, fill=num_feeds), shape=21, size=2, stroke=0.3, color="black")  +
    geom_point(data=as.data.frame(cbind(x=c(1,2), y=as.data.frame(g3[[1]][[4]])[(5:6)+6*(mynum-1),"y"])) %>% mutate(subcluster="Ld3", num_feeds=c("1BF", "2BF"), gene=top25[mynum]) %>% left_join(annotations[, c("gene_name", "description")], by=c("gene"="gene_name")), aes(x=num_feeds, y=y, group=subcluster, fill=num_feeds), shape=22, size=2, stroke=0.3, color="black")  }

gx <- gx + theme_classic()
#gx <- gx + coord_cartesian(ylim=c(3.25,3.67))
#gx <- gx + coord_cartesian(ylim=c(0,5.7))
gxx <- ggplot_gtable(ggplot_build(gx))
stripr <- which(grepl('strip-t', gxx$layout$name))
fills2 <- c(colLd1,colLd2,colLd3)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', gxx$grobs[[i]]$grobs[[1]]$childrenOrder))
  gxx$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills2[k]
  k <- k+1
}

gxx$grobs[[33]]$grobs[[1]]$children[[2]]$children[[1]]$label <- "1"
gxx$grobs[[34]]$grobs[[1]]$children[[2]]$children[[1]]$label <- "2"
gxx$grobs[[35]]$grobs[[1]]$children[[2]]$children[[1]]$label <- "3"

grid::grid.draw(gxx)









cbind(g$data[[1]]$y[g$data[[1]]$biorep %in% "R4" & g$data[[1]]$num_feeds %in% "1BF"], g$data[[1]]$y[g$data[[1]]$biorep %in% "R4" & g$data[[1]]$num_feeds %in% "2BF"])


ggplot(vln_df %>% filter(gene %in% "PF3D7-1016900"), aes(fill = num_feeds, alpha=biorep,
                    weight = (10^6)*count/(RepSubclusterMemSize),
                    x=subcluster)) + 
  geom_bar(position = position_dodge(preserve = "single")) + 
  scale_fill_manual(values=c("goldenrod1","seagreen3")) +
  scale_alpha_manual(values=c(0.6,1)) +
  labs(y="Gene transcripts / Total TPM\n") 

ggsave("fM_1-25_bars_3subclus_weighted_by_numcells.pdf", width=10, height=10)
ggplot(vln_sub, aes(fill = num_feeds, alpha=biorep,
                    weight = (10^3)*count/(RepSubclusterMemSize),
                    x=subcluster)) + 
  geom_bar(position = position_dodge(preserve = "single")) + 
  facet_wrap(.~factor(gene, levels=top25)) + 
  scale_fill_manual(values=c("goldenrod1","seagreen3")) +
  scale_alpha_manual(values=c(0.6,1)) +
  labs(y="Gene transcripts / 1k cells\n")
dev.off()

ggsave("fM_1-25_bars_3subclus_weighted_by_TPM.pdf", width=10, height=10)
ggplot(vln_sub, aes(fill = num_feeds, alpha=biorep,
                    weight = (10^6)*count/(RepSubclusterCountSum),
                    x=subcluster)) + 
  geom_bar(position = position_dodge(preserve = "single")) + 
  facet_wrap(.~factor(gene, levels=top25)) + 
  scale_fill_manual(values=c("goldenrod1","seagreen3")) +
  scale_alpha_manual(values=c(0.6,1)) +
  labs(y="Gene transcripts / Total TPM\n")
dev.off()

ggsave("fM_1-25_bars_3subclus_no_weighting.pdf", width=10, height=10)
ggplot(vln_sub, aes(fill = num_feeds, alpha=biorep,
                    weight = count,
                    x=subcluster)) + 
  geom_bar(position = position_dodge(preserve = "single")) + 
  facet_wrap(.~factor(gene, levels=top25)) + 
  scale_fill_manual(values=c("goldenrod1","seagreen3")) +
  scale_alpha_manual(values=c(0.6,1)) +
  labs(y="Gene transcripts (no weighting)\n")
dev.off()


ggsave("fM_1-25_bars_2clus_weighted_by_numcells.pdf", width=10, height=10)
ggplot(vln_sub, aes(fill = num_feeds, alpha=biorep,
                    weight = (10^3)*count/(RepClusterMemSize),
                    x=cluster)) + 
  geom_bar(position = position_dodge(preserve = "single")) + 
  facet_wrap(.~factor(gene, levels=top25)) + 
  scale_fill_manual(values=c("goldenrod1","seagreen3")) +
  scale_alpha_manual(values=c(0.6,1)) +
  labs(y="Gene transcripts / 1k cells\n")
dev.off()

ggsave("fM_1-25_bars_2clus_weighted_by_TPM.pdf", width=10, height=10)
ggplot(vln_sub, aes(fill = num_feeds, alpha=biorep,
                    weight = (10^6)*count/(RepClusterCountSum),
                    x=cluster)) + 
  geom_bar(position = position_dodge(preserve = "single")) + 
  facet_wrap(.~factor(gene, levels=top25)) + 
  scale_fill_manual(values=c("goldenrod1","seagreen3")) +
  scale_alpha_manual(values=c(0.6,1)) +
  labs(y="Gene transcripts / Total TPM\n")
dev.off()

ggsave("fM_1-25_bars_2clus_no_weighting.pdf", width=10, height=10)
  ggplot(vln_sub, aes(fill = num_feeds, alpha=biorep,
                    weight = count,
                    x=cluster)) + 
  geom_bar(position = position_dodge(preserve = "single")) + 
  facet_wrap(.~factor(gene, levels=top25)) + 
  scale_fill_manual(values=c("goldenrod1","seagreen3")) +
  scale_alpha_manual(values=c(0.6,1)) +
  labs(y="Gene transcripts (no weighting)\n")
dev.off()


ggplot(vln_sub, aes(fill = num_feeds, alpha=biorep,
                    weight = count/RepSubclusterMemSize,
                    x=subcluster)) + 
  geom_bar(position = position_dodge(preserve = "single")) + 
  facet_wrap(.~gene) + 
  scale_fill_manual(values=c("goldenrod1","seagreen3")) +
  scale_alpha_manual(values=c(0.6,1)) +
  labs(y="Gene transcripts / cell count\n")


# Plot cell clusters
DimPlot(seurat_integrated, group.by="leiden_res00015_clusters",
        reduction = "umap",
        label = TRUE,
        label.size = 10,
        pt.size=1.5)
# Plot cell clusters
DimPlot(seurat_integrated, group.by="leiden_res00015_clusters",
        reduction = "pca",
        label = TRUE,
        label.size = 10,
        pt.size=1.5)


seurat_integrated <- IntegrateLayers(object = seurat_integrated, method = JointPCAIntegration, normalization.method="SCT", assay="SCT",
                # k.weight=64,
                orig.reduction = "sct.pca", new.reduction = "integrated.sct.jointpca",
                verbose = TRUE)
# Compute UMAP embedding
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:9, reduction = "integrated.sct.jointpca", reduction.name="integrated.sct.jointpca.umap")

# Plot gene expression on UMAP
FeaturePlot(seurat_integrated,
            features="PF3D7-1032700",
            reduction = "integrated.sct.rpca.umap", order=T, dims=c(1,2),
            label = FALSE,
            pt.size=2.5)
# Plot gene expression on UMAP
FeaturePlot(seurat_integrated,
            features="PF3D7-1032700",
            reduction = "integrated.jointpca.umap", order=T, dims=c(1,2),
            label = FALSE,
            pt.size=2.5)
# Plot gene expression on UMAP
FeaturePlot(seurat_integrated,
            features="PF3D7-1032700",
            reduction = "integrated.sct.jointpca.umap", order=T, dims=c(1,2),
            label = FALSE,
            pt.size=2.5)
# Plot gene expression on UMAP
FeaturePlot(seurat_integrated,
            features="PF3D7-1032700",
            reduction = "pca", order=T, dims=c(1,2),
            label = FALSE,
            pt.size=2.5)
# Plot gene expression on UMAP
FeaturePlot(seurat_integrated,
            features="PF3D7-1032700",
            reduction = "umap", order=T, dims=c(1,2),
            label = FALSE,
            pt.size=2.5)


# Select the RNA counts slot to be the default assay and normalize RNA data

DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE) ### this is the function which seems buggy in some versions of Seurat


### Take a look at cluster assignments in umap at broad resolution (integrated_snn_res.0.3)
### Using finer resolution integrated_snn_res.0.4 etc can also be explored

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:2,
                                   reduction="umap")

# Determine the clusters for various resolutions                                
# Run graph-based clustering
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8))


# Plot cell clusters
DimPlot(seurat_integrated, group.by="RNA_snn_res.0.05",
        reduction = "umap",
        label = TRUE,
        label.size = 10,
        pt.size=1.5)


# Plot gene expression on UMAP
FeaturePlot(seurat_integrated, 
            reduction = "integrated.rpca.umap", 
            features = "nUMI",
            pt.size = 2.5, label.size=0, 
            order = TRUE,
            label = TRUE)

### Alternatively look at leiden-based clustering results created above

# Plot cell clusters
DimPlot(seurat_integrated, group.by="leiden_res00015_clusters",
        reduction = "pca", dimensions=2:3,
        label = TRUE,
        label.size = 15,
        pt.size=1.5)


# Plot gene expression on UMAP
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("PF3D7-1032700" ), 
            order = TRUE,
            min.cutoff = 'q10', pt.size = 1.6)


### You can color samples by other metadata columns for example "inject" or "num_feeds"
### check colnames(seurat_integrated@meta.data) for possible grouping variables

# Plot cell clusters
DimPlot(seurat_integrated,  group.by="inject",
        reduction = "umap",
        label = FALSE,
        label.size = 5,
        pt.size=1.5)

#########################################################################################################################################
################ Seurat FindMarkers DE section ##########################################################################################
#########################################################################################################################################

adj_pval_to_use=0.05
log2FC_to_use=0.584962501

Seurat_FM_upin2BF_only_d11_R4_R5 <- get_wilcox_DE_upingroup2(replicates_to_use=c("R4","R5"),
                                                          dpbm_to_use="d11", num_feeds_to_use = c("1BF", "2BF"), adj_pval_to_use=adj_pval_to_use, log2FC_to_use=log2FC_to_use,
                                                          groupingvar="num_feeds", group1="1BF", group2="2BF")
Seurat_FM_upin2BF_only_d11_R1_R2_R3 <- get_wilcox_DE_upingroup2(replicates_to_use=c("R1","R2","R3"),
                                                             dpbm_to_use="d11", num_feeds_to_use = c("1BF", "2BF"), adj_pval_to_use=adj_pval_to_use, log2FC_to_use=log2FC_to_use,
                                                             groupingvar="num_feeds", group1="1BF", group2="2BF")
Seurat_FM_upin2BF_only_tot <- get_wilcox_DE_upingroup2(replicates_to_use=c("R1","R2","R3","R4","R5"),
                                                    dpbm_to_use="d11", num_feeds_to_use = c("1BF", "2BF"), adj_pval_to_use=adj_pval_to_use, log2FC_to_use=log2FC_to_use,
                                                    groupingvar="num_feeds", group1="1BF", group2="2BF")

#########################################################################################################################################
################ scran psuedobulk DE section ############################################################################################
#########################################################################################################################################

FDR_to_use=0.05
lnFC_to_use=0.4054651081 ### log2FC = log(1.5) = 0.4054651081
pval_to_use=0.05

#### First using the broad, 4-cluster resolution integrated_snn_res.0.3

get_pseubul_DE_broad_res(FDR_to_use=FDR_to_use, lnFC_to_use = lnFC_to_use, pval_to_use = pval_to_use)

### The function returns an object for each cluster, for example r0_UPREG_2BF.broad for cluster 0
### These object are contained within a list called pseubul_DE_broad_res_results_up2bf and a list called pseubul_DE_broad_res_results_up1bf

### we can concatenate significant DE gene names across clusters 

scran_pseubul_FM_upin2BF <- as.data.frame(unique(c(pseubul_DE_broad_res_results_up2bf$clus1_signal$gene, pseubul_DE_broad_res_results_up2bf$clus2_signal$gene, pseubul_DE_broad_res_results_up2bf$clus3_signal$gene)))
colnames(scran_pseubul_FM_upin2BF) <- "gene"

scran_pseubul_FM_upin1BF <- as.data.frame(unique(c(pseubul_DE_broad_res_results_up1bf$clus1_signal$gene, pseubul_DE_broad_res_results_up1bf$clus2_signal$gene, pseubul_DE_broad_res_results_up1bf$clus3_signal$gene)))
colnames(scran_pseubul_FM_upin1BF) <- "gene"

### you can print a graphical summary pdf which highlights expression changes across clusters for genes with significant upreg in 2bf
### name it as you wish...

make_PDF_pseubul_DE_broad_res(filename="pseubul_DE_broad_res_july31_2024.pdf")

#### Next lets try using the higher, 6-cluster resolution integrated_snn_res.0.4

# get_pseubul_DE_hi_res(FDR_to_use=FDR_to_use, lnFC_to_use = lnFC_to_use, pval_to_use = pval_to_use)
# 
# make_PDF_pseubul_DE_hi_res(filename="name_it_as_you_wish2.pdf")


### only using repset 45
only_d11 <- subset(x = seurat_integrated, subset= (repset=="R45"))
M <- SetIdent(only_d11, value = "num_feeds")
DefaultAssay(M) <- "RNA"
merged_broad_repset45 <- as.SingleCellExperiment(only_d11)
merged_broad_repset45$num_feeds <- merged_broad_repset45$num_feeds=="2BF"
abundances_broad <- unclass(table(merged_broad_repset45$leiden_res00015_clusters, merged_broad_repset45$sample))
colLabels(merged_broad_repset45) <- factor(merged_broad_repset45$leiden_res00015_clusters)
m.out <- findMarkers(merged_broad_repset45, block=merged_broad_repset45$orig.ident, groups=merged_broad_repset45$num_feeds, assay.type="logcounts")
mdf <- as.data.frame(m.out@listData)
mdf$gene <- rownames(mdf)
up2bf_repset45 <- mdf[mdf$FALSE.FDR<0.05 & (mdf$FALSE.logFC.TRUE < (-0.25)),] %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
scran_fM_up2BF_repset45 <- up2bf_repset45[,c(7,8,10,11,12)] %>% arrange(-TRUE.logFC.FALSE)
up1bf_repset45 <- mdf[mdf$FALSE.FDR<0.05 & (mdf$FALSE.logFC.TRUE > (0.25)),] %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
scran_fM_up1BF_repset45 <- up1bf_repset45[,c(7,8,10,11,12)] %>% arrange(TRUE.logFC.FALSE)


summed_broad_repset45 <- aggregateAcrossCells(merged_broad_repset45, use.assay.type = 1,
                                               id=colData(merged_broad_repset45)[,c("leiden_res00015_clusters", "sample")])

summed.filt.broad_repset45 <- summed_broad_repset45[,summed_broad_repset45$ncells >= 10]

de.results.broad_repset45 <- pseudoBulkDGE(summed.filt.broad_repset45,
                                            label=summed.filt.broad_repset45$leiden_res00015_clusters,
                                            design=~factor(orig.ident) + num_feeds,
                                            coef="num_feedsTRUE",
                                            condition=summed.filt.broad_repset45$num_feeds)

r1.broad_repset45 <- as.data.frame(de.results.broad_repset45[["Ld1"]])
r1.broad_repset45.noNA <- r1.broad_repset45[complete.cases(r1.broad_repset45), ]
r2.broad_repset45 <- as.data.frame(de.results.broad_repset45[["Ld2"]])
r2.broad_repset45.noNA <- r2.broad_repset45[complete.cases(r2.broad_repset45), ]
r3.broad_repset45 <- as.data.frame(de.results.broad_repset45[["Ld3"]])
r3.broad_repset45.noNA <- r3.broad_repset45[complete.cases(r3.broad_repset45), ]

FDR_to_use=0.05
log2FC_to_use = log2(4/3)
pval_to_use=0.05

r1_UPREG_2BF.broad_repset45 <- r1.broad_repset45.noNA[r1.broad_repset45.noNA$PValue<pval_to_use & r1.broad_repset45.noNA$FDR<FDR_to_use & r1.broad_repset45.noNA$logFC>(log2FC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

r1_DOWNREG_2BF.broad_repset45 <- r1.broad_repset45.noNA[r1.broad_repset45.noNA$PValue<pval_to_use & r1.broad_repset45.noNA$FDR<FDR_to_use & r1.broad_repset45.noNA$logFC<(-log2FC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

r2_UPREG_2BF.broad_repset45 <- r2.broad_repset45.noNA[r2.broad_repset45.noNA$PValue<pval_to_use & r2.broad_repset45.noNA$FDR<FDR_to_use & r2.broad_repset45.noNA$logFC>(log2FC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

r2_DOWNREG_2BF.broad_repset45 <- r2.broad_repset45.noNA[r2.broad_repset45.noNA$PValue<pval_to_use & r2.broad_repset45.noNA$FDR<FDR_to_use & r2.broad_repset45.noNA$logFC<(-log2FC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

r3_UPREG_2BF.broad_repset45 <- r3.broad_repset45.noNA[r3.broad_repset45.noNA$PValue<pval_to_use & r3.broad_repset45.noNA$FDR<FDR_to_use & r3.broad_repset45.noNA$logFC>(log2FC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

r3_DOWNREG_2BF.broad_repset45 <- r3.broad_repset45.noNA[r3.broad_repset45.noNA$PValue<pval_to_use & r3.broad_repset45.noNA$FDR<FDR_to_use & r3.broad_repset45.noNA$logFC<(-log2FC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

for (i in c(1,2,3)){
  if (dim(list(r1_UPREG_2BF.broad_repset45,r2_UPREG_2BF.broad_repset45,r3_UPREG_2BF.broad_repset45)[[i]])[1]==0){
    print(paste("FYI, cluster ", i, " has no significant 2bf upreg"))}}
for (i in c(1,2,3)){
  if (dim(list(r1_DOWNREG_2BF.broad_repset45,r2_DOWNREG_2BF.broad_repset45,r3_DOWNREG_2BF.broad_repset45)[[i]])[1]==0){
    print(paste("FYI, cluster ", i, " has no significant 1bf upreg"))}}

pseubul_DE_broad_repset45_res_results_up2bf <- list(r1_UPREG_2BF.broad_repset45,r2_UPREG_2BF.broad_repset45,r3_UPREG_2BF.broad_repset45)
names(pseubul_DE_broad_repset45_res_results_up2bf) <- c("clus1_signal","clus2_signal","clus3_signal")



pseubul_DE_broad_repset45_res_results_up1bf <- list(r1_DOWNREG_2BF.broad_repset45,r2_DOWNREG_2BF.broad_repset45,r3_DOWNREG_2BF.broad_repset45)
names(pseubul_DE_broad_repset45_res_results_up1bf) <- c("clus1_signal","clus2_signal","clus3_signal")

bind_rows(pseubul_DE_broad_repset45_res_results_up2bf, .id="id") %>% arrange(-logFC)
bind_rows(pseubul_DE_broad_repset45_res_results_up1bf, .id="id") %>% arrange(logFC)

### only using repset 45 and large clusters

# counts <- GetAssayData(seurat_integrated, layer="counts", assay="RNA")  
# genes.filter <- rownames(counts[-(which(rownames(counts) %in% annotations$gene_name[annotations$ebi_biotype %in% "rRNA"])),])
# counts.sub <- counts[genes.filter,]
# seurat_integrated_norRNA  <- CreateSeuratObject(counts=counts.sub)
# seurat_integrated_norRNA <- NormalizeData(seurat_integrated_norRNA, verbose = FALSE)
# seurat_integrated_norRNA <- FindVariableFeatures(seurat_integrated_norRNA)
# seurat_integrated_norRNA <- ScaleData(seurat_integrated_norRNA)
# seurat_integrated_norRNA@meta.data$cells <- rownames(seurat_integrated_norRNA@meta.data)
# seurat_integrated_norRNA@meta.data$num_feeds <- str_extract(seurat_integrated_norRNA@meta.data$cells, ".BF")
# seurat_integrated_norRNA@meta.data$leiden_res00015_clusters <- seurat_integrated@meta.data$leiden_res00015_clusters[match(seurat_integrated_norRNA@meta.data$cells, seurat_integrated@meta.data$cells)]
# seurat_integrated_norRNA@meta.data$leiden_res00015_mergedclusters <- seurat_integrated@meta.data$leiden_res00015_mergedclusters[match(seurat_integrated_norRNA@meta.data$cells, seurat_integrated@meta.data$cells)]
# seurat_integrated_norRNA@meta.data$repset  <- seurat_integrated@meta.data$repset[match(seurat_integrated_norRNA@meta.data$cells, seurat_integrated@meta.data$cells)]
# seurat_integrated_norRNA@meta.data$sample  <- seurat_integrated@meta.data$sample[match(seurat_integrated_norRNA@meta.data$cells, seurat_integrated@meta.data$cells)]


only_d11 <- subset(x = seurat_integrated, subset= (repset=="R45" & leiden_res00015_mergedclusters %in% "Ld2+3"))
M <- SetIdent(only_d11, value = "num_feeds")
DefaultAssay(M) <- "RNA"
merged_broad_repset45 <- as.SingleCellExperiment(only_d11)
merged_broad_repset45$num_feeds <- merged_broad_repset45$num_feeds=="2BF"
abundances_broad <- unclass(table(merged_broad_repset45$leiden_res00015_mergedclusters, merged_broad_repset45$sample))
colLabels(merged_broad_repset45) <- factor(merged_broad_repset45$leiden_res00015_mergedclusters)
m.out <- findMarkers(merged_broad_repset45, block=merged_broad_repset45$orig.ident, groups=merged_broad_repset45$num_feeds, assay.type="logcounts")
mdf <- as.data.frame(m.out@listData)
mdf$gene <- rownames(mdf)
up2bf_repset45 <- mdf[mdf$FALSE.FDR<0.05 & (mdf$FALSE.logFC.TRUE < (-0.25)),] %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
scran_fM_up2BF_repset45 <- up2bf_repset45[,c(7,8,10,11,12)] %>% arrange(-TRUE.logFC.FALSE)
up1bf_repset45 <- mdf[mdf$FALSE.FDR<0.05 & (mdf$FALSE.logFC.TRUE > (0.25)),] %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
scran_fM_up1BF_repset45 <- up1bf_repset45[,c(7,8,10,11,12)] %>% arrange(TRUE.logFC.FALSE)


ggg.compl <- mdf

ggg.compl$DE <- "non-significant"
#ggg_compl$avg_log2FC_base1BF <- (-1 * ggg_compl$logFC)

# if logFC > 0.25 and pvalue < 0.05, set as "UP" 
ggg.compl$DE[ggg.compl$TRUE.logFC.FALSE > (log(1.75)) & ggg.compl$TRUE.FDR < 1e-6] <- ">1.75 FC in 2BF"

ggg.compl$gene_symbol <- rownames(ggg.compl)

ggg.compl$delabel <- NA
ggg.compl$delabel[ggg.compl$DE %in% c(">1.75 FC in 1BF", ">1.75 FC in 2BF")] <- ggg.compl$gene_symbol[ggg.compl$DE %in% c(">1.75 FC in 1BF", ">1.75 FC in 2BF")]

ggplot(data=ggg.compl, aes(x=TRUE.logFC.FALSE, y=TRUE.FDR, fill=DE, label=delabel)) + 
  geom_point(shape=21, size=2.4) + 
  theme_minimal() +
  geom_text_repel(aes(col=DE), size=3.2, max.overlaps = 50) +
  scale_fill_manual(values=c("goldenrod1","seagreen3", "lightgoldenrod1", "palegreen1", "grey")) +
  scale_color_manual(values=c("goldenrod1","seagreen3","black","black","black")) +
  geom_vline(xintercept=c(-log(1.75), log(1.75)), col="black", linetype = "longdash", size=0.01) +
  geom_hline(yintercept=1e-6, col="black", linetype = "longdash", size=0.01) +
  scale_x_continuous(breaks=seq(-3,3, 0.25)) +
  scale_y_continuous(trans=scales::trans_new("sq", transform=function(x) -log10(x),
                                             inverse=function(x) 10^(x)),
                     breaks=c(0,1e-6,1e-12,1e-18,1e-24,1e-30,1e-36,1e-42,1e-48,1e-54,1e-60,1e-66,1e-72,1e-78),
                     labels=c(0,1e-6,1e-12,1e-18,1e-24,1e-30,1e-36,1e-42,1e-48,1e-54,1e-60,1e-66,1e-72,1e-78)) +
  labs(x="\nlnFC 2BF vs. 1BF", y="FDR\n") + 
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 




summed_broad_repset45 <- aggregateAcrossCells(merged_broad_repset45, use.assay.type = 1,
                                              id=colData(merged_broad_repset45)[,c("leiden_res00015_mergedclusters", "sample")])

summed.filt.broad_repset45 <- summed_broad_repset45[,summed_broad_repset45$ncells >= 10]
de.results.broad_repset45 <- pseudoBulkDGE(summed.filt.broad_repset45,
                                           label=summed.filt.broad_repset45$leiden_res00015_mergedclusters,
                                           design=~factor(orig.ident) + num_feeds,
                                           coef="num_feedsTRUE",
                                           condition=summed.filt.broad_repset45$num_feeds)

r1.broad_repset45 <- as.data.frame(de.results.broad_repset45[["Ld1"]])
r1.broad_repset45.noNA <- r1.broad_repset45[complete.cases(r1.broad_repset45), ]
r23.broad_repset45 <- as.data.frame(de.results.broad_repset45[["Ld2+3"]])
r23.broad_repset45.noNA <- r23.broad_repset45[complete.cases(r23.broad_repset45), ]

pval_to_use = 0.05
FDR_to_use = 0.05
log2FC_to_use = log2(1.25)


r1_UPREG_2BF.broad_repset45 <- r1.broad_repset45.noNA[r1.broad_repset45.noNA$PValue<pval_to_use & r1.broad_repset45.noNA$FDR<FDR_to_use & r1.broad_repset45.noNA$logFC>(log2FC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

r1_DOWNREG_2BF.broad_repset45 <- r1.broad_repset45.noNA[r1.broad_repset45.noNA$PValue<pval_to_use & r1.broad_repset45.noNA$FDR<FDR_to_use & r1.broad_repset45.noNA$logFC<(-log2FC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

r23_UPREG_2BF.broad_repset45  <- r23.broad_repset45.noNA[r23.broad_repset45.noNA$PValue<pval_to_use & r23.broad_repset45.noNA$FDR<FDR_to_use & r23.broad_repset45.noNA$logFC>(log2FC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

r23_DOWNREG_2BF.broad_repset45 <- r23.broad_repset45.noNA[r23.broad_repset45.noNA$PValue<pval_to_use & r23.broad_repset45.noNA$FDR<FDR_to_use & r23.broad_repset45.noNA$logFC<(-log2FC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))


for (i in c(1,2)){
  if (dim(list(r1_UPREG_2BF.broad_repset45,r23_UPREG_2BF.broad_repset45)[[i]])[1]==0){
    print(paste("FYI, cluster ", i, " has no significant 2bf upreg"))}}
for (i in c(1,2)){
  if (dim(list(r1_DOWNREG_2BF.broad_repset45,r23_DOWNREG_2BF.broad_repset45)[[i]])[1]==0){
    print(paste("FYI, cluster ", i, " has no significant 1bf upreg"))}}

pseubul_DE_broad_repset45_res_results_up2bf <- list(r1_UPREG_2BF.broad_repset45,r23_UPREG_2BF.broad_repset45)
names(pseubul_DE_broad_repset45_res_results_up2bf) <- c("clus1_signal","clus23_signal")

pseubul_DE_broad_repset45_res_results_up1bf <- list(r1_DOWNREG_2BF.broad_repset45,r23_DOWNREG_2BF.broad_repset45)
names(pseubul_DE_broad_repset45_res_results_up1bf) <- c("clus1_signal","clus23_signal")

bind_rows(pseubul_DE_broad_repset45_res_results_up2bf, .id="id") %>% arrange(-logFC)
bind_rows(pseubul_DE_broad_repset45_res_results_up1bf, .id="id") %>% arrange(logFC)

scran_fM_up2BF_repset45 %>% arrange(-TRUE.logFC.FALSE)
hist(scran_fM_up2BF_repset45$TRUE.logFC.FALSE, breaks=30)


filt_fMgenes1bf <- scran_fM_up1BF_repset45$gene[scran_fM_up1BF_repset45$TRUE.logFC.FALSE < (-log(4/3))]
filt_fMgenes2bf <- scran_fM_up2BF_repset45$gene[scran_fM_up2BF_repset45$TRUE.logFC.FALSE > (log(4/3))]

#########################################################################################################################################
################ pseudotime analysis section ###########################################################################################
#########################################################################################################################################

# start Fresh

rm(list=ls())

setwd("C:/Users/philipp/OneDrive - Harvard University/Desktop/single_cell_rnaseq")

library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)
library(scater)
library(scran)
library(batchelor)
library(edgeR)
library(bluster)
library(grid)
library(ggplot2)
library(slingshot)
library(monocle3)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(magrittr)
library(corrplot)
library(gplots)
library(rstatix)
library(ggforce)
library(data.table)
library(slingshot)

seurat_integrated <- readRDS("results/5reps_RPCA_9dimUMAP_400umi_125gene_only10x.rds")

table(seurat_integrated@meta.data$day_stage)

 seurat_integrated@meta.data$stage <- NA
 seurat_integrated@meta.data$exp <- "NC10x"
 seurat_integrated@meta.data$exp[grepl("mca", seurat_integrated@meta.data$cells)] <- "Howick"
# 
 seurat_integrated@meta.data$exp_numfeeds <- paste(seurat_integrated@meta.data$exp,seurat_integrated@meta.data$num_feeds)
# 
 xx <- read.table("C:/Users/philipp/Downloads/mcameta.txt", header=T)
# 
# 
seurat_integrated@meta.data$stage <- xx$stage[match(seurat_integrated@meta.data$cells, paste("mca_", xx$cell, sep=""))]
 seurat_integrated@meta.data$stage <- xx$stage[match(seurat_integrated@meta.data$cells, xx$cell)]
 seurat_integrated@meta.data$day <- xx$day[match(seurat_integrated@meta.data$cells, paste("mca_", xx$cell, sep=""))]
 seurat_integrated@meta.data$day <- xx$day[match(seurat_integrated@meta.data$cells, xx$cell)]
 seurat_integrated@meta.data$stage[grepl("^R", seurat_integrated@meta.data$cells)] <- "sgSpz"
 seurat_integrated@meta.data$day[grepl("_d11_", seurat_integrated@meta.data$cells)] <- "day11"
 seurat_integrated@meta.data$day[grepl("_d15_", seurat_integrated@meta.data$cells)] <- "day15"
# 
 seurat_integrated@meta.data$day_stage <- paste(seurat_integrated@meta.data$day, seurat_integrated@meta.data$stage, sep="_")

# Compute UMAP embedding
 seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:9, reduction = "integrated.rpca", reduction.name="umap")
x <- seurat_integrated
DefaultAssay(x) <- "RNA"
x <- JoinLayers(x)


### run this to create some prerequisite objects for pseudotemporal visualizations...

# I used to have this all wrapped up as make_hq_cell_set_with_ps_tim() function but it seemed buggy so rather run the verbose code below

cds <- as.cell_data_set(x)
cds <- cluster_cells(cds, resolution=0.00015)

# clus <- read.table("clusters_to_use.txt", header=T)
# cds@clusters$UMAP$clusters == clus$x

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE, group_label_size=10, reduction_method = "UMAP")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

# integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
# cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_principal_points = TRUE,
           group_label_size=10)

### we are rooting trajectory on leiden cluster 3

cds <- order_cells(cds,
                   #root_cells = names(cds@clusters$UMAP$clusters[cds@clusters$UMAP$clusters %in% "7"]),
                   root_pr_nodes = "Y_68",
                   reduction_method="UMAP")


# Load package: slingshot
library(slingshot)

sce.sling <- slingshot(as.SingleCellExperiment(seurat_integrated), reducedDim='UMAP')
head(sce.sling$slingPseudotime_1)
embedded <- embedCurves(sce.sling, "UMAP")
embedded <- slingCurves(embedded)[[1]] # only 1 path.
embedded <- data.frame(embedded$s[embedded$ord,])

gridExtra::grid.arrange(plotPCA(sce.sling, colour_by="leiden_res00015_clusters", point_size=1.25) +
                          theme(legend.position = "none"),
                        plotUMAP(sce.sling, colour_by="leiden_res00015_clusters", point_size=1.25) +
                          theme(legend.position = "none"),
  plotUMAP(sce.sling, colour_by="slingPseudotime_1", point_size=1.25) +
    geom_path(data=embedded, aes(x=umap_1, y=umap_2), size=0.8) + 
    geom_text(label="\u25BA", data=tail(head(embedded,2000),1), mapping=aes(x=umap_1, y=umap_2), size=4, angle=6) +
    geom_text(label="\u25BA", data=tail(head(embedded,4060),1), mapping=aes(x=umap_1, y=umap_2), size=4, angle=-50) +
    #theme(legend.position = "none") + 
    theme(legend.key.height=unit(3, "cm")) +
    #scale_color_gradient2(low = "turquoise", mid = "blue", high = "hotpink", midpoint = 11,
    scale_color_gradientn(colors=colorRampPalette(c("turquoise", "blue", "hotpink"))(20),
      breaks=c(seq(min(sce.sling$slingPseudotime_1),max(sce.sling$slingPseudotime_1, length.out=10)))))