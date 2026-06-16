
get_conserved <- function(clus_num=clus_num, groupvar=groupvar){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = clus_num,
                       grouping.var = groupvar,
                       only.pos = TRUE,
                       logfc.threshold = 0.25) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = clus_num, .)
}


get_wilcox_DE_upingroup2 <- function(replicates_to_use=replicates_to_use, dpbm_to_use=dpbm_to_use, adj_pval_to_use=adj_pval_to_use,
                                  log2FC_to_use=log2FC_to_use, groupingvar=groupingvar, group1=group1, group2=group2, num_feeds_to_use=num_feeds_to_use){
  only_d11_R4_R5 <- subset(x = seurat_integrated, subset= ((orig.ident %in% replicates_to_use) & dpbm==dpbm_to_use & num_feeds %in% num_feeds_to_use))
  M_only_d11_R4_R5 <- SetIdent(only_d11_R4_R5, value = groupingvar)
  DefaultAssay(M_only_d11_R4_R5) <- "RNA"
  FindMarkers(M_only_d11_R4_R5, ident.1 = group1, ident.2 = group2, assay="RNA") %>% filter(p_val_adj<adj_pval_to_use & avg_log2FC<(-log2FC_to_use)) %>%  rownames_to_column("gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
}


get_wilcox_DE_upingroup1 <- function(replicates_to_use=replicates_to_use, dpbm_to_use=dpbm_to_use, adj_pval_to_use=adj_pval_to_use,
                                  log2FC_to_use=log2FC_to_use, groupingvar=groupingvar, group1=group1, group2=group2, num_feeds_to_use=num_feeds_to_use){
  only_d11_R4_R5 <- subset(x = seurat_integrated, subset= ((orig.ident %in% replicates_to_use) & dpbm==dpbm_to_use & num_feeds %in% num_feeds_to_use))
  M_only_d11_R4_R5 <- SetIdent(only_d11_R4_R5, value = groupingvar)
  DefaultAssay(M_only_d11_R4_R5) <- "RNA"
  FindMarkers(M_only_d11_R4_R5, ident.1 = group1, ident.2 = group2) %>% filter(p_val_adj<adj_pval_to_use & avg_log2FC>log2FC_to_use) %>%  rownames_to_column("gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
}


get_pseubul_DE_broad_res <- function(FDR_to_use=FDR_to_use, pval_to_use=pval_to_use, lnFC_to_use=lnFC_to_use){
only_d11 <- subset(x = seurat_integrated)
M <- SetIdent(only_d11, value = "num_feeds")
DefaultAssay(M) <- "RNA"
merged_broad <- as.SingleCellExperiment(only_d11)
merged_broad$num_feeds <- merged_broad$num_feeds=="2BF"
abundances_broad <- unclass(table(merged_broad$leiden_res00015_clusters, merged_broad$sample))
colLabels(merged_broad) <- factor(merged_broad$leiden_res00015_clusters)

merged_broad <<- merged_broad
summed_broad <- aggregateAcrossCells(merged_broad, use.assay.type = 1,
                               id=colData(merged_broad)[,c("leiden_res00015_clusters", "sample")])

summed.filt.broad <<- summed_broad[,summed_broad$ncells >= 10]
de.results.broad <- pseudoBulkDGE(summed.filt.broad, 
                            label=summed.filt.broad$leiden_res00015_clusters,
                            design=~factor(orig.ident) + num_feeds,
                            coef="num_feedsTRUE",
                            condition=summed.filt.broad$num_feeds)


r1.broad <- as.data.frame(de.results.broad[["Ld1"]])
r1.broad.noNA <- r1.broad[complete.cases(r1.broad), ]
r2.broad <- as.data.frame(de.results.broad[["Ld2"]])
r2.broad.noNA <- r2.broad[complete.cases(r2.broad), ]
r3.broad <- as.data.frame(de.results.broad[["Ld3"]])
r3.broad.noNA <- r3.broad[complete.cases(r3.broad), ]


r1_UPREG_2BF.broad <<- r1.broad.noNA[r1.broad.noNA$PValue<pval_to_use & r1.broad.noNA$FDR<FDR_to_use & r1.broad.noNA$logFC>(lnFC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

r1_DOWNREG_2BF.broad <<- r1.broad.noNA[r1.broad.noNA$PValue<pval_to_use & r1.broad.noNA$FDR<FDR_to_use & r1.broad.noNA$logFC<(-lnFC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

r2_UPREG_2BF.broad <<- r2.broad.noNA[r2.broad.noNA$PValue<pval_to_use & r2.broad.noNA$FDR<FDR_to_use & r2.broad.noNA$logFC>(lnFC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

r2_DOWNREG_2BF.broad <<- r2.broad.noNA[r2.broad.noNA$PValue<pval_to_use & r2.broad.noNA$FDR<FDR_to_use & r2.broad.noNA$logFC<(-lnFC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

r3_UPREG_2BF.broad <<- r3.broad.noNA[r3.broad.noNA$PValue<pval_to_use & r3.broad.noNA$FDR<FDR_to_use & r3.broad.noNA$logFC>(lnFC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

r3_DOWNREG_2BF.broad <<- r3.broad.noNA[r3.broad.noNA$PValue<pval_to_use & r3.broad.noNA$FDR<FDR_to_use & r3.broad.noNA$logFC<(-lnFC_to_use),] %>%  rownames_to_column("gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))

for (i in c(1,2,3)){
  if (dim(list(r1_UPREG_2BF.broad,r2_UPREG_2BF.broad,r3_UPREG_2BF.broad)[[i]])[1]==0){
    print(paste("FYI, cluster ", i, " has no significant 2bf upreg"))}}

pseubul_DE_broad_res_results_up2bf <- list(r1_UPREG_2BF.broad,r2_UPREG_2BF.broad,r3_UPREG_2BF.broad)
names(pseubul_DE_broad_res_results_up2bf) <- c("clus1_signal","clus2_signal","clus3_signal")

pseubul_DE_broad_res_results_up2bf <<- pseubul_DE_broad_res_results_up2bf

for (i in c(1,2,3)){
  if (dim(list(r1_DOWNREG_2BF.broad,r2_DOWNREG_2BF.broad,r3_DOWNREG_2BF.broad)[[i]])[1]==0){
    print(paste("FYI, cluster ", i, " has no significant 1bf upreg"))}}

pseubul_DE_broad_res_results_up1bf <- list(r1_DOWNREG_2BF.broad,r2_DOWNREG_2BF.broad,r3_DOWNREG_2BF.broad)
names(pseubul_DE_broad_res_results_up1bf) <- c("clus1_signal","clus2_signal","clus3_signal")

pseubul_DE_broad_res_results_up1bf <<- pseubul_DE_broad_res_results_up1bf

}


#################

get_pseubul_DE_hi_res <- function(FDR_to_use=FDR_to_use, pval_to_use=pval_to_use, lnFC_to_use=lnFC_to_use){
  only_d11 <- subset(x = seurat_integrated)
  M <- SetIdent(only_d11, value = "num_feeds")
  DefaultAssay(M) <- "RNA"
  merged_hi <- as.SingleCellExperiment(only_d11)
  merged_hi$num_feeds <- merged_hi$num_feeds=="2BF"
  abundances_hi <- unclass(table(merged_hi$leiden_res00015_clusters, merged_hi$sample))
  colLabels(merged_hi) <- factor(merged_hi$leiden_res00015_clusters)
  
  merged_hi <<- merged_hi
  
  table(Cluster=colLabels(merged_hi), Sample=merged_hi$sample)
  summed_hi <- aggregateAcrossCells(merged_hi, use.assay.type = 1,
                                       id=colData(merged_hi)[,c("leiden_res00015_clusters", "sample")])
  summed.filt.hi <<- summed_hi[,summed_hi$ncells >= 10]
  de.results.hi <- pseudoBulkDGE(summed.filt.hi, 
                                    label=summed.filt.hi$leiden_res00015_clusters,
                                    design=~factor(orig.ident) + num_feeds,
                                    coef="num_feedsTRUE",
                                    condition=summed.filt.hi$num_feeds)
  
  r0.hi <- as.data.frame(de.results.hi[["0"]])
  r0.hi.noNA <- r0.hi[complete.cases(r0.hi), ]
  r1.hi <- as.data.frame(de.results.hi[["1"]])
  r1.hi.noNA <- r1.hi[complete.cases(r1.hi), ]
  r2.hi <- as.data.frame(de.results.hi[["2"]])
  r2.hi.noNA <- r2.hi[complete.cases(r2.hi), ]
  r3.hi <- as.data.frame(de.results.hi[["3"]])
  r3.hi.noNA <- r3.hi[complete.cases(r3.hi), ]
  r4.hi <- as.data.frame(de.results.hi[["4"]])
  r4.hi.noNA <- r4.hi[complete.cases(r4.hi), ]
  r5.hi <- as.data.frame(de.results.hi[["5"]])
  r5.hi.noNA <- r5.hi[complete.cases(r5.hi), ]
  
  r0_UPREG_2BF.hi <<- r0.hi.noNA[r0.hi.noNA$PValue<pval_to_use & r0.hi.noNA$FDR<FDR_to_use & r0.hi.noNA$logFC>(lnFC_to_use),] %>%  rownames_to_column("gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
  
  r0_DOWNREG_2BF.hi <<- r0.hi.noNA[r0.hi.noNA$PValue<pval_to_use & r0.hi.noNA$FDR<FDR_to_use & r0.hi.noNA$logFC<(-lnFC_to_use),] %>%  rownames_to_column("gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
  
  r1_UPREG_2BF.hi <<- r1.hi.noNA[r1.hi.noNA$PValue<pval_to_use & r1.hi.noNA$FDR<FDR_to_use & r1.hi.noNA$logFC>(lnFC_to_use),] %>%  rownames_to_column("gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
  
  r1_DOWNREG_2BF.hi <<- r1.hi.noNA[r1.hi.noNA$PValue<pval_to_use & r1.hi.noNA$FDR<FDR_to_use & r1.hi.noNA$logFC<(-lnFC_to_use),] %>%  rownames_to_column("gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
  
  r2_UPREG_2BF.hi <<- r2.hi.noNA[r2.hi.noNA$PValue<pval_to_use & r2.hi.noNA$FDR<FDR_to_use & r2.hi.noNA$logFC>(lnFC_to_use),] %>%  rownames_to_column("gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
  
  r2_DOWNREG_2BF.hi <<- r2.hi.noNA[r2.hi.noNA$PValue<pval_to_use & r2.hi.noNA$FDR<FDR_to_use & r2.hi.noNA$logFC<(-lnFC_to_use),] %>%  rownames_to_column("gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
  
  r3_UPREG_2BF.hi <<- r3.hi.noNA[r3.hi.noNA$PValue<pval_to_use & r3.hi.noNA$FDR<FDR_to_use & r3.hi.noNA$logFC>(lnFC_to_use),] %>%  rownames_to_column("gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
  
  r3_DOWNREG_2BF.hi <<- r3.hi.noNA[r3.hi.noNA$PValue<pval_to_use & r3.hi.noNA$FDR<FDR_to_use & r3.hi.noNA$logFC<(-lnFC_to_use),] %>%  rownames_to_column("gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
  
  r4_UPREG_2BF.hi <<- r4.hi.noNA[r4.hi.noNA$PValue<pval_to_use & r4.hi.noNA$FDR<FDR_to_use & r4.hi.noNA$logFC>(lnFC_to_use),] %>%  rownames_to_column("gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
  
  r4_DOWNREG_2BF.hi <<- r4.hi.noNA[r4.hi.noNA$PValue<pval_to_use & r4.hi.noNA$FDR<FDR_to_use & r4.hi.noNA$logFC<(-lnFC_to_use),] %>%  rownames_to_column("gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
  
  r5_UPREG_2BF.hi <<- r5.hi.noNA[r5.hi.noNA$PValue<pval_to_use & r5.hi.noNA$FDR<FDR_to_use & r5.hi.noNA$logFC>(lnFC_to_use),] %>%  rownames_to_column("gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
  
  r5_DOWNREG_2BF.hi <<- r5.hi.noNA[r5.hi.noNA$PValue<pval_to_use & r5.hi.noNA$FDR<FDR_to_use & r5.hi.noNA$logFC<(-lnFC_to_use),] %>%  rownames_to_column("gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
  
  
  for (i in c(1,2,3,4,5,6)){
    if (dim(list(r0_UPREG_2BF.hi,r1_UPREG_2BF.hi,r2_UPREG_2BF.hi,r3_UPREG_2BF.hi,r4_UPREG_2BF.hi,r5_UPREG_2BF.hi)[[i]])[1]==0){
      print(paste("FYI, cluster ", i-1, " has no significant signals between treatments"))}}
  
  pseubul_DE_hi_res_results_up2bf <- list(r0_UPREG_2BF.hi,r1_UPREG_2BF.hi,r2_UPREG_2BF.hi,r3_UPREG_2BF.hi,r4_UPREG_2BF.hi,r5_UPREG_2BF.hi)
  names(pseubul_DE_hi_res_results_up2bf) <- c("clus0_signal","clus1_signal","clus2_signal","clus3_signal","clus4_signal","clus5_signal")
  
  pseubul_DE_hi_res_results_up2bf <<- pseubul_DE_hi_res_results_up2bf
  
  for (i in c(1,2,3,4,5,6)){
    if (dim(list(r0_DOWNREG_2BF.hi,r1_DOWNREG_2BF.hi,r2_DOWNREG_2BF.hi,r3_DOWNREG_2BF.hi,r4_DOWNREG_2BF.hi,r5_DOWNREG_2BF.hi)[[i]])[1]==0){
      print(paste("FYI, cluster ", i-1, " has no significant signals between treatments"))}}
  
  pseubul_DE_hi_res_results_up1bf <- list(r0_DOWNREG_2BF.hi,r1_DOWNREG_2BF.hi,r2_DOWNREG_2BF.hi,r3_DOWNREG_2BF.hi,r4_DOWNREG_2BF.hi,r5_DOWNREG_2BF.hi)
  names(pseubul_DE_hi_res_results_up1bf) <- c("clus0_signal","clus1_signal","clus2_signal","clus3_signal","clus4_signal","clus5_signal")
  
  pseubul_DE_hi_res_results_up1bf <<- pseubul_DE_hi_res_results_up1bf
  
}



make_PDF_pseubul_DE_broad_res <- function(filename=filename){
  
  rownames(r1_UPREG_2BF.broad) <- r1_UPREG_2BF.broad$gene
  rownames(r2_UPREG_2BF.broad) <- r2_UPREG_2BF.broad$gene
  rownames(r3_UPREG_2BF.broad) <- r3_UPREG_2BF.broad$gene
  
  UPREG_2BF.broad_vector <- c(rownames(r2_UPREG_2BF.broad[order(-r2_UPREG_2BF.broad$logFC),]),
                              rownames(r3_UPREG_2BF.broad[order(-r3_UPREG_2BF.broad$logFC),]),
                              rownames(r1_UPREG_2BF.broad[order(-r1_UPREG_2BF.broad$logFC),]))
  
  UPREG_2BF.broad_multiple_occurrences <- UPREG_2BF.broad_vector[duplicated(UPREG_2BF.broad_vector)]
  UPREG_2BF.broad_single_occurrences <- UPREG_2BF.broad_vector[!duplicated(UPREG_2BF.broad_vector)]
  
  mygenes <- c(unique(UPREG_2BF.broad_vector[UPREG_2BF.broad_vector %in% UPREG_2BF.broad_multiple_occurrences]),
               UPREG_2BF.broad_vector[UPREG_2BF.broad_vector %in% UPREG_2BF.broad_single_occurrences])
  
  summed.filt.broad$leiden_res00015_clusters_F <- factor(summed.filt.broad$leiden_res00015_clusters, levels=c("Ld1","Ld2","Ld3"))
  pdf(file=filename, height=3, width=7)
  for (i in mygenes){
    myplot <- plotExpression(logNormCounts(summed.filt.broad), features=i, add_legend=T, point_size=5, show_violin=T, show_median=F,
                             x="num_feeds", colour_by="sample",
                             xlab=paste("", "2 blood feeds?", "", "logFC", gsub("\\=  ","\\= NS ", paste("Ld1 = ", r2_UPREG_2BF.broad[rownames(r1_UPREG_2BF.broad)==i,]$logFC, " ", 
                                                                                                             "Ld2 = ", r3_UPREG_2BF.broad[rownames(r2_UPREG_2BF.broad)==i,]$logFC, " ", 
                                                                                                             "Ld3 = ", r1_UPREG_2BF.broad[rownames(r3_UPREG_2BF.broad)==i,]$logFC, " ", sep="")), "",
                                        "FDR", gsub("\\=  ","\\= NS ", paste("Ld1 = ", r2_UPREG_2BF.broad[rownames(r1_UPREG_2BF.broad)==i,]$FDR, " ", 
                                                                                 "Ld2 = ", r3_UPREG_2BF.broad[rownames(r2_UPREG_2BF.broad)==i,]$FDR, " ", 
                                                                                 "Ld3 = ", r1_UPREG_2BF.broad[rownames(r3_UPREG_2BF.broad)==i,]$FDR, " ", sep="")), "", i, sep="\n"),
                             other_fields="leiden_res00015_clusters_F") + facet_grid(~leiden_res00015_clusters_F) + theme(axis.title.x = element_text(size = 5)) +
      scale_color_manual(values=c("lightgreen","lightgreen","grey80","grey80","violet","violet","pink","pink","lightblue","lightblue")) +
      theme(legend.key.size = unit(0.01, 'cm')) + guides(color = guide_legend(override.aes = list(size = 3)))
    
    
    
    
    g <- ggplot_gtable(ggplot_build(myplot))
    stripr <- which(grepl('strip-t', g$layout$name))
    fills <- paste(c(sum(as.integer(rownames(r1_UPREG_2BF.broad)==i)),
                     sum(as.integer(rownames(r2_UPREG_2BF.broad)==i)),
                     sum(as.integer(rownames(r3_UPREG_2BF.broad)==i))), collapse=" ")
    fills <- gsub(0, "grey", fills)
    fills <- gsub(1, "yellow", fills)
    fills2 <- strsplit(fills, " ")[[1]]
    k <- 1
    for (i in stripr) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills2[k]
      k <- k+1
    }
    
    
    grid::grid.draw(g)
    grid.newpage()
  }
  
  dev.off()
  
}



make_PDF_pseubul_DE_hi_res <- function(filename=filename){
  
  rownames(r0_UPREG_2BF.hi) <- r0_UPREG_2BF.hi$gene
  rownames(r1_UPREG_2BF.hi) <- r1_UPREG_2BF.hi$gene
  rownames(r2_UPREG_2BF.hi) <- r2_UPREG_2BF.hi$gene
  rownames(r3_UPREG_2BF.hi) <- r3_UPREG_2BF.hi$gene
  rownames(r4_UPREG_2BF.hi) <- r4_UPREG_2BF.hi$gene
  rownames(r5_UPREG_2BF.hi) <- r5_UPREG_2BF.hi$gene
  
  UPREG_2BF.hi_vector <- c(rownames(r5_UPREG_2BF.hi[order(-r5_UPREG_2BF.hi$logFC),]),
                              rownames(r0_UPREG_2BF.hi[order(-r5_UPREG_2BF.hi$logFC),]),
                              rownames(r1_UPREG_2BF.hi[order(-r1_UPREG_2BF.hi$logFC),]),
                              rownames(r3_UPREG_2BF.hi[order(-r3_UPREG_2BF.hi$logFC),]),
                              rownames(r4_UPREG_2BF.hi[order(-r4_UPREG_2BF.hi$logFC),]),
                              rownames(r2_UPREG_2BF.hi[order(-r2_UPREG_2BF.hi$logFC),]))
  
  UPREG_2BF.hi_multiple_occurrences <- UPREG_2BF.hi_vector[duplicated(UPREG_2BF.hi_vector)]
  UPREG_2BF.hi_single_occurrences <- UPREG_2BF.hi_vector[!duplicated(UPREG_2BF.hi_vector)]
  
  mygenes <- c(unique(UPREG_2BF.hi_vector[UPREG_2BF.hi_vector %in% UPREG_2BF.hi_multiple_occurrences]),
               UPREG_2BF.hi_vector[UPREG_2BF.hi_vector %in% UPREG_2BF.hi_single_occurrences])
  
  summed.filt.hi$leiden_res00015_clusters_F <- factor(summed.filt.hi$leiden_res00015_clusters, levels=c("2","4","3","1","0","5"))
  pdf(file=filename, height=3, width=7)
  for (i in mygenes){
    myplot <- plotExpression(logNormCounts(summed.filt.hi), features=i, add_legend=T, point_size=5, show_violin=T, show_median=F,
                             x="num_feeds", colour_by="sample",
                             xlab=paste("", "2 blood feeds?", "", "logFC", gsub("\\= \\,","\\= NS\\,", paste("c2 = ", r2_UPREG_2BF.hi[rownames(r2_UPREG_2BF.hi)==i,]$logFC, ", ", 
                                                                                                             "c4 = ", r4_UPREG_2BF.hi[rownames(r4_UPREG_2BF.hi)==i,]$logFC, ", ", 
                                                                                                             "c3 = ", r3_UPREG_2BF.hi[rownames(r3_UPREG_2BF.hi)==i,]$logFC, ", ",
                                                                                                             "c1 = ", r1_UPREG_2BF.hi[rownames(r1_UPREG_2BF.hi)==i,]$logFC, ", ", 
                                                                                                             "c0 = ", r0_UPREG_2BF.hi[rownames(r0_UPREG_2BF.hi)==i,]$logFC, ", ", 
                                                                                                             "c5 = ", r5_UPREG_2BF.hi[rownames(r5_UPREG_2BF.hi)==i,]$logFC, sep="")), "",
                                        "FDR", gsub("\\= \\,","\\= NS\\,", paste("c2 = ", r2_UPREG_2BF.hi[rownames(r2_UPREG_2BF.hi)==i,]$FDR, ", ", 
                                                                                 "c4 = ", r4_UPREG_2BF.hi[rownames(r4_UPREG_2BF.hi)==i,]$FDR, ", ", 
                                                                                 "c3 = ", r3_UPREG_2BF.hi[rownames(r3_UPREG_2BF.hi)==i,]$FDR, ", ",
                                                                                 "c1 = ", r1_UPREG_2BF.hi[rownames(r1_UPREG_2BF.hi)==i,]$FDR, ", ",
                                                                                 "c0 = ", r0_UPREG_2BF.hi[rownames(r0_UPREG_2BF.hi)==i,]$FDR, ", ",
                                                                                 "c5 = ", r5_UPREG_2BF.hi[rownames(r5_UPREG_2BF.hi)==i,]$FDR, sep="")), "", i, sep="\n"),
                             other_fields="leiden_res00015_clusters_F") + facet_grid(~leiden_res00015_clusters_F) + theme(axis.title.x = element_text(size = 5)) +
      scale_color_manual(values=c("lightgreen","lightgreen","grey80","grey80","violet","violet","pink","pink","lightblue","lightblue")) +
      theme(legend.key.size = unit(0.01, 'cm')) + guides(color = guide_legend(override.aes = list(size = 3)))
    
    
    
    
    g <- ggplot_gtable(ggplot_build(myplot))
    stripr <- which(grepl('strip-t', g$layout$name))
    fills <- paste(c(sum(as.integer(rownames(r2_UPREG_2BF.hi)==i)),
                     sum(as.integer(rownames(r4_UPREG_2BF.hi)==i)),
                     sum(as.integer(rownames(r3_UPREG_2BF.hi)==i)),
                     sum(as.integer(rownames(r1_UPREG_2BF.hi)==i)),
                     sum(as.integer(rownames(r0_UPREG_2BF.hi)==i)),
                     sum(as.integer(rownames(r5_UPREG_2BF.hi)==i))), collapse=" ")
    fills <- gsub(0, "grey", fills)
    fills <- gsub(1, "yellow", fills)
    fills2 <- strsplit(fills, " ")[[1]]
    k <- 1
    for (i in stripr) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills2[k]
      k <- k+1
    }
    
    
    grid::grid.draw(g)
    grid.newpage()
  }
  
  dev.off()
  
}
