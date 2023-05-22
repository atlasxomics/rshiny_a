library("ArchR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Rnorvegicus.UCSC.rn6")
library("data.table")
library("dplyr")
library("GenomicRanges")
library("plyr")
library("readr")  
library("Seurat")

# globals ---------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

archrproject <- args[1]

proj <- loadArchRProject(path = archrproject)

# clusters gsm -----------------------------------------------------------------

markersGS_clusters <- getMarkerFeatures( 
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest",
)

heatmapGS_cluster <- plotMarkerHeatmap(
  seMarker = markersGS_clusters, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.20",
  plotLog2FC = TRUE,
  transpose = F,
  returnMatrix = TRUE
)

saveRDS(markersGS_clusters,"/root/outs/markersGS_clusters.rds")
write.csv(heatmapGS_cluster,"/root/outs/genes_per_cluster_hm.csv")

# samples gsm ------------------------------------------------------------------

markersGS_samples <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest",
)

proj <- addImputeWeights(proj, reducedDims= "Harmony")

heatmapGS_samples <- plotMarkerHeatmap(
  seMarker = markersGS_samples, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.20",
  plotLog2FC = TRUE,
  transpose = F,
  returnMatrix = TRUE
)

saveRDS(markersGS_samples, "/root/outs/markersGS_sample.rds")
write.csv(heatmapGS_samples, "/root/outs/genes_per_sample_hm.csv")

# conditions gsm ---------------------------------------------------------------

markersGS_condition <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Condition",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest",
)

heatmapGS_condition <- plotMarkerHeatmap(
  seMarker = markersGS_condition, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.20",
  plotLog2FC = TRUE,
  transpose = F,
  returnMatrix = TRUE
)

saveRDS(markersGS_condition,"/root/outs/markersGS_condition.rds")
write.csv(heatmapGS_condition,"/root/outs/genes_per_condition_hm.csv")

# volcano plots ----------------------------------------------------------------

ncells <- length(proj$cellNames)
markerList <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Condition",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = ncells,
  normBy = "none",
  testMethod = "ttest"
)

# remove bad clusters
req_DF <- as.data.frame(getCellColData(proj))
distr_df <- table(req_DF$Clusters, req_DF$Condition)
distr <- as.data.frame.matrix(round(prop.table(as.matrix(distr_df),1),2))

lst <- list()
for(i in 1:nrow(distr)) {
  row <- distr[i,]
  if (
    sum(unname(unlist(row))>= 0.85) == 1) {
    rownames(row) -> lst[[i]]
  }
}
not_req_list <- unlist(lst)

req_clusters <- unique(proj$Clusters)
req_clusters <- req_clusters[order(as.numeric(gsub("C","",req_clusters)))]

req_clusters <- req_clusters[which(!req_clusters%in%not_req_list)]

# marker genes by condition ----------------------------------------------------

markerList_C <- list()
proj_C <- list()

for (i in seq_along(req_clusters)) {
  
  idxSample <- BiocGenerics::which(proj$Clusters == req_clusters[i])
  
  cellsSample <- proj$cellNames[idxSample]
  proj_C[i]  <- proj[cellsSample,]
  
  ncells[i] <- length(proj_C[[i]]$cellNames)
  markerList_C[[i]] <- getMarkerFeatures(
    ArchRProj = proj_C[[i]],                                                 
    useMatrix = "GeneScoreMatrix",                                                  
    groupBy = "Condition",                                                  
    bias = c("TSSEnrichment", "log10(nFrags)"),                                                 
    maxCells = ncells[[i]],                                                 
    normBy = "none",                                                  
    testMethod = "ttest"                                                  
  )                                                 
}
names(markerList_C) <- req_clusters

# lets find empty genes
gsm <- getMatrixFromProject(proj)
gsm_mat <- assay(getMatrixFromProject(proj),"GeneScoreMatrix")
empty_gene_idx <- which(rowSums((gsm_mat))==0)
empty_gene <- rowData(gsm)$name[empty_gene_idx]

# volcano data Conditions ------------------------------------------------------

markerList_df1 <- assay(markerList, "Log2FC")
markerList_df2 <- assay(markerList, "Pval")
markerList_df3 <- assay(markerList, "FDR")
markerList_df <- cbind(markerList_df1, markerList_df2, markerList_df3)
markerList_df$genes<- rowData(markerList)$name
markerList_df$cluster <- rep("All",length(rownames(markerList_df)))

markerList_df <- markerList_df[,c(1,3,5,7,8)]                                                     
colnames(markerList_df) <- c(                                                     
  "avg_log2FC",                                                     
  "p_val",                                                      
  "p_val_adj",                                                      
  "gene",                                                     
  "cluster"                                                     
)  

# volcano data Conditions ------------------------------------------------------

markerList_df1_C <- list()
markerList_df2_C <- list()
markerList_df3_C <- list()
markerList_df_C  <- list()

for (i in seq_along(req_clusters)){
  
  cluster <- req_clusters[i]
  
  markerList_df1_C[[i]] <- assay(markerList_C[[i]], "Log2FC")
  markerList_df2_C[[i]] <- assay(markerList_C[[i]], "Pval")
  markerList_df3_C[[i]] <- assay(markerList_C[[i]], "FDR")
  markerList_df_C[[i]] <- cbind(
    markerList_df1_C[[i]],
    markerList_df2_C[[i]],
    markerList_df3_C[[i]]
  )
  markerList_df_C[[i]]$genes<- rowData(markerList_C[[i]])$name
  markerList_df_C[[i]]$cluster <- rep(
    cluster, 
    length(rownames(markerList_df_C[[i]]))
  )

  markerList_df_C[[i]] <- markerList_df_C[[i]][,c(1,3,5,7,8)]
  
  colnames(markerList_df_C[[i]]) <- c(
    "avg_log2FC",
    "p_val",
    "p_val_adj",
    "gene",
    "cluster"
  )
}
names(markerList_df_C) <- req_clusters

# merge all data frames
markersGS_merged_df <- do.call("rbind", markerList_df_C)

# also data frame for all clusters together needs to be added
markersGS_merged_df <- rbind(markerList_df, markersGS_merged_df)

# remove empty genes
markersGS_merged_df <- markersGS_merged_df[which(!markersGS_merged_df$gene%in%empty_gene),]

# remove na values
markersGS_merged_df <- na.omit(markersGS_merged_df)

# remove FDR equal to 0
markersGS_merged_df <- markersGS_merged_df[which(!markersGS_merged_df$p_val_adj== 0),]

# make logfc limiation between 1 and -1
markersGS_merged_df <- markersGS_merged_df[which(abs(markersGS_merged_df$avg_log2FC)< 1.2),]

markersGS_merged_df$Significance <- ifelse(
  markersGS_merged_df$p_val_adj < 10^-1,
  ifelse(
    markersGS_merged_df$avg_log2FC > 0.0,
    colnames(markerList)[1],
    colnames(markerList)[2]),
  'Not siginficant'
  )

de <- markersGS_merged_df
write.table(
  de,
  "/root/outs/inpMarkers.txt",
  sep = '\t',
  quote = F,
  row.names = F
)
inpMarkers = fread("/root/outs//inpMarkers.txt") 

# extract UMAP -----------------------------------------------------------------

UMAPHarmony <-getEmbedding(
  ArchRProj = proj,
  embedding = "UMAP",
  returnDF = TRUE
)

write.csv(UMAPHarmony,"/root/outs/UMAPHarmony.csv")

# heatmap genes by cluster -----------------------------------------------------

hm_per_clust <- read.csv("/root/outs/genes_per_cluster_hm.csv")

nClust <- length(unique(proj@cellColData@listData[["Clusters"]]))

df = list()
for (i in seq_along(1:nClust)){
  df[[i]] <- hm_per_clust[,c(1,i+1)]
  
  #select top 5 values by group
  df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:5,1]
  
}
final <- do.call(rbind, df)

req_genes1 <- unlist(df)
req_genes1 <- req_genes1[!duplicated(req_genes1)]

write.csv(req_genes1,"/root/outs/req_genes1.csv")

# heatmap genes by condition ---------------------------------------------------

hm_per_cond <- read.csv("/root/outs/genes_per_condition_hm.csv")

nConds <- length(unique(proj@cellColData@listData[["Condition"]]))

df = list()
for (i in seq_along(1:nConds)){
  df[[i]] <- hm_per_cond[,c(1,i+1)]
  
  # select top 20 values by condition
  df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:20,1]
}

final <- do.call(rbind, df)

req_genes2 <- unlist(df)
req_genes2 <- req_genes2[!duplicated(req_genes2)]

write.csv(req_genes2,"/root/outs/req_genes2.csv")

# heatmap genes by sample -- ---------------------------------------------------

hm_per_sample <- read.csv("/root/outs/genes_per_sample_hm.csv")

nSamples = length(proj@sampleColData@rownames)

df = list()
for (i in seq_along(1:nSamples)){
  df[[i]] <- hm_per_sample[,c(1,i+1)]
  
  # select top 10 values by condition
  df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:10,1]
}

final <- do.call(rbind, df)

req_genes3 <- unlist(df)
req_genes3 <- req_genes3[!duplicated(req_genes3)]

write.csv(req_genes3,"/root/outs/req_genes3.csv")
