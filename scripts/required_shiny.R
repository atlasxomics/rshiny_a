library("ArchR")
library('BSgenome')
library("BSgenome.Mmusculus.UCSC.mm10")
library("circlize")
library("chromVARmotifs")
library("ComplexHeatmap")
library("data.table")
library("dplyr")
library("GenomicRanges")
library("ggrepel")
library("ggseqlogo")
library("plyr")
library("readr")  
library("Seurat")
library("seqLogo")
library("qdap")


setwd('~/latch/rshiny_a/scripts/')

proj3 <- loadArchRProject(path = "./craft-test_ArchRProject/")

proj3 <- addHarmony(
  ArchRProj = proj3,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

proj3 <- addClusters(
  input = proj3,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters",
  resolution = c(1), 
  force = TRUE
)

proj3 <- addUMAP(
  ArchRProj = proj3,
  reducedDims = "Harmony",
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.0,  
  metric = "cosine",
  force = TRUE
)

# clusters gsm -----------------------------------------------------------------

markersGS_clusters <- getMarkerFeatures( # what if only one cluster? error?
  ArchRProj = proj3, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest"
)

heatmapGS_cluster <- plotMarkerHeatmap(
  seMarker = markersGS_clusters, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.20",
  plotLog2FC = TRUE,
  transpose = F,
  returnMatrix = TRUE
)

saveRDS(markersGS_clusters,"shiny2/markersGS_clusters.rds")
write.csv(heatmapGS_cluster,"shiny2/genes_per_cluster_hm.csv")

# samples gsm ------------------------------------------------------------------

markersGS_samples <- getMarkerFeatures(
  ArchRProj = proj3, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest",
)

proj3 <- addImputeWeights(proj3, reducedDims= "Harmony")

heatmapGS_samples <- plotMarkerHeatmap(
  seMarker = markersGS_samples, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.20",
  plotLog2FC = TRUE,
  transpose = F,
  returnMatrix = TRUE
)

saveRDS(markersGS_samples, "shiny2/markersGS_sample.rds")
write.csv(heatmapGS_samples, "shiny2/genes_per_sample_hm.csv")

# conditions gsm ---------------------------------------------------------------

markersGS_condition <- getMarkerFeatures(
  ArchRProj = proj3, 
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

saveRDS(markersGS_condition,"shiny2/markersGS_condition.rds")
write.csv(heatmapGS_condition,"shiny2/genes_per_condition_hm.csv")

# volcano plots ----------------------------------------------------------------

req_DF <- as.data.frame(getCellColData(proj3))
ncells <- length(proj3$cellNames)
markerList <- getMarkerFeatures(
  ArchRProj = proj3,
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Condition",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = ncells,
  normBy = "none",
  testMethod = "ttest")

# remove clusters with 
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

req_clusters <- unique(proj3$Clusters)
req_clusters <- req_clusters[order(as.numeric(gsub("C","",req_clusters)))]

req_clusters <- req_clusters[which(!req_clusters%in%not_req_list)]

# marker genes by condition ----------------------------------------------------

markerList_C <- list()
proj3_C <- list()

for (i in seq_along(req_clusters)) {
  
  idxSample <- BiocGenerics::which(proj3$Clusters == req_clusters[i])
  
  cellsSample <- proj3$cellNames[idxSample]
  proj3_C[i]  <- proj3[cellsSample,]
  
  ncells[i] <- length(proj3_C[[i]]$cellNames)
  
  # this needs to be in a trycatch                #
  markerList_C[[i]] <- getMarkerFeatures(         #                                        
    ArchRProj = proj3_C[[i]],                     #                            
    useMatrix = "GeneScoreMatrix",                #                                  
    groupBy = "Condition",                        #                          
    bias = c("TSSEnrichment", "log10(nFrags)"),   #                                              
    maxCells = ncells[[i]],                       #                          
    normBy = "none",                              #                    
    testMethod = "ttest"                          #                        
  )                                               #  
}
names(markerList_C) <- req_clusters

# lets find empty genes
gsm <- getMatrixFromProject(proj3)
gsm_mat <- assay(getMatrixFromProject(proj3),"GeneScoreMatrix")
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

markerList_df

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
  markerList_df_C[[i]]$cluster <- rep(cluster, 
                                      length(rownames(markerList_df_C[[i]])))
    
  # What if conditions == 3??
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
markersGS_merged_df <- markersGS_merged_df[
  which(!markersGS_merged_df$gene%in%empty_gene),
]

# remove na values
markersGS_merged_df <- na.omit(markersGS_merged_df)

# remove FDR equal to 0
markersGS_merged_df <- markersGS_merged_df[
  which(!markersGS_merged_df$p_val_adj == 0),
]

# make logfc limiation between 1 and -1
markersGS_merged_df <- markersGS_merged_df[
  which(abs(markersGS_merged_df$avg_log2FC)< 1.2),
]

markersGS_merged_df$Significance = ifelse(
  markersGS_merged_df$p_val_adj < 10^-1 ,
  ifelse(
    markersGS_merged_df$avg_log2FC > 0.0 ,
    colnames(markerList)[1],
    colnames(markerList)[2]),
  'Not siginficant'
)

de <- markersGS_merged_df
write.table(
  de,
  "shiny2/inpMarkers.txt",
  sep = '\t',
  quote = F,
  row.names = F
)
inpMarkers = fread("shiny2/inpMarkers.txt")

# plot volcano -----------------------------------------------------------------

targeted <- 'C4'

minfdr = minfdr = 0.05
minfdr2 = 10^-(2/3 *(-log10(min(inpMarkers[cluster == targeted]$p_val_adj))))

# How to correctly assign conditions?                                      # 
ggData = inpMarkers[cluster == targeted]                                   #       
ggData$Significance = ifelse(ggData$p_val_adj < minfdr,                    #                     
                             ifelse(ggData$avg_log2FC > 0.0,"WT","Lupus"), #                                         
                             'Not siginficant')                            #             

ggData$Significance <- factor(
  ggData$Significance,
  levels = c('WT','Lupus','Not siginficant')
)


ggData[p_val_adj < 1e-300]$p_val_adj = 1e-300
ggData$log10fdr = -log10(ggData$p_val_adj)
ggplot(ggData, aes(avg_log2FC, log10fdr)) + 
  geom_point() + 
  ylab("-log10(FDR)") +
  geom_point(aes(color = Significance)) +
  scale_color_manual(values = c("red","blue","grey"))  +
  geom_text_repel(data = subset(ggData, p_val_adj < 10^-9 ), aes(label = gene))

# extract UMAP -----------------------------------------------------------------

UMAPHarmony <-getEmbedding(
  ArchRProj = proj3,
  embedding = "UMAP",
  returnDF = TRUE
)

write.csv(UMAPHarmony,"shiny2/UMAPHarmony.csv")

# motifs -----------------------------------------------------------------------

proj3 <- addGroupCoverages(ArchRProj = proj3, groupBy = "Clusters")

pathToMacs2 <- findMacs2()
proj3 <- addReproduciblePeakSet(
  ArchRProj = proj3, 
  groupBy = "Clusters", 
  pathToMacs2 = pathToMacs2,
  genomeSize = 3.0e+09, # put in dict
  force = TRUE
)

proj3 <- addMotifAnnotations(
  ArchRProj = proj3,
  motifSet = "cisbp",
  name = "Motif",
  force = TRUE
)

PWMs <- getPeakAnnotation(proj3, "Motif")$motifs
PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
}

ProbMatrices <- lapply(PWMs, PWMatrixToProbMatrix)
lapply(ProbMatrices, colSums) %>% range

saveRDS(ProbMatrices,"shiny2/seqlogo.rds")
        
# heatmap genes by cluster -----------------------------------------------------

hm_per_clust <- read.csv("shiny2/genes_per_cluster_hm.csv")

nClust <- length(unique(proj3@cellColData@listData[["Clusters"]]))

df = list()
for (i in seq_along(1:nClust)){
  df[[i]] <- hm_per_clust[,c(1,i+1)]
  
  #select top 5 values by group
  df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:5,1]
  
}
final <- do.call(rbind, df)

req_genes1 <- unlist(df)
req_genes1<- req_genes1[!duplicated(req_genes1)]

write.csv(req_genes1,"shiny2/req_genes1.csv")

# heatmap genes by condition ---------------------------------------------------

hm_per_cond <- read.csv("shiny2/genes_per_condition_hm.csv")

nConds <- length(unique(proj3@cellColData@listData[["Condition"]]))

df = list()
for (i in seq_along(1:nConds)){
  df[[i]] <- hm_per_cond[,c(1,i+1)]
  
  # select top 20 values by condition
  df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:20,1]
}

final <- do.call(rbind, df)

req_genes2 <- unlist(df)
req_genes2 <- req_genes2[!duplicated(req_genes2)]

write.csv(req_genes2,"shiny2/req_genes2.csv")

# heatmap genes by sample -- ---------------------------------------------------

hm_per_sample <- read.csv("shiny2/genes_per_sample_hm.csv")

nSamples = length(proj3@sampleColData@rownames)

df = list()
for (i in seq_along(1:nSamples)){
  df[[i]] <- hm_per_sample[,c(1,i+1)]
  
  # select top 10 values by condition
  df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:10,1]
}

final <- do.call(rbind, df)

req_genes3 <- unlist(df)
req_genes3<- req_genes3[!duplicated(req_genes3)]

write.csv(req_genes3,"shiny2/req_genes3.csv")


# pal <- paletteDiscrete(values = getCellColData(proj3)$Condition)
# 
# options(repr.plot.width=12, repr.plot.height=7)
# 
# req_DF <- as.data.frame(getCellColData(proj3))
# 
# req_table <- melt(table(req_DF$Clusters, req_DF$Condition))
# colnames(req_table) <- c("Cluster","Treatment","%cells in knn clusters")
# req_table$Cluster <- factor(
#   req_table$Cluster,
#   levels = (
#     unique(
#       req_table[order(as.numeric(gsub("C","",req_table$Cluster))),]$Cluster)
#     )
#   )
# 
# ggplot(req_table, aes(fill=Treatment, y=`%cells in knn clusters`, x= Cluster)) +
#   geom_bar(stat="identity", position = "fill") +
#   theme_classic() +
#   theme(text = element_text(size = 25)) +theme(axis.title.x=element_blank(),) +
#   scale_fill_manual(values= (pal))
