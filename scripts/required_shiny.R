library("ArchR")
library("GenomicRanges")
library('BSgenome')
library("org.Rn.eg.db")
library("dplyr")                                    # Load dplyr package
library("plyr")                                     # Load plyr package
library("readr")  
library("qdap")
library('Seurat')
library("ComplexHeatmap")
library("circlize")
library(data.table)

setwd('~/fastq/archr/')

proj3 <- loadArchRProject(path = "./another_saving_test_ArchRProject/")

proj3 <- addIterativeLSI(
  ArchRProj = proj3,
  useMatrix = 'TileMatrix',
  name = 'IterativeLSI',
  iterations = 2, 
  clusterParams = list(
    resolution = c(1), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

proj3 <- addClusters(
  input = proj3,
  method = 'Seurat',
  name = 'Clusters',
  resolution = c(1), 
  force = TRUE
)

proj3 <- addUMAP(
  ArchRProj = proj3, 
  name = 'UMAP', 
  nNeighbors = 30, 
  minDist = 0.0,  
  metric = 'cosine',
  force = TRUE
)

# clusters ---------------------------------------------------------------------
markersGS_clusters <- getMarkerFeatures(
  ArchRProj = proj3, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest",#"wilcoxon"
  
)

heatmapGS_cluster <- plotMarkerHeatmap(
  seMarker = markersGS_clusters, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.20",
  plotLog2FC = TRUE,
  transpose = F,
  returnMatrix = TRUE
)

saveRDS(markers_clusters,"markersGS_clusters.rds")
write.csv(heatmapGS,"genes_per_cluster_hm.csv")

# samples  ---------------------------------------------------------------------

markersGS_samples <- getMarkerFeatures(
  ArchRProj = proj3, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest",#"wilcoxon"
)

proj3 <- addImputeWeights(proj3, reducedDims= "Harmony")                     

options(repr.plot.width=7, repr.plot.height=7)

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.20", plotLog2FC = TRUE,
  #   labelMarkers = seleceted_markers,
  transpose = F,  returnMatrix = FALSE
)

heatmapGS

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.20", plotLog2FC = TRUE,
  #   labelMarkers = seleceted_markers,
  transpose = F,  returnMatrix = TRUE
)
write.csv(heatmapGS,"genes_per_sample_hm.csv")      


saveRDS(markersGS,"markersGS_sample.rds")



