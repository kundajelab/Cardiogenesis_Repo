library(BiocManager)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(ggplot2)
library(TFBSTools)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(edgeR)
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38
addArchRThreads(22)


pathFragments <- "path to the fragment files "
inputFiles <- list.files(pathFragments, pattern = ".gz$", full.names = TRUE)
names(inputFiles) <- gsub(".fragments.tsv.gz", "", list.files(pathFragments, pattern = ".gz$"))

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  filterTSS = 6,
  filterFrags = <setfragment number>, 
  sampleNames = names(inputFiles),
  geneAnnotation = geneAnno,
  genomeAnnotation = genomeAnno,
  force = TRUE
)

proj_1<- ArchRProject(
  ArrowFiles = ArrowFiles, 
  geneAnnotation = geneAnno,
  genomeAnnotation = genomeAnno,
  outputDirectory = "path to archr project"
)

proj_1

req_names <- read.csv('differentiation_names.csv')
req_names1<-req_names$X
proj<-subsetCells(ArchRProj=proj_1,cellNames=req_names1)

proj <- addIterativeLSI(
  ArchRProj = proj, 
  useMatrix = "TileMatrix",iterations = 4,force=TRUE
)


proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI",force=TRUE
)


proj <- addClusters(input = proj, reducedDims = "IterativeLSI", resolution = 1.5,force=TRUE)

plotList <- list()
plotList[[1]] <- plotEmbedding(ArchRProj = proj, colorBy = "colData", name = "Sample")
plotList[[2]] <- plotEmbedding(ArchRProj = proj, colorBy = "colData", name = "Clusters", plotParams = list(labelMeans=TRUE))
plotPDF(plotList = plotList, name = "UMAP-Samples-Clusters", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)
plotList



coldata<-getCellColData(proj)# Annotation of clusters 
head(coldata)


#
coldata[coldata$Clusters=='C9','Clusters1']='iPSC'
coldata[coldata$Clusters=='C10','Clusters1']='iPSC-Mes'
coldata[coldata$Clusters=='C12','Clusters1']='i-Mes'
coldata[coldata$Clusters=='C14','Clusters1']='i-CP'
coldata[coldata$Clusters=='C13','Clusters1']='i-Mes-CP'
#
coldata[coldata$Clusters=='C11','Clusters1']='i-Mes-like'
coldata[coldata$Clusters=='C3','Clusters1']='i-MyoF-like'
coldata[coldata$Clusters=='C16','Clusters1']='i-pCM'
coldata[coldata$Clusters=='C15','Clusters1']='i-CM'
coldata[coldata$Clusters=='C1','Clusters1']='i-EC'
coldata[coldata$Clusters=='C2','Clusters1']='i-EC'
#
coldata[coldata$Clusters=='C6','Clusters1']='i-EPC'
coldata[coldata$Clusters=='C4','Clusters1']='i-EPC'
coldata[coldata$Clusters=='C5','Clusters1']='i-SMC'
coldata[coldata$Clusters=='C8','Clusters1']='i-SMC'
coldata[coldata$Clusters=='C7','Clusters1']='i-CF'
coldata[(coldata$Sample=='H5_v2') & (coldata$Clusters=='C13') ,'Clusters1']='i-CM'

proj$Clusters1<-coldata$Clusters1


#paletteDiscrete
#paletteDiscrete(values = factor, set = "bear", reverse = FALSE)
plotList <- list()
#plotEmbedding(ArchRProj = proj, colorBy = "colData",name = "Sample")
plotList[[1]] <-plotEmbedding(ArchRProj = proj,pal=c("H5_v2"="#E22929","H4"="#FCB31A","EC_v2"="#FF8B74",
                                    "EPC_v2"="#60824F","SMC_v2"="#00CC00","CF_v2"="#00A08A","H3_v2"="#FBDF72",
                                    "H2_v2"="#5BB1CB","H1_v2"="#084081"),colorBy = "colData",name = "Sample")
plotList[[2]] <-plotEmbedding(ArchRProj = proj,pal=c("H5_v2"="#E22929","H4"="#FCB31A","EC_v2"="#FF8B74",
                                    "EPC_v2"="#60824F","SMC_v2"="#00CC00","CF_v2"="#00A08A","H3_v2"="#FBDF72",
                                    "H2_v2"="#5BB1CB","H1_v2"="#084081"),colorBy = "colData",name = "Sample", labelMeans=FALSE)

#plotEmbedding(ArchRProj = proj, pal=ArchRPalettes$bear,colorBy = "colData",name = "Sample")
plotPDF(plotList = plotList, name = "UMAP-Samples-Clusters", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)
plotList

#plotting gene scores on umap
proj <- addImputeWeights(ArchRProj = proj)

markerGenes  <- c("POU5F1","NANOG","MESP1","MESP2","ISL1","TBX5","TNNT2","TTN","MYL2","HAND1","HAND2","PECAM1","CD34","CDH5","TCF21",
    "NOTCH1","NOTCH4","DCN","COL1A1","LUM","MYH11","SOX17","FOXA1","TAGLN", "PDGFRB", "MYH6", "MYH7", "MYL7"
  )


#Plot the UMAP Embedding with Marker Genes Overlayed w/ Imputation
plotList <- list()
plotList[[1]] <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = markerGenes, imputeWeights = getImputeWeights(proj))
plotPDF(plotList = plotList, name = "UMAP-Marker-Gene-Scores-w-Imputation", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)

#Plot the UMAP Embedding with Marker Genes Overlayed w/o Imputation
plotList <- list()
plotList[[1]] <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = markerGenes, imputeWeights = getImputeWeights(proj))
plotPDF(plotList = plotList, name = "UMAP-Marker-Gene-integration-w-Imputation", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)
plotList

#getting reproducible peak 
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters1",force=TRUE)
#Call Reproducible Peaks w/ Macs2 (~5-10 minutes)
proj <- addReproduciblePeakSet(ArchRProj = proj,groupBy = "Clusters1",force=TRUE)
#Add Peak Matrix
proj <- addPeakMatrix(ArchRProj = proj)



#Identify Marker Peaks
markersPeaks <- markerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "Clusters1")

#Visualize Markers as a heatmap
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 1",
  pal = paletteContinuous(set = "blueYellow")
)
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 12, ArchRProj = proj, addDOC = FALSE)

#addMotifAnnotations
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force = TRUE)


#Identify Motif Enrichments
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
heatmapEM <- enrichHeatmap(enrichMotifs,pal = paletteContinuous(set = "solarExtra", n = 100))
plotPDF(heatmapEM, name = "Motifs-Enrich-Heatmap", width = 8, height = 10, ArchRProj = proj, addDOC = FALSE)


#Adding deviations
proj <- addBgdPeaks(proj,force = TRUE)
#Add chromVAR Deviations (~20-25 min if using CisBP Motif Set)
proj <- addDeviationsMatrix(ArchRProj = proj, peakAnnotation = "Motif",threads=15,force = TRUE)










