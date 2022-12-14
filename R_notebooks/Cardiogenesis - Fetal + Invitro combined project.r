library(BiocManager)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(ggplot2)
library(TFBSTools)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(SeuratData)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(edgeR)
library(metaMA)
library(cicero)
library(ggbiplot)
library(sctransform)
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38
addArchRThreads(28)


pathFragments <- "path to all fragment files from both invitro differentiation and fetal heart"
inputFiles <- list.files(pathFragments, pattern = ".gz$", full.names = TRUE)
names(inputFiles) <- gsub(".fragments.tsv.gz", "", list.files(pathFragments, pattern = ".gz$"))

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  filterTSS = 6,
  filterFrags = 1000, 
  sampleNames = names(inputFiles),
  geneAnnotation = geneAnno,
  genomeAnnotation = genomeAnno,
  force = FALSE
)

proj_1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  geneAnnotation = geneAnno,
  genomeAnnotation = genomeAnno,
  outputDirectory = "path to the combined archr project"
)

proj_1

invivo_cells<-read.csv('path to precleaned fetal barcodes')
invitro_cells<-read.csv('path to precleaned invitro differentiation barcodes')

head(invivo_cells)
head(invitro_cells)


req_names<-c(as.character(invitro_cells$X),as.character(invivo_cells$X))
length(req_names)
head(req_names)

proj<-subsetArchRProject(ArchRProj=proj_1,cells=req_names,
                        outputDirectory='path to subsetted project')

proj

proj <- addIterativeLSI(
  ArchRProj = proj, 
  useMatrix = "TileMatrix",force=TRUE,iterations = 4
)


proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI",force=TRUE
)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI",  resolution = 0.6,force=TRUE)

plotList <- list()
plotList[[1]] <- plotEmbedding(ArchRProj = proj, colorBy = "colData", name = "Sample")
plotList[[2]] <- plotEmbedding(ArchRProj = proj, colorBy = "colData", name = "Clusters", plotParams = list(labelMeans=TRUE))


#getting reproducible peak 
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters")

#Call Reproducible Peaks w/ Macs2 (~5-10 minutes)
proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = "Clusters")

#Add Peak Matrix
proj <- addPeakMatrix(ArchRProj = proj)



proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")

###chromvar

proj <- addBgdPeaks(proj,force = TRUE)
#Add chromVAR Deviations (~20-25 min if using CisBP Motif Set)
proj <- addDeviationsMatrix(ArchRProj = proj, peakAnnotation = "Motif")



