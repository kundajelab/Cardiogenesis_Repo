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
library(ArchR)
data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38
addArchRThreads(25)


library(Matrix)
library(SummarizedExperiment)
#library(tidyverse)
library(uwot)
library(edgeR)
library(FNN)
library(matrixStats)
library(Rcpp)
set.seed(1)

sparseRowVariances <- function (m){
    rM <- Matrix::rowMeans(m)
    rV <- computeSparseRowVariances(m@i + 1, m@x, rM, ncol(m))
    return(rV)
}

#Helper function for summing sparse matrix groups
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
        else {
            rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}

sparseMatTTest <- function(mat1, mat2, m0 = 0){
	#Get Population Values
	n1 <- ncol(mat1)
	n2 <- ncol(mat2)
	n <- n1 + n2
	#Sparse Row Means
	m1 <- Matrix::rowMeans(mat1, na.rm=TRUE)
	m2 <- Matrix::rowMeans(mat2, na.rm=TRUE)
	#Sparse Row Variances
	v1 <- ArchR:::computeSparseRowVariances(mat1@i + 1, mat1@x, m1, n1)
	v2 <- ArchR:::computeSparseRowVariances(mat2@i + 1, mat2@x, m2, n2)
	#Calculate T Statistic
	se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2) )
    tstat <- (m1-m2-m0)/se
	#tstat <- sqrt((n1 * n2) / n) / sqrt((n1-1)/(n-2)*v1 + (n2-1)/(n-2)*v2)
	pvalue <- 2*pt(-abs(tstat), n - 2)
	fdr <- p.adjust(pvalue, method = "fdr")
	out <- data.frame(fdr = fdr, pval = pvalue, tstat = tstat, mean1 = m1, mean2 = m2, var1 = v1, var2 = v2, n1 = n1, n2 = n2)
	return(out)
}

fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
 for(i in seq_along(fn)){
  tryCatch({
   eval(parse(text=paste0(fn[i], '<-ArchR:::', fn[i])))
  }, error = function(x){
  })
 }

# Code below adapted from ArchR function
projectLSI <- function(mat_se = NULL, LSI = NULL){  
    require(Matrix)
    set.seed(LSI$seed)

    subset_rows <- paste(rowData(mat_se)$seqnames, rowData(mat_se)$start) %in% paste(LSI$LSIFeatures$seqnames, LSI$LSIFeatures$start)
    mat <- assay(mat_se)
    mat <- mat[subset_rows,]

    #Get Same Features--whats stored here in lsi isnt exactly whats needed, so I added the lines above this to subset
    mat <- mat[LSI$idx,]

    #Binarize Matrix
    if(LSI$binarize){
        mat@x[mat@x > 0] <- 1       
    }
    
    #TF
    colSm <- Matrix::colSums(mat)
    if(any(colSm == 0)){
      exclude <- which(colSm==0)
      mat <- mat[,-exclude]
      colSm <- colSm[-exclude]
    }
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))

    #Adapted from Stuart et al.

    #IDF
    idf   <- as(LSI$nCol / LSI$rowSm, "sparseVector")

    #TF-IDF
    mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

    #Log transform TF-IDF
    mat@x <- log(mat@x * LSI$scaleTo + 1) 

    gc()

    #Clean Up Matrix
    idxNA <- Matrix::which(is.na(mat),arr.ind=TRUE)
    if(length(idxNA) > 0){
        mat[idxNA] <- 0
    }

    #Calc V
    V <- Matrix::t(mat) %*% LSI$svd$u %*% Matrix::diag(1/LSI$svd$d)

    #LSI Diagonal
    svdDiag <- matrix(0, nrow=LSI$nDimensions, ncol=LSI$nDimensions)
    diag(svdDiag) <- LSI$svd$d
    matSVD <- Matrix::t(svdDiag %*% Matrix::t(V))
    matSVD <- as.matrix(matSVD)
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))
    matSVD
}


# Load fetal project
proj_featal_invivo <- loadArchRProject(path = "path to the fetal archr project")
#Load invitro project
proj_all_invitro_peaks <- loadArchRProject(path = "path to the invitro archr project")
#load the combined project 
all_combined<-loadArchRProject('path to the combined fetal + invitro archr project')


# Load saved lsi
lsi <- getReducedDims(proj_featal_invivo, reducedDims = "IterativeLSI", returnMatrix = FALSE)

# Load Saved UMAP Manifold
umap <- getEmbedding(proj_featal_invivo, embedding = "UMAP", returnDF = FALSE)
umapManifold <- uwot::load_uwot(umap$params$uwotModel[1])


#read the matrix from the differentiation and get the matSVD returned.

sampleName='invitro sample name' #eg. H5_V2 ( Day 30 cardiomyocytes)
mat_se <- getMatrixFromArrow(
  ArrowFile = paste0("path to the invitro project arraowfiles", sampleName, ".arrow"),
  useMatrix = "TileMatrix",
  #if we need to subset and project only a part of these cells
  #cellNames = rownames(proj_all_invitro_peaks[ (proj_all_invitro_peaks$Clusters1=='M_Cardiac Fiborblast') & (proj_all_invitro_peaks$Sample=='SMC_v2'),]),
  useSeqnames = NULL,
  ArchRProj = proj_all_invitro_peaks,
  verbose = TRUE,
  binarize = TRUE
)

lsiProjection <- projectLSI(mat_se, lsi)
#UMAP Projection
set.seed(1)
umapProjection <- uwot::umap_transform(as.matrix(lsiProjection)[,1:30], umapManifold, verbose = TRUE)
#Plot Projection
refDF <- data.frame(row.names = proj_featal_invivo$CellNames , X1 = umapManifold$embedding[,1], X2 = umapManifold$embedding[,2], Type = "reference")
proDF <- data.frame(row.names = proj_all_invitro_peaks$CellNames, X1 = umapProjection[,1], X2 = umapProjection[,2], Type = "smc")
projectionDF <- rbind(refDF, proDF)
#
plotParams <- list()
plotParams$x <- projectionDF[, 'X1']
plotParams$y <- projectionDF[, 'X2']
plotParams$title <- " Colored by Clusters"
plotParams$baseSize <- 10
plotParams$xlabel <- "UMAP Dimension 1"
plotParams$ylabel <- "UMAP Dimension 2"
plotParams$color <- as.character(projectionDF$Type)
plotParams$size <- 0.2
plotParams$randomize <- TRUE
plotParams$pal <- c("reference"="#E0ECFF","smc"="#00CC00")
plotParams$labelMeans <- FALSE
gg<-do.call(ggPoint,plotParams)
gg
#
#plotPDF(gg, name = "smc_cardiomyocytes_invivo_invitro", width = 8, height = 8, ArchRProj = proj_all_invitro_peaks, addDOC = FALSE)



#Getting the peakset from the combined project to perform differential analysis
se<-getMatrixFromProject(all_combined,useMatrix = "PeakMatrix")

###################################### Enhancer and Tf differentials ######################################

#Enhancers differential 
#differentials between the nearest neighbours  
#Input Parameters
input_knn <- 25

#LSI-SVD
svdReference <- as.data.frame(lsi$matSVD) #loaded lsi
svdDisease <- as.data.frame(as.matrix(lsiProjection)) # defined from projectLSI

#KNN Nearest Neighbor using FNN #find 25 nn cells
library(FNN)
knnDisease <- get.knnx(
    data = svdReference,
    query = svdDisease,
    k = input_knn)


head(knnDisease$nn.index[,1])
uniqueIdx <- unique(as.vector(knnDisease$nn.index))
length(uniqueIdx)

idxReference <- rownames(svdReference)[uniqueIdx]
idxDisease <- colnames(mat_se)
idxDisease<-setdiff(idxDisease,setdiff(idxDisease,  rownames(all_combined)))

if(length(idxReference) > length(idxDisease)){
    idxReference <- sample(idxReference, length(idxDisease))
}else{
    idxDisease <- sample(idxDisease, length(idxReference))
}
promoterPeaks <- subjectHits(findOverlaps(resize(getTSS(all_combined), 500 * 2 + 1), getPeakSet(all_combined), ignore.strand=TRUE))
matHealthy <- assay(se[,idxReference])
matDisease <- assay(se[,idxDisease])
#Normalize to scaleTo
matNormDisease <- t(t(matDisease)/Matrix::colSums(matDisease[promoterPeaks,])) * 5000
matNormHealthy <- t(t(matHealthy)/Matrix::colSums(matHealthy[promoterPeaks,])) * 5000
#T-Test Comparisons
dfTT <- sparseMatTTest(matNormDisease, matNormHealthy)
dfTT$feature <- rownames(matNormDisease)
dfTT$log2Mean <- log2(rowMeans(cbind(dfTT$mean1, dfTT$mean2)) + 10^-4)
dfTT$log2FC <- log2((dfTT$mean1 + 10^-4)/(dfTT$mean2 + 10^-4))
plotDiff <- data.frame(row.names=row.names(dfTT),log2Mean=dfTT$log2Mean,log2FC=dfTT$log2FC,FDR=dfTT$fdr)
#plotDiff <- plotDiff[complete.cases(plotDiff),]
plotDiff$type <- "not-differential"
plotDiff$type[plotDiff$log2FC >= 1 & plotDiff$FDR <= 0.05] <- "up-regulated"
plotDiff$type[plotDiff$log2FC <= -1 & plotDiff$FDR <= 0.05] <- "do-regulated"



name_col<-'cell name'
tt<-data.frame(plotDiff$log2FC)
names(tt)<-name_col
#
tt1<-data.frame(plotDiff$FDR)
names(tt1)<-name_col
#
tt2<-data.frame(plotDiff$log2Mean)
names(tt2)<-name_col

#Obtain the differntial TFs in within the differential enhancers
enrich_se <- SummarizedExperiment(assays =  SimpleList(Log2FC = tt, FDR=tt1, log2Mean=tt2 ) )
rowData(enrich_se)<-data.frame(getPeakSet(all_combined))[c('seqnames','idx','start','end')]
metadata(enrich_se)$Params$useMatrix<-"PeakMatrix"
motifsUp <- peakAnnoEnrichment(
    seMarker = enrich_se,
    ArchRProj = all_combined,
    peakAnnotation = "Motif",
    cutOff = "Log2FC >= 0.5 & FDR <= 0.1" 
 )


motifsDo <- peakAnnoEnrichment(
    seMarker = enrich_se,
    ArchRProj = all_combined,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5" 
 )

#
df_up <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df_up <- df_up[order(df_up$mlog10Padj, decreasing = TRUE),]
df_up$rank <- seq_len(nrow(df_up))
ggUp <- ggplot(df_up, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp

df_down <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df_down <- df_down[order(df_down$mlog10Padj, decreasing = TRUE),]
df_down$rank <- seq_len(nrow(df_down))
ggDo <- ggplot(df_down, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo




