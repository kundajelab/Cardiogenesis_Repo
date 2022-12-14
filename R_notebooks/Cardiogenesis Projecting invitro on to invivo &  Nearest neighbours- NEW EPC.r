library(BiocManager)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(ggplot2)
library(TFBSTools)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(edgeR)
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


# Load normal project 
proj_featal_invivo <- loadArchRProject(path = "ARCHR_INVIVO_PROJECTPATH")
# Load the new EPC project
proj_all_invitro_peaks <- loadArchRProject(path = "ARCHR_NEWEPC_PROJECTPATH")
#load the combined object 
all_combined<-loadArchRProject('ARCHR_combined_INVIVO_NEWEPC_PROJECTPATH')


# Load saved lsi
lsi <- getReducedDims(proj_featal_invivo, reducedDims = "IterativeLSI", returnMatrix = FALSE)

# Load Saved UMAP Manifold
umap <- getEmbedding(proj_featal_invivo, embedding = "UMAP", returnDF = FALSE)
umapManifold <- uwot::load_uwot(umap$params$uwotModel[1])



#read the matrix from the differentiation and get the matSVD returned.
sampleName='EPC'
mat_se <- getMatrixFromArrow(
  ArrowFile = paste0("ARCHRPROJECTPATH/ArrowFiles/", sampleName, ".arrow"),  
  useMatrix = "TileMatrix",
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
plotParams$pal <- c("reference"="#E0ECFF","smc"="#60824F")
plotParams$labelMeans <- FALSE
gg<-do.call(ggPoint,plotParams)
gg



#Loading in the PeakMatrix to identify differential enhancers from the in vivo and new EPC combined ArchR Project
se<-getMatrixFromProject(all_combined,useMatrix = "PeakMatrix")

###################################### differentials ######################################

#Enhancers differential 
#differentials between the nearest neighbours  
#enhancers and enhacer and gene links 
#Input Parameters
input_knn <- 25

#LSI-SVD
svdReference <- as.data.frame(lsi$matSVD) #loaded lsi
svdDisease <- as.data.frame(as.matrix(lsiProjection)) # defined from projectLSI

#KNN Nearest Neighbor using FNN #find 25 nn cells
library(FNN)
set.seed(1)
knnDisease <- get.knnx(
    data = svdReference,
    query = svdDisease,
    k = input_knn)


head(knnDisease$nn.index[,1])
uniqueIdx <- unique(as.vector(knnDisease$nn.index))
length(uniqueIdx)

#Reference cells for testing
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
plotDiff$type <- "not-differential"
plotDiff$type[plotDiff$log2FC > 1 & plotDiff$FDR < 0.05] <- "up-regulated"
plotDiff$type[plotDiff$log2FC < -1 & plotDiff$FDR < 0.05] <- "do-regulated"



#Annotating Cell cluster assignments for nearest neighbour in vivo cells
cellcoldata<-getCellColData(proj_featal_invivo)
reqnames<-rownames(lsi$matSVD[as.vector(knnDisease$nn.index),])
temp_df<-data.frame(cellcoldata[reqnames,'Clusters2'])
names(temp_df)<-'clusters'
Cluterfreq<-data.frame(table(temp_df)/dim(temp_df)[1])
Cluterfreq$Rank<-rank(a$Freq)



#Setting the color scheme as per the original cellular annotations
b<-c("B_Atrial_Cardiomyocytes"="#FFAC53","C_Ventricular_Cardiomyocytes"="#F15F30","A_Myocardium"="#FFE500"
,"HA_Cardiac_Fibroblast"="#3B9AB2","F_Cardiac_fibroblast_progenitors"="#0AD7D3","K_vSMC"="#A8DDB5","L_Pericytes"="#79FFFF"
,"N_Neural_Crest"="#D0CD47","H_Early_Cardiac_fibroblast"="#90D5E4","D_Early_OFT_SMC"="#AAD962","M_Undifferentiated_epicardium"="#208A42"
,"E_valveFibroblast"="#74A9FF","G_valveFibroblast_late"="#74A9CF","J_Vasculatur_development"="#74C476","O_Endocaridum"="#C06CAB","P_Endocaridum_transitioning"="#C06CFF","T_Venal_endothelium"="#EC6BB1"
,"Q_Endocaridum_unidentifed"="#89288F","S_Capillary_endothelium"="#E6C2DC","R_Arterial_endothelium"="#F37B7D")
colordf<-data.frame(b)
colordf$temp_df<-rownames(colordf)
names(colordf)<-c('colors','temp_df')
head(colordf)


#Creating the nearest neighbour cell cluster frequency 
new_diff<-Cluterfreq
new_diff$type<-'new'
new_diff1<-merge(new_diff,colordf,by='temp_df')


#Loading the previous nearest neighbour cell cluster frequency and plotting
old_diff<-read.csv('olddifferentiation_nearestneighbournumbers.csv')
old_diff$X <- NULL
old_diff$type<-'old'
old_diff1<-merge(old_diff,colordf,by='temp_df')
newdf<-rbind(new_diff1,old_diff1)


pdf(file = "NearestNeighboursplot.pdf", width = 8, height = 6) 

p1 <- ggplot(newdf) +
    geom_col(aes(x = type, y = Freq, fill = temp_df ) ) +scale_fill_manual(values = c("B_Atrial_Cardiomyocytes"="#FFAC53","C_Ventricular_Cardiomyocytes"="#F15F30","A_Myocardium"="#FFE500"
,"HA_Cardiac_Fibroblast"="#3B9AB2","F_Cardiac_fibroblast_progenitors"="#0AD7D3","K_vSMC"="#A8DDB5","L_Pericytes"="#79FFFF"
,"N_Neural_Crest"="#D0CD47","H_Early_Cardiac_fibroblast"="#90D5E4","D_Early_OFT_SMC"="#AAD962","M_Undifferentiated_epicardium"="#208A42"
,"E_valveFibroblast"="#74A9FF","G_valveFibroblast_late"="#74A9CF","J_Vasculatur_development"="#74C476","O_Endocaridum"="#C06CAB","P_Endocaridum_transitioning"="#C06CFF","T_Venal_endothelium"="#EC6BB1"
,"Q_Endocaridum_unidentifed"="#89288F","S_Capillary_endothelium"="#E6C2DC","R_Arterial_endothelium"="#F37B7D")) + theme_classic()

p1

dev.off()






