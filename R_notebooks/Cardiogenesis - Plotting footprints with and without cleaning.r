library(BiocManager)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(ggplot2)
library(TFBSTools)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
#library(SeuratData)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(edgeR)

data("geneAnnoHg38")
data("genomeAnnoHg38")
geneAnno <- geneAnnoHg38
genomeAnno <- genomeAnnoHg38
addArchRThreads(12)


#this load archr dev mode
fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
  for(i in seq_along(fn)){
    tryCatch({
      eval(parse(text=paste0(fn[i], '<-ArchR:::', fn[i])))
    }, error = function(x){
    })
  }

##Loading ArchR Projects as described in previous notebooks 

mef_cleaned=read.table('getting the cleaned motifs as a bed file for Mef',sep='\t',header=TRUE)
names(mef_cleaned)<-c('seqnames','start','end','tfname')
mef_cleaned$strand='*'
#just a dummy score
mef_cleaned$score=1.2
head(mef_cleaned)
mef_cleaned<-GRanges(mef_cleaned)


#modfiying the getFootprints
getFootprints<-function (ArchRProj = NULL, positions = NULL, plotName = "Plot-Footprints", 
    groupBy = "Clusters", useGroups = NULL, flank = 250, minCells = 25, 
    nTop = NULL, threads = getArchRThreads(), verbose = TRUE, 
    logFile = createLogFile("getFootprints")) 
{
    .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
    .validInput(input = positions, name = "positions", valid = c("grangeslist"))
    .validInput(input = plotName, name = "plotName", valid = c("character"))
    .validInput(input = groupBy, name = "groupBy", valid = c("character"))
    .validInput(input = useGroups, name = "useGroups", valid = c("character", 
        "null"))
    .validInput(input = flank, name = "flank", valid = c("integer"))
    .validInput(input = minCells, name = "minCells", valid = c("integer"))
    .validInput(input = nTop, name = "nTop", valid = c("integer", 
        "null"))
    .validInput(input = threads, name = "threads", valid = c("integer"))
    .validInput(input = verbose, name = "verbose", valid = c("boolean"))
    .validInput(input = logFile, name = "logFile", valid = c("character"))
    tstart <- Sys.time()
    .startLogging(logFile = logFile)
    .logThis(mget(names(formals()), sys.frame(sys.nframe())), 
        "Input-Parameters", logFile = logFile)
    coverageMetadata <- .getCoverageMetadata(ArchRProj = ArchRProj, 
        groupBy = groupBy, minCells = minCells)
    coverageParams <- .getCoverageParams(ArchRProj = ArchRProj, 
        groupBy = groupBy)
    kmerLength <- coverageParams$kmerLength
    .logThis(coverageMetadata, "coverageMetadata", logFile = logFile)
    .logThis(coverageParams, "coverageParams", logFile = logFile)
    if (!is.null(useGroups)) {
        if (sum(coverageMetadata[, 1] %in% useGroups) == 0) {
            stop("No Groups found matching useGroups!")
        }
        coverageMetadata <- coverageMetadata[coverageMetadata[, 
            1] %in% useGroups, ]
    }
    genome <- getGenome(ArchRProj)
    .requirePackage(genome)
    .requirePackage("Biostrings", source = "bioc")
    BSgenome <- eval(parse(text = genome))
    BSgenome <- validBSgenome(BSgenome)
    .logDiffTime("Computing Kmer Bias Table", tstart, verbose = verbose, 
        logFile = logFile)
    kmerTableList <- .kmerPositionFrequency(featureList = positions, 
        genome = BSgenome, flank = flank, k = kmerLength, threads = 1, 
        verbose = FALSE, logFile = logFile)
    .logDiffTime("Computing Footprints", tstart, verbose = verbose, 
        logFile = logFile)
    footprintList <- .computeFootprints(featureList = positions, 
        coverageFiles = coverageMetadata$File, flank = flank, 
        threads = threads, verbose = FALSE, logFile = logFile)
    .logDiffTime("Computing Footprints Bias", tstart, verbose = verbose, 
        logFile = logFile)
    footprintBiasList <- .computeFootprintsBias(kmerTableList = kmerTableList, 
        coverageFiles = coverageMetadata$File, threads = threads, 
        verbose = FALSE)
    .logDiffTime("Summarizing Footprints", tstart, verbose = verbose, 
        logFile = logFile)
    footAssay <- lapply(seq_along(positions), function(x) {
        footMat <- lapply(seq_along(footprintList), function(y) {
            footprintList[[y]][, x]
        }) %>% Reduce("cbind", .)
        colnames(footMat) <- coverageMetadata$Name
        biasMat <- lapply(seq_along(footprintBiasList), function(y) {
            footprintBiasList[[y]][, x]
        }) %>% Reduce("cbind", .)
        colnames(biasMat) <- coverageMetadata$Name
        rbind(footMat, biasMat)
    }) %>% SimpleList
    names(footAssay) <- names(positions)
    rm(footprintList, footprintBiasList)
    gc()
    rowData <- DataFrame(x = c(seq(-flank, flank), seq(-flank, 
        flank)), type = c(rep("footprint", flank * 2 + 1), rep("bias", 
        flank * 2 + 1)))
    se <- SummarizedExperiment::SummarizedExperiment(assays = footAssay, 
        colData = coverageMetadata, rowData = rowData)
    metadata(se)$Params <- SimpleList(kmerLength = kmerLength, 
        flank = flank, date = Sys.Date())
    return(se)
}

motifPositions <- getPositions(proj)
motifs <- c("MEF2C")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs

mef<-grep('MEF2C', names(motifPositions), value = TRUE)
mef
                              
mef_cleaned=read.table('/oak/stanford/groups/akundaje/laks/illuminafiles/cardiomyogenesis/archr_input_output/aggregated_motif_files/VCM_GATA4_stringentcleaned.csv',
                      sep='\t',header=TRUE)
names(mef_cleaned)<-c('seqnames','start','end','tfname')
mef_cleaned$strand='*'
#dummy score
mef_cleaned$score=1.2
head(mef_cleaned)
mef_cleaned<-GRanges(mef_cleaned)

                              
                              


seFoot <- getFootprints(
  ArchRProj = proj, 
  #if we need to the precleaned motifs, just uncomment the following line and comment the line after that 
  #positions = motifPositions[markerMotifs],
  positions=GRangesList(mef_cleaned),
  groupBy='Clusters1'
    #useGroups=c('Ventricular_Cardiomyocytes')    
)

#renaming the se object. It wont work without giving it a name
names(assays(seFoot))<-'MEF2C'

.groupMeans<-function (mat = NULL, groups = NULL, na.rm = TRUE, sparse = FALSE) 
{
    #print('inside group means')
    #print(groups)
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowMeans(mat[, which(groups == x), drop = F], 
                na.rm = na.rm)
        }
        else {
            rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    #print('done with the group means')
    #gm<-matrix(gm)
    #print(unique(groups))
    colnames(gm) <- unique(groups)
    return(gm)
}

.groupSds<-function (mat = NULL, groups = NULL, na.rm = TRUE, sparse = FALSE) 
{
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gs <- lapply(unique(groups), function(x) {
        if (sparse) {
            matrixStats::rowSds(as.matrix(mat[, which(groups == 
                x), drop = F]), na.rm = na.rm)
        }
        else {
            matrixStats::rowSds(mat[, which(groups == x), drop = F], 
                na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    #gs<-matrix(gs)
    colnames(gs) <- unique(groups)
    return(gs)
}

.ggFootprint<-function (seFoot = NULL, name = NULL, pal = NULL, smoothWindow = NULL, 
    flank = NULL, flankNorm = NULL, baseSize = 6, normMethod = NULL, 
    logFile = NULL) 
{
    errorList <- list()
    rowDF <- SummarizedExperiment::rowData(seFoot)
    footMat <- .getAssay(seFoot[BiocGenerics::which(rowDF[, 2] == 
        "footprint"), ], name)
    biasMat <- .getAssay(seFoot[BiocGenerics::which(rowDF[, 2] == 
        "bias"), ], name)
    footDF <- rowDF[BiocGenerics::which(rowDF[, 2] == "footprint"), 
        ]
    biasDF <- rowDF[BiocGenerics::which(rowDF[, 2] == "bias"), 
        ]
    errorList$footMat <- footMat
    errorList$biasMat <- biasMat
    errorList$footDF <- footDF
    errorList$biasDF <- biasDF
    if (!is.null(smoothWindow)) {
        .logMessage("Applying smoothing window to footprint", 
            logFile = logFile)
        footMat <- apply(footMat, 2, function(x) .centerRollMean(x, 
            smoothWindow))
        biasMat <- apply(biasMat, 2, function(x) .centerRollMean(x, 
            smoothWindow))
    }
    .logMessage("Normalizing by flanking regions", logFile = logFile)
    idx <- which(abs(footDF$x) >= flank - flankNorm)
    footMat <- t(t(footMat)/colMeans(footMat[idx, , drop = FALSE]))
    biasMat <- t(t(biasMat)/colMeans(biasMat[idx, , drop = FALSE]))
    errorList$footMatNorm <- footMat
    errorList$biasMatNorm <- footMat
    if (tolower(normMethod) == "none") {
        title <- ""
    }
    else if (tolower(normMethod) == "subtract") {
        title <- "Tn5 Bias Subtracted\n"
        footMat <- footMat - biasMat
        #print('Here in ggfootprint Subtract')
    }
    else if (tolower(normMethod) == "divide") {
        title <- "Tn5 Bias Divided\n"
        footMat <- footMat/biasMat
    }
    else {
        stop("normMethod not recognized!")
    }
    .logMessage(paste0("NormMethod = ", normMethod), logFile = logFile)
    #print('Here in ggfootprint block done ')
    print(SummarizedExperiment::colData(seFoot)$Group)
    print(dim(footMat))
    footMatMean <- .groupMeans(footMat, SummarizedExperiment::colData(seFoot)$Group)
    #footMatMean <- .groupMeans(footMat, 'Ventricular_Cardiomyocytes')
    footMatSd <- .groupSds(footMat, SummarizedExperiment::colData(seFoot)$Group)
    biasMatMean <- .groupMeans(biasMat, SummarizedExperiment::colData(seFoot)$Group)
    biasMatSd <- .groupSds(biasMat, SummarizedExperiment::colData(seFoot)$Group)
    smoothFoot <- rowMaxs(apply(footMatMean, 2, function(x) .centerRollMean(x, 
        11)))
    errorList$footMatMean <- footMatMean
    errorList$footMatSd <- footMatSd
    errorList$biasMatMean <- biasMatMean
    errorList$biasMatSd <- biasMatSd
    errorList$smoothFoot <- smoothFoot
    plotIdx <- seq_len(nrow(footMatMean))
    plotFootDF <- lapply(seq_len(ncol(footMatMean)), function(x) {
        data.frame(x = footDF$x, mean = footMatMean[, x], sd = footMatSd[, 
            x], group = colnames(footMatMean)[x])[plotIdx, , 
            drop = FALSE]
    }) %>% Reduce("rbind", .)
    plotFootDF$group <- factor(paste0(plotFootDF$group), levels = unique(gtools::mixedsort(paste0(plotFootDF$group))))
    plotBiasDF <- lapply(seq_len(ncol(biasMatMean)), function(x) {
        data.frame(x = biasDF$x, mean = biasMatMean[, x], sd = biasMatSd[, 
            x], group = colnames(biasMatMean)[x])[plotIdx, , 
            drop = FALSE]
    }) %>% Reduce("rbind", .)
    plotBiasDF$group <- factor(paste0(plotBiasDF$group), levels = unique(gtools::mixedsort(paste0(plotBiasDF$group))))
    errorList$plotFootDF <- plotFootDF
    errorList$plotBiasDF <- plotBiasDF
    out <- tryCatch({
        if (is.null(pal)) {
            pal <- paletteDiscrete(values = gtools::mixedsort(SummarizedExperiment::colData(seFoot)$Group))
        }
        plotMax <- plotFootDF[order(plotFootDF$mean, decreasing = TRUE), 
            ]
        plotMax <- plotMax[abs(plotMax$x) > 20 & abs(plotMax$x) < 
            50, ]
        plotMax <- plotMax[!duplicated(plotMax$group), ]
        plotMax <- plotMax[seq_len(ceiling(nrow(plotMax)/4)), 
            ]
        plotMax$x <- 25
        ggFoot <- ggplot(plotFootDF, aes(x = x, y = mean, color = group)) + 
            geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, 
                linetype = NA, fill = group), alpha = 0.4) + 
            geom_line() + scale_color_manual(values = pal) + 
            scale_fill_manual(values = pal) + xlab("Distance to motif center (bp)") + 
            coord_cartesian(expand = FALSE, ylim = c(quantile(plotFootDF$mean, 
                1e-04), 1.15 * quantile(smoothFoot, 0.999)), 
                xlim = c(min(plotFootDF$x), max(plotFootDF$x))) + 
            theme_ArchR(baseSize = baseSize) + ggtitle(name) + 
            guides(fill = FALSE) + guides(color = FALSE) + ylab(paste0(title, 
            "Normalized Insertions")) + ggrepel::geom_label_repel(data = plotMax, 
            aes(label = group), size = 3, xlim = c(75, NA))
        ggBias <- ggplot(plotBiasDF, aes(x = x, y = mean, color = group)) + 
            geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, 
                linetype = NA, fill = group), alpha = 0.4) + 
            geom_line() + scale_color_manual(values = pal) + 
            scale_fill_manual(values = pal) + xlab("Distance to motif center (bp)") + 
            coord_cartesian(expand = FALSE, ylim = c(quantile(plotBiasDF$mean, 
                1e-04), 1.05 * quantile(plotBiasDF$mean, 0.999)), 
                xlim = c(min(plotBiasDF$x), max(plotBiasDF$x))) + 
            theme_ArchR(baseSize = baseSize) + ylab("Tn5-Bias Normalized Insertions") + 
            theme(legend.position = "bottom", legend.box.background = element_rect(color = NA))
        ggAlignPlots(ggFoot, .ggSmallLegend(ggBias), sizes = c(2, 
            1), draw = FALSE)
    }, error = function(e) {
        .logError(e, fn = ".ggFootprint", info = name, errorList = errorList, 
            logFile = logFile)
    })
    out
}


plotFootprints<-function (seFoot = NULL, names = NULL, pal = NULL, flank = 250, 
    flankNorm = 50, normMethod = "Subtract", smoothWindow = NULL, 
    baseSize = 6, plot = TRUE, ArchRProj = NULL, plotName = paste0("Plot-Footprints-", 
        normMethod), height = 6, width = 4, addDOC = TRUE, force = FALSE, 
    logFile = createLogFile("plotFootprints")) 
{
    .validInput(input = seFoot, name = "seFoot", valid = c("SummarizedExperiment"))
    .validInput(input = names, name = "names", valid = c("character", 
        "null"))
    .validInput(input = pal, name = "pal", valid = c("palette", 
        "null"))
    .validInput(input = flank, name = "flank", valid = c("integer"))
    .validInput(input = flankNorm, name = "flankNorm", valid = c("integer"))
    .validInput(input = normMethod, name = "normMethod", valid = c("character"))
    .validInput(input = smoothWindow, name = "smoothWindow", 
        valid = c("integer", "null"))
    .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
    .validInput(input = plot, name = "plot", valid = c("boolean"))
    .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
    .validInput(input = plotName, name = "plotName", valid = c("character"))
    .validInput(input = height, name = "height", valid = c("numeric"))
    .validInput(input = width, name = "width", valid = c("numeric"))
    .validInput(input = addDOC, name = "addDOC", valid = c("boolean"))
    .validInput(input = force, name = "force", valid = c("boolean"))
    .validInput(input = logFile, name = "logFile", valid = c("character"))
    tstart <- Sys.time()
    .startLogging(logFile = logFile)
    if (is.null(names)) {
        names <- names(assays(seFoot))
    }
    if (length(names) > 25) {
        if (!plot) {
            if (force) {
                .logMessage("Plotting more than 25 footprints can create large storage of ggplots. Continuing since force = TRUE", 
                  verbose = TRUE, logFile = logFile)
            }
            else {
                .logStop("Plotting more than 25 footprints can create large storage of ggplots. Stopping since force = FALSE", 
                  logFile = logFile)
            }
        }
    }
    if (plot) {
        name <- gsub("\\.pdf", "", plotName)
        if (is.null(ArchRProj)) {
            outDir <- "Plots"
        }
        else {
            ArchRProj <- .validArchRProject(ArchRProj)
            outDir <- file.path(getOutputDirectory(ArchRProj), 
                "Plots")
        }
        dir.create(outDir, showWarnings = FALSE)
        if (addDOC) {
            doc <- gsub(":", "-", stringr::str_split(Sys.time(), 
                pattern = " ", simplify = TRUE)[1, 2])
            filename <- file.path(outDir, paste0(name, "_Date-", 
                Sys.Date(), "_Time-", doc, ".pdf"))
        }
        else {
            filename <- file.path(outDir, paste0(name, ".pdf"))
        }
        pdf(filename, width = width, height = height, useDingbats = FALSE)
    }
    ggList <- lapply(seq_along(names), function(x) {
        print(x)
        print(names[1])
        .logDiffTime(sprintf("Plotting Footprint : %s (%s of %s)", 
            names[x], x, length(names)), t1 = tstart, logFile = logFile, 
            verbose = TRUE)
        gg <- .ggFootprint(seFoot = seFoot, name = names[x], 
            pal = pal, smoothWindow = smoothWindow, flank = flank, 
            flankNorm = flankNorm, baseSize = baseSize, normMethod = normMethod, 
            logFile = logFile)
        print(gg)
        if (plot) {
            if (x != 1) {
                grid::grid.newpage()
            }
            grid::grid.draw(gg)
            return(0)
        }
        else {
            return(gg)
        }
    })
    .endLogging(logFile = logFile)
    if (!plot) {
        names(ggList) <- names
        ggList
    }
    else {
        dev.off()
        return(invisible(0))
    }
}

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Divide",
  plotName = "MEF2C-Footprints-stringentcleaned",
  addDOC = FALSE
  ,
    smoothWindow = 15,
    plot=TRUE,
    pal=c(
"Atrial_Cardiomyocytes"="#FFAC53","Ventricular_Cardiomyocytes"="#F15F30","Myocardium"="#FFE500"
,"Cardiac_Fibroblast"="#3B9AB2","Cardiac_fibroblast_progenitors"="#0AD7D3","vSMC"="#A8DDB5","Pericytes"="#79FFFF"
,"Neural_Crest"="#D0CD47","Early_Cardiac_fibroblast"="#90D5E4","Early_OFT_SMC"="#AAD962","Undifferentiated_epicardium"="#208A42"
,"Endocardial_cushion"="#74A9FF","Endocardial_cushion_late"="#74A9CF","Vasculatur_development"="#74C476","Endocaridum"="#C06CAB","Endocaridum_transitioning"="#C06CFF","Venal_endothelium"="#EC6BB1"
,"Endocaridum_unidentifed"="#89288F","Capillary_endothelium"="#E6C2DC","Arterial_endothelium"="#F37B7D"
)
)



