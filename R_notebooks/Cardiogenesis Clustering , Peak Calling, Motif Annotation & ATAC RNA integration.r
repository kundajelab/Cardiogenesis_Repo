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


fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
  for(i in seq_along(fn)){
    tryCatch({
      eval(parse(text=paste0(fn[i], '<-ArchR:::', fn[i])))
    }, error = function(x){
    })
  }

pathFragments <- "path to fragment files"
inputFiles <- list.files(pathFragments, pattern = ".gz$", full.names = TRUE)
names(inputFiles) <- gsub(".fragments.tsv.gz", "", list.files(pathFragments, pattern = ".gz$"))

#
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  filterTSS = 6,
  filterFrags = 1000, 
  sampleNames = names(inputFiles),
  geneAnnotation = geneAnno,
  genomeAnnotation = genomeAnno
  force = TRUE
)

proj_1<- ArchRProject(
  ArrowFiles = ArrowFiles, 
  geneAnnotation = geneAnno,
  genomeAnnotation = genomeAnno,
  outputDirectory = "project path"
)

proj_1

###PreIdentified Macrophage barcodes are removed from the main analysis
req_barcodes <- read.csv('Final_barcodes_fetal_heart_NOMICROPHAGES.csv')
rownames(req_barcodes)<-req_barcodes$X
proj<-subsetCells(ArchRProj=proj_1,cellNames=rownames(req_barcodes))


proj <- addIterativeLSI(
  ArchRProj = proj, 
  useMatrix = "TileMatrix",force=TRUE,iterations = 4
)


proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI",force=TRUE
)



proj <- addClusters(input = proj, reducedDims = "IterativeLSI", resolution =2,force=TRUE)

plotList <- list()
plotList[[1]] <- plotEmbedding(ArchRProj = proj, colorBy = "colData", name = "Sample", labelMeans=FALSE)
plotList[[2]] <- plotEmbedding(ArchRProj = proj, colorBy = "colData", name = "Clusters2",keepAxis=FALSE, labelMeans=TRUE)
#plotPDF(plotList = plotList, name = "UMAP-Samples-Clusters", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)
plotList

coldata<-getCellColData

#Renaming clusters to more meaningful names and collapsing clusters of same cell types
coldata[coldata$Sample=='F19_v2','Sample1']='PCW19'
coldata[coldata$Sample=='F8_v2','Sample1']='PCW8'
coldata[coldata$Sample=='F6_v2','Sample1']='PCW6'
#
head(coldata)

#
coldata[coldata$Clusters=='C5','Clusters1']='Myocardium'
coldata[coldata$Clusters=='C1','Clusters1']='Ventricular_Cardiomyocytes'
coldata[coldata$Clusters=='C2','Clusters1']='Ventricular_Cardiomyocytes'
coldata[coldata$Clusters=='C4','Clusters1']='Ventricular_Cardiomyocytes'
coldata[coldata$Clusters=='C3','Clusters1']='Atrial_Cardiomyocytes'
#
coldata[coldata$Clusters=='C11','Clusters1']='Endocaridum'
coldata[coldata$Clusters=='C6','Clusters1']='Endocaridum_transitioning'
coldata[coldata$Clusters=='C7','Clusters1']='Endocaridum_unidentifed'
coldata[coldata$Clusters=='C8','Clusters1']='Arterial_endothelium'
coldata[coldata$Clusters=='C9','Clusters1']='Capillary_endothelium'
coldata[coldata$Clusters=='C10','Clusters1']='Venal_endothelium'
#
coldata[coldata$Clusters=='C21','Clusters1']='Neural_Crest'
coldata[coldata$Clusters=='C12','Clusters1']='Cardiac_fibroblast_progenitors'
coldata[coldata$Clusters=='C13','Clusters1']='Endocardial_cushion'
coldata[coldata$Clusters=='C14','Clusters1']='Early_OFT_SMC'
coldata[coldata$Clusters=='C22','Clusters1']='Endocardial_cushion_late'
coldata[coldata$Clusters=='C20','Clusters1']='Early_Cardiac_fibroblast'
coldata[coldata$Clusters=='C23','Clusters1']='Vasculatur_development'
coldata[coldata$Clusters=='C17','Clusters1']='Pericytes' 
coldata[coldata$Clusters=='C18','Clusters1']='vSMC'
coldata[coldata$Clusters=='C15','Clusters1']='Cardiac_Fibroblast'
coldata[coldata$Clusters=='C16','Clusters1']='Cardiac_Fibroblast'
coldata[coldata$Clusters=='C19','Clusters1']='Undifferentiated_epicardium'
#

#
coldata[coldata$Clusters=='C5','Clusters2']='A_Myocardium'
coldata[coldata$Clusters=='C1','Clusters2']='C_Ventricular_Cardiomyocytes'
coldata[coldata$Clusters=='C2','Clusters2']='C_Ventricular_Cardiomyocytes'
coldata[coldata$Clusters=='C4','Clusters2']='C_Ventricular_Cardiomyocytes'
coldata[coldata$Clusters=='C3','Clusters2']='B_Atrial_Cardiomyocytes'
#
coldata[coldata$Clusters=='C11','Clusters2']='O_Endocaridum'
coldata[coldata$Clusters=='C6','Clusters2']='P_Endocaridum_transitioning'
coldata[coldata$Clusters=='C7','Clusters2']='Q_Endocaridum_unidentifed'
coldata[coldata$Clusters=='C8','Clusters2']='R_Arterial_endothelium'
coldata[coldata$Clusters=='C9','Clusters2']='S_Capillary_endothelium'
coldata[coldata$Clusters=='C10','Clusters2']='T_Venal_endothelium'
#
coldata[coldata$Clusters=='C21','Clusters2']='N_Neural_Crest'
coldata[coldata$Clusters=='C12','Clusters2']='F_Cardiac_fibroblast_progenitors'
coldata[coldata$Clusters=='C13','Clusters2']='E_valveFibroblast'
coldata[coldata$Clusters=='C14','Clusters2']='D_Early_OFT_SMC'
coldata[coldata$Clusters=='C22','Clusters2']='G_valveFibroblast_late'
coldata[coldata$Clusters=='C20','Clusters2']='H_Early_Cardiac_fibroblast'
coldata[coldata$Clusters=='C23','Clusters2']='J_Vasculatur_development'
coldata[coldata$Clusters=='C17','Clusters2']='L_Pericytes' 
coldata[coldata$Clusters=='C18','Clusters2']='K_vSMC'
coldata[coldata$Clusters=='C15','Clusters2']='HA_Cardiac_Fibroblast'
coldata[coldata$Clusters=='C16','Clusters2']='HA_Cardiac_Fibroblast'
coldata[coldata$Clusters=='C19','Clusters2']='M_Undifferentiated_epicardium'

proj$Sample1<-coldata$Sample1
proj$Clusters1<-coldata$Clusters1
proj$Clusters2<-coldata$Clusters2

# # #"F6_v2"="#F37B7D","F8_v2"="#FCBF6E","F19_v2"="#90D5E4"), name = "Sample", labelMeans=FALSE
## # #
plotList <- list()
plotList[[1]]<-plotEmbedding(ArchRProj = proj, colorBy = "colData",pal=c("F6_v2"="#4EB3D3","F8_v2"="#F59899","F19_v2"="#6CD3A7"), name = "Sample", labelMeans=FALSE)
plotList[[2]]<-plotEmbedding(ArchRProj = proj, colorBy = "colData", name = "Clusters2",pal=c(
"B_Atrial_Cardiomyocytes"="#FFAC53","C_Ventricular_Cardiomyocytes"="#F15F30","A_Myocardium"="#FFE500"
,"HA_Cardiac_Fibroblast"="#3B9AB2","F_Cardiac_fibroblast_progenitors"="#0AD7D3","K_vSMC"="#A8DDB5","L_Pericytes"="#79FFFF"
,"N_Neural_Crest"="#D0CD47","H_Early_Cardiac_fibroblast"="#90D5E4","D_Early_OFT_SMC"="#AAD962","M_Undifferentiated_epicardium"="#208A42"
,"E_valveFibroblast"="#74A9FF","G_valveFibroblast_late"="#74A9CF","J_Vasculatur_development"="#74C476","O_Endocaridum"="#C06CAB","P_Endocaridum_transitioning"="#C06CFF","T_Venal_endothelium"="#EC6BB1"
,"Q_Endocaridum_unidentifed"="#89288F","S_Capillary_endothelium"="#E6C2DC","R_Arterial_endothelium"="#F37B7D"
) ,keepAxis=FALSE, labelMeans=TRUE)
plotList[[3]]<-plotEmbedding(ArchRProj = proj, colorBy = "colData", name = "Clusters2",pal=c(
"B_Atrial_Cardiomyocytes"="#FFAC53","C_Ventricular_Cardiomyocytes"="#F15F30","A_Myocardium"="#FFE500"
,"HA_Cardiac_Fibroblast"="#3B9AB2","F_Cardiac_fibroblast_progenitors"="#0AD7D3","K_vSMC"="#A8DDB5","L_Pericytes"="#79FFFF"
,"N_Neural_Crest"="#D0CD47","H_Early_Cardiac_fibroblast"="#90D5E4","D_Early_OFT_SMC"="#AAD962","M_Undifferentiated_epicardium"="#208A42"
,"E_valveFibroblast"="#74A9FF","G_valveFibroblast_late"="#74A9CF","J_Vasculatur_development"="#74C476","O_Endocaridum"="#C06CAB","P_Endocaridum_transitioning"="#C06CFF","T_Venal_endothelium"="#EC6BB1"
,"Q_Endocaridum_unidentifed"="#89288F","S_Capillary_endothelium"="#E6C2DC","R_Arterial_endothelium"="#F37B7D"
) ,keepAxis=FALSE, labelMeans=FALSE)
plotPDF(plotList = plotList, name = "UMAP-Samples-Clusters_LABELLING", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)
plotList

#Getting Marker Genes based on gene score for clusters
markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters1",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

#Visualizing gene scores in UMAP
proj <- addImputeWeights(ArchRProj = proj)

markerGenes  <- c("TNNT2","ACTN2","MYH6","MYH7","DCN","MYH11","PDGFRA","PDGFRB","TCF21","TAGLN","COL1A1","PRDM16","MYH11","TNNT2","TTN","GATA4","MEF2C","TBX5","WT1"
                  ,"TBX18","PECAM1","CDH5","CDH11","GJA5","NPPA","MYL6","MYL7","MYL2","NPPB","SELE","DCN","LUM","COL1A1","CA4","UNC5B","CD34"
      )


#Plot the UMAP Embedding with Marker Genes Overlayed w/ Imputation
p1<-plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix",pal = paletteContinuous(set = "blueYellow"), name = markerGenes, imputeWeights = getImputeWeights(proj))
plotPDF(plotList = p1, name = "GeneScoresSupplemental", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)
p1

#getting reproducible peak 
proj <- addGroupCoverages(ArchRProj = proj, force = TRUE,groupBy = "Clusters1")
#Call Reproducible Peaks w/ Macs2 (~5-10 minutes)
proj <- addReproduciblePeakSet(ArchRProj = proj,force = TRUE, groupBy = "Clusters1")
#Add Peak Matrix
proj <- addPeakMatrix(ArchRProj = proj,force = TRUE)



#adding the local accessibility around promoters for genes
addGeneScoreMatrix(proj,useTSS = TRUE,matrixName='promoter_Accessibility_1000',extendUpstream = c(1000,1000),
                   extendDownstream = c(1000, 1000))

#addMotifAnnotations
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force = TRUE)



#Identify Motif Enrichments
markersPeaks <- markerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "Clusters2")

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
#    
 )
plotPDF(heatmapEM, name = "Motifs-Enrich-Heatmap", width = 8, height = 12, ArchRProj = proj, addDOC = FALSE)


heatmapEM <- enrichHeatmap(enrichMotifs,pal = paletteContinuous(set = "solarExtra", n = 100))
heatmapEM


options(repr.plot.width = 18, repr.plot.height = 12)
heatmapEM

proj <- addBgdPeaks(proj,force = TRUE)
#Add chromVAR Deviations (~20-25 min if using CisBP Motif Set)
proj <- addDeviationsMatrix(ArchRProj = proj, peakAnnotation = "Motif",force = TRUE)
proj <- addImputeWeights(ArchRProj = proj)

getFeatures(proj, select ="NKX25|MESP1|MESP2|HEY1|HEY2|MEF2C|GATA4|TBX5|TCF21|TBX18|TBX1|WT1|WNT1|BMP1|PTCH1|SOX9|NFATC2|SRF|MEF2A|MEF2D|TEAD4|PRDM6|NFIC|ETV1|SOX11|SOX17|SRF|GATA5|TWIST2|HOXA2|ATF2|MEOX1|NR4A1|FOXP1|SNAI2", useMatrix = "MotifMatrix")


#plotting the tf deviations 
markerMotifs <- c('MEOX1_396','TBX5_781','SOX17_764','NKX25_550','TBX5_781','TCF21_39','SRF_641')
 
p1<-plotEmbedding(ArchRProj = proj, colorBy = "MotifMatrix", name = paste0("z:",markerMotifs), imputeWeights = getImputeWeights(proj))
p1


#adding the scRNA using CCA analysis
gene_integration_matrix<-readRDS('path to seurat rna rds object')
umap_rna<-Embeddings(object = gene_integration_matrix, reduction = "umap")
gene_integration_matrix



Idents(gene_integration_matrix)<-'new_manual_annotation'

gene_integration_matrix
new_all_rna<-gene_integration_matrix

gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$Sample==1,'time']="PCW6"
gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$Sample==2,'time']="PCW6"
gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$Sample==3,'time']="PCW19"
gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$Sample==4,'time']="PCW8"
gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$Sample==5,'time']="PCW12"

ventricular_rna<-c('Ventricular_myocytes')
atrial_rna<-c('Atrial_myocytes')
pericytes_rna<-c('Pericytes')
vSMC_rna<-c('vSMC')
Endocaridal_cushion_rna<-c('Endocardial_cushion_late')
Cardiac_fib_rna<-c('Cardiac_Fibroblast')
OFT_SMC_rna<-c('OFT_SMC_Progenitors')
Endocardium_rna<-c('Endocardium')
venal_rna<-c('Venal_endothelium')
capillary_rna<-c('Capillary_endothelium')
arterial_rna<-c('Arterial_endothelium')
lymph_ec_rna<-c('Lymph_ec')
Epicardium_rna<-c('Epicaridum')
Neuronal_rna<-c('Neuronal_1')
##
ventricular_atac<-c('Ventricular_Cardiomyocytes','Myocardium')
atrial_atac<-c('Atrial_Cardiomyocytes')
pericytes_atac<-c('Pericytes')
vSMC_atac<-c('vSMC')
Endocaridal_cushion_atac<-c('Endocardial_cushion','Endocardial_cushion_late')
Cardiac_fib_atac<-c('Early_Cardiac_fibroblast','Cardiac_fibroblast_progenitors','Cardiac_Fibroblast')
OFT_SMC_atac<-c('Vasculatur_development','Early_OFT_SMC')
Endocardium_atac<-c('Endocaridum','Endocaridum_transitioning')
venal_atac<-c('Venal_endothelium')
capillary_atac<-c('Capillary_endothelium')
arterial_atac<-c('Arterial_endothelium')
lymph_ec_atac<-c('Endocaridum_unidentifed')
Epicardium_atac<-c('Undifferentiated_epicardium')
Neuronal_atac<-c('Neural_Crest')


groupList <- SimpleList(
    Ventricular_myocytes = SimpleList(
        ATAC = proj$cellNames[proj$Clusters1 %in% ventricular_atac],
        RNA = rownames(gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$new_manual_annotation %in% ventricular_rna,])
    ),
    Atrial_myocytes = SimpleList(
        ATAC = proj$cellNames[proj$Clusters1 %in% atrial_atac],
        RNA = rownames(gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$new_manual_annotation %in% atrial_rna,])
    )  
    ,
    pericytes = SimpleList(
        ATAC = proj$cellNames[proj$Clusters1 %in% pericytes_atac],
        RNA = rownames(gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$new_manual_annotation %in% pericytes_rna,])
    )  
    ,
    vSMC = SimpleList(
        ATAC = proj$cellNames[proj$Clusters1 %in% vSMC_atac],
        RNA = rownames(gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$new_manual_annotation %in% vSMC_rna,])
    )  
    ,
    Endocardial_cushion = SimpleList(
        ATAC = proj$cellNames[proj$Clusters1 %in% Endocaridal_cushion_atac],
        RNA = rownames(gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$new_manual_annotation %in% Endocaridal_cushion_rna,])
    )  
    ,
    Cardiac_fib = SimpleList(
        ATAC = proj$cellNames[proj$Clusters1 %in% Cardiac_fib_atac],
        RNA = rownames(gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$new_manual_annotation %in% Cardiac_fib_rna,])
    )  
    ,
    OFT_smc = SimpleList(
        ATAC = proj$cellNames[proj$Clusters1 %in% OFT_SMC_atac],
        RNA = rownames(gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$new_manual_annotation %in% OFT_SMC_rna,])
    )  
    ,
    Endocardium_cells = SimpleList(
        ATAC = proj$cellNames[proj$Clusters1 %in% Endocardium_atac],
        RNA = rownames(gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$new_manual_annotation %in% Endocardium_rna,])
    ) 
    ,
    Arterial_cells = SimpleList(
        ATAC = proj$cellNames[proj$Clusters1 %in% arterial_atac],
        RNA = rownames(gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$new_manual_annotation %in% arterial_rna,])
    )  
    ,
    Capillary_cells = SimpleList(
        ATAC = proj$cellNames[proj$Clusters1 %in% capillary_atac],
        RNA = rownames(gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$new_manual_annotation %in% capillary_rna,])
    )  
    ,
    Venal_cells = SimpleList(
        ATAC = proj$cellNames[proj$Clusters1 %in% venal_atac],
        RNA = rownames(gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$new_manual_annotation %in% venal_rna,])
    )  
    ,
    lymph_cells = SimpleList(
        ATAC = proj$cellNames[proj$Clusters1 %in% lymph_ec_atac],
        RNA = rownames(gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$new_manual_annotation %in% lymph_ec_rna,])
    )  
    ,
    Epicardium = SimpleList(
        ATAC = proj$cellNames[proj$Clusters1 %in% c(Epicardium_atac,Neuronal_atac)],
        RNA = rownames(gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$new_manual_annotation %in% c(Epicardium_rna,Neuronal_rna),])
    ) 
    #,
    #Neuronal = SimpleList(
    #    ATAC = proj$cellNames[proj$Clusters1 %in% Neuronal_atac],
    #    RNA = rownames(gene_integration_matrix@meta.data[gene_integration_matrix@meta.data$manual_annotation %in% Neuronal_rna,])
    #) 
)

proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix_new",
    reducedDims = "IterativeLSI",
    seRNA = gene_integration_matrix,
    addToArrow = TRUE,
    groupATAC="Clusters1",
    groupList = groupList,
    groupRNA = "new_manual_annotation",
    nameCell = "predictedCell_new",
    nameGroup = "predictedGroup_new",
    nameScore = "predictedScore_new",
    sampleCellsRNA=30000,
    sampleCellsATAC = 30000,
    #threads = 1,
    force = TRUE
)

cm<-prop.table(table(proj$Clusters2, proj$predictedGroup_new))
a<-pheatmap(cm[c('C_Ventricular_Cardiomyocytes','A_Myocardium','B_Atrial_Cardiomyocytes','HA_Cardiac_Fibroblast',
   'H_Early_Cardiac_fibroblast','F_Cardiac_fibroblast_progenitors','E_valveFibroblast','G_valveFibroblast_late',
   'M_Undifferentiated_epicardium','D_Early_OFT_SMC','J_Vasculatur_development','K_vSMC','L_Pericytes',
   'N_Neural_Crest','O_Endocaridum','T_Venal_endothelium','Q_Endocaridum_unidentifed','R_Arterial_endothelium',
   'S_Capillary_endothelium'),c('Ventricular_myocytes','Atrial_myocytes','Cardiac_Fibroblast','Endocardial_cushion_late','Epicaridum',
 'OFT_SMC_Progenitors','vSMC','Pericytes','Neuronal_1','Endocardium','Venal_endothelium','Lymph_ec','Arterial_endothelium',
 'Capillary_endothelium')],scale='row',cluster_cols = FALSE,cluster_rows = FALSE,treeheight_row = 0, treeheight_col = 0)
plotPDF(a, name = "RNA-ATAC-Mapping", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)
