This folder ArchR notebooks for all single cell analysis presented manuscript and the in vivo to in vitro nearest neighbour comparisons.

1. Cardiogenesis - Plotting footprints with and without cleaning.ipynb - Notebook contains the modified footprinting function adapted from ArchR for plotting the naive overlap motif and BPNet cleaned motif footprints.
2. Cardiogenesis - Fetal + Invitro combined project.ipynb - Notebook contains the ArchR code for creating the in vivo and in vitro combined ArchR project, followed by peak calling and motif calling. 
	a) Fragment files are available on GEO to create arrow files
	b) fetal_hearts_names_Final.csv - Post filtering in vivo barcodes required for subsetting the ArchR project
	c) differentiation_names_Final.csv - Post filtering in vitro barcodes required for subsetting the ArchR project
3. Cardiogenesis - Invitro cells clustering, peak calling, motif annotation.ipynb - Notebook contains the ArchR code for creating the in vitro (iPSC derived cardiac cells) ArchR project, followed by peak calling, motif calling and integrating scRNA data from Friedman et al.
	a) differentiation_names.csv - Post filtering in vitro barcodes required for subsetting the ArchR project
	b) The cleaned scRNA object is available in zenodo - friedman_final_cleaned.rds https://zenodo.org/record/7392252
4. Cardiogenesis Clustering , Peak Calling, Motif Annotation & ATAC RNA integration.ipynb - Notebook contains the ArchR code for creating the in vivo (fetal cardiac cells) ArchR project, followed by peak calling, motif calling and integrating scRNA data from roughly matched public datasources.
	a) Final_barcodes_fetal_heart_NOMICROPHAGES.csv - Post filtering in vivo barcodes required for subsetting the ArchR project
	b) The cleaned scRNA object is available in zenodo - human_6_8_12and19_merged_final.rds https://zenodo.org/record/7392252
5. Cardiogenesis - Projecting invitro cells on invivo cells & DIfferent enhancers and TFs.ipynb - Notebook contains the code to perform the differential analysis between fetal heart in vivo cells and iPSC derived cardiac cell types (in vitro cells), nearest neighbour differential analysis and motif enrichment.
6. Cardiogenesis_Projecting invitro on to invivo & Trajectories & Nearest neighbours- NEW EPC.ipynb - Notebook contains the code to perform the differential analysis between fetal heart in vivo cells and new iPSC derived Epicardial cells (EPC), nearest neighbour analysis and differential enhancer.
