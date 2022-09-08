This folder scripts and objects required to obtain optimal transport inferred scATAC trajectories.

1.Genescore_foropt_fetalheart1.loom
The Genescore_foropt_fetalheart1.loom contains the seurat object of cell by gene score matrix. The cell by gene score matrix is extracted from ArchR using getMatrixFromProject() and converted to seurat object and stored as the loom object

2.Cell_details.csv - This file contains the metadata for all cellular barcodes in the fetal atlas

3.adata_obs_for_peaks.csv - This file contains the metadata for cellular barcodes including initial cell growth rate estimates computed using the gene score matrix ( Cardiogenesis - Optimal Transport constructing the growth rates .ipynb)

4.coords_df.csv - This file contains the 2 dimensional UMAP coordinates produced through iterativeLSI using ArchR

5.major_cell_set.gmt - Dictionary contain clusters names as keys and cell barcodes as values

6.GeneScore_OT_g.txt- Cell growth rate estimates for cell barcode across the 3 iterations of OT model

7.Creating the initial growth rates
	a) Cardiogenesis - Optimal Transport constructing the growth rates .ipynb notebook contains scripts to calculate the initial growth rate estimates
	b) cellcycle_apoptosis.gmx - contains the set of cell cycle and apoptosis gene sets from the original manusciprt. (Geoffrey Schiebinger, Jian Shu, Marcin Tabaka, Brian Cleary, et. al.)
	c) growth_gs_init.txt is the file with the initial growth estimates

8.Constructing optimal transport infered trajectories
	a) Optimal Transport GeneScore -Cardiogenesis trajectories - .ipynb notebook contains scripts to calculate the trajectories
	b) Cell metadata for the adata object are initiated with the precomputed cell growth rates as described in point 7.
	c) Computed OT model weight (h5ad) objects are stored in tmaps/GeneScore_OT
	
