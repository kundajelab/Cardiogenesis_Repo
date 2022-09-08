Cardiogenesis Repo This repo contains scripts for the single cell analysis and BPNet analysis for the Cardiogenesis dataset

1. The single cell analysis code is based off ArchR (https://www.archrproject.com/index.html) and the notebooks contain modifications to standard archr functions when they are used in the manuscript. 
2. The BPNet code base uses two other repos : KerasAC (https://zenodo.org/record/4248179#.X8skj5NKiF0) and seqdataloader (https://zenodo.org/record/3771365#.X8skqZNKiF0) as part of the data processing and model training scripts.
3. BPNet models are available at Zenodo at https://doi.org/10.5281/zenodo.6789181
4. The folders are arranged for different analysis presented in the manuscript
	a) BPNet - contains all scripts and files for training the BPNet models and scoring the congenital heart disorder mutations
	b) Motif_Filtering - contains all scripts and files required for cleaning the naive motif calls using BPNet models
	c) OptimalTransport - contains all scripts and files required for trajectory analysis from the fetal atlas
	d) R_notebooks - contains scripts and code for all ArchR based single cell analysis and the in vivo to in vitro nearest neighbour analysis
