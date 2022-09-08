This folder contains the files required for the Congenital Heart Disorder mutation analysis

1. Folder "mutations" contains all subsetted mutations that overlap cluster ATAC peaks. The mutations are overlapped against the 1000bp peak regions that are used for train the BPNet models

2. Notebook "Scoring mutations-Cardiogenesis .ipynb" contains the script for scoring the overlapped CHD mutations using 5 fold BPNet models.

3. chd_final_0base.csv & chd_control_final_0base.csv contain the mutations from CHD cases and SSD controls as obtained from Richter et al. (2020)

4. The trained weight files for the BPNet models are available at https://doi.org/10.5281/zenodo.6789181

