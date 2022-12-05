Description of the scripts and notebooks used for motif filtering:

1. create_shap_inputs.ipynb : Create one-hot encoded sequences corresponding to cell type-specific peaks as inputs for DeepSHAP
2. create_shap_shuf_inputs.ipynb : Create one-hot encoded sequences for dinucleotide-shuffled peaks sequences
3. shap.counts.py : Run DeepSHAP using the counts head of the model on observed peak sequences
4. shap.counts.sh : Wrapper around shap.counts.py to run DeepSHAP for all cross-validation folds for each cell type
5. shuffled.shap.counts.py : Run DeepSHAP using the counts head of the model on dinucleotide-shuffled peak sequences
6. shuffled.shap.counts.sh :  Wrapper around shuffled.shap.counts.sh to run DeepSHAP for all cross-validation folds for each cell type
7. score_motifs.py : Get base-pair level counts predictions from the model at motif hits using the observed reference sequence
8. score_motifs.sh : Wrapper around score_motifs.py to get predictions for each cell type
9. score_motifs.shuf.py : Get base-pair level counts predictions from the model at motif hits with the motif positions replaced by dinucleotide-shuffled sequences
10. score_motifs.shuf.sh : Wrapper around score_motifs.shuf.py to get predictions with multiple random seeds for each cell type
11. get_mean_shap_scores.ipynb : Get mean DeepSHAP scores for each consensus peak using each cell type-specific model across all cross-validation folds
12. get_mean_predictions_celltype.ipynb : Get mean base-pair level counts predictions for each cell type-specific peak from the models trained on that cell type across all cross-validation folds
13. get_mean_predictions_consensus.ipynb : Get mean base-pair level counts predictions for each consensus peak from the models trained on each cell type across all cross-validation folds
14. get_shap_scores.ipynb : Use DeepSHAP scores for each motif in each cell type to call cell type-specific active motifs
15. get_motif_ism_scores.ipynb : Use motif-ISM scores for each motif in each cell type to call cell type-specific active motifs
16. merge_shap_ism_motifs.ipynb : Obtain the union of active motifs called using DeepSHAP and motif-ISM

