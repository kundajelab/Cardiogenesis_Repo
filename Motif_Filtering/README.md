Description of the scripts and notebooks used for motif filtering:

1. create_shap_inputs.ipynb : Create one-hot encoded sequences corresponding to cell type-specific peaks as inputs for DeepSHAP
2. create_shap_shuf_inputs.ipynb : Create one-hot encoded sequences for dinucleotide-shuffled peaks sequences
3. shap.counts.py : Run DeepSHAP using the counts head of the model on observed peak sequences
4. shap.counts.sh : Wrapper around shap.counts.py to run DeepSHAP for all cross-validation folds for each cell type
5. shuffled.shap.counts.py : Run DeepSHAP using the counts head of the model on dinucleotide-shuffled peak sequences
6. shuffled.shap.counts.sh :  Wrapper around shuffled.shap.counts.sh to run DeepSHAP for all cross-validation folds for each cell type
