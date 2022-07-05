#!/bin/bash

cluster=$1
gpu=$2
model_dir=/home/soumyak/cardiogenesis/models
onehot_dir=/home/soumyak/cardiogenesis/onehot_peaks
shap_dir=/home/soumyak/cardiogenesis/deepshap

CUDA_VISIBLE_DEVICES=$gpu

[[ -d $shap_dir/$cluster ]] || mkdir $shap_dir/$cluster

for fold in 0 1 2 3 4
do
    echo $cluster
    echo Fold-$fold

    python shap.counts.py --model_path $model_dir/$cluster.$fold \
                          --onehot $onehot_dir/${cluster}_onehot.npy \
                          --out_prefix $shap_dir/$cluster/fold$fold
done

