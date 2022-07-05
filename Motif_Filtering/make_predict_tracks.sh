#!/bin/bash

fold=0
bed_dir=/oak/stanford/groups/akundaje/projects/cardiogenesis/peaks_renamed
#pred_dir=/oak/stanford/groups/akundaje/laks/illuminafiles/cardiomyogenesis/BPNET_deepshaps/predict_consensus/mean
#bw_dir=/oak/stanford/groups/akundaje/laks/illuminafiles/cardiomyogenesis/BPNET_deepshaps/predict_consensus/bigwigs
pred_dir=/oak/stanford/groups/akundaje/laks/illuminafiles/cardiomyogenesis/BPNET_deepshaps/predict_celltype/mean
bw_dir=/oak/stanford/groups/akundaje/laks/illuminafiles/cardiomyogenesis/BPNET_deepshaps/predict_celltype/bigwigs
ref=/oak/stanford/groups/akundaje/refs/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
chrom_sizes=/oak/stanford/groups/akundaje/soumyak/refs/hg38/hg38.chrom.sizes

for cluster in all atrial cfmature cfwk8 epiicardium myocardium smcwk19 smcwk8 valvelate ventricular arteries capillary cfwk6 endocardium lymphec neuralcrest pericytes smcwk6 valveearly veins
do
    echo $cluster
    python make_predict_bw.py -bed $bed_dir/peaks_${cluster}_bpnet.csv \
                              -npy $pred_dir/$cluster.mean.preds.npy \
                              -fasta $ref \
                              -c $chrom_sizes \
                              --type count_shap \
                              --tqdm 1 \
                              -o $bw_dir/$cluster.mean.preds.bw \
                              -s $bw_dir/$cluster.mean.preds.stats.txt
done

