#!/bin/bash

fold=0
bed_dir=/oak/stanford/groups/akundaje/projects/cardiogenesis/peaks_renamed
shap_dir=/oak/stanford/groups/akundaje/laks/illuminafiles/cardiomyogenesis/BPNET_deepshaps/deepshap_consensus/mean
bw_dir=/oak/stanford/groups/akundaje/laks/illuminafiles/cardiomyogenesis/BPNET_deepshaps/deepshap_consensus/dynseq
ref=/oak/stanford/groups/akundaje/refs/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
chrom_sizes=/oak/stanford/groups/akundaje/soumyak/refs/hg38/hg38.chrom.sizes

#for cluster in all atrial cfmature cfwk8 epiicardium myocardium smcwk19 smcwk8 valvelate ventricular arteries capillary cfwk6 endocardium lymphec neuralcrest pericytes smcwk6 valveearly veins
for cluster in veins
do
    echo $cluster
    python make_importance_bw.py -bed $bed_dir/peaks_all_bpnet.csv \
                                 -npy $shap_dir/$cluster.mean.count.shap.npy \
                                 -fasta $ref \
                                 -c $chrom_sizes \
                                 --type count_shap \
                                 --tqdm 1 \
                                 -o $bw_dir/$cluster.mean.count.shap.bw \
                                 -s $bw_dir/$cluster.mean.count.shap.stats.txt
done

