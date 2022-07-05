#!/bin/bash

cluster=$1
gpu=$2
ism_dir=/home/soumyak/cardiogenesis/motif_ism
table=$ism_dir/motif.hits.bed
fasta=/home/soumyak/refs/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
model_dir=/home/soumyak/cardiogenesis/models

echo $gpu
CUDA_VISIBLE_DEVICES=$gpu
echo $CUDA_VISIBLE_DEVICES

[[ -d $ism_dir/$cluster ]] || mkdir $ism_dir/$cluster

for fold in 0
do

    for shuf in 1 2
    do

        echo $cluster
        echo Fold-$fold
        echo Shuf-$shuf

        [[ -d $ism_dir/$cluster/fold$fold ]] || mkdir $ism_dir/$cluster/fold$fold

        python score_motifs.shuf.py --model_path $model_dir/$cluster.$fold \
                                    --table $table \
                                    --fasta $fasta \
                                    --index $shuf \
                                    --out_prefix $ism_dir/$cluster/fold$fold
    done
done

