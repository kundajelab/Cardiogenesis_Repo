#!/bin/bash
fold=4
outdir=.
model=c7
seed=1234
gpu=2
title='bias corrected '
./c7_train.sh $fold $gpu $model $seed $outdir
./c7_predict.sh $fold $gpu $model $seed $outdir 
./c7_score.sh $outdir $model $fold $title 
