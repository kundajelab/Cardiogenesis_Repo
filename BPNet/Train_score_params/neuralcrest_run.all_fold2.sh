#!/bin/bash
fold=2
outdir=.
model=c7
seed=1234
gpu=1
title='bias corrected '
./c7_train.sh $fold $gpu $model $seed $outdir
./c7_predict.sh $fold $gpu $model $seed $outdir 
./c7_score.sh $outdir $model $fold $title 
