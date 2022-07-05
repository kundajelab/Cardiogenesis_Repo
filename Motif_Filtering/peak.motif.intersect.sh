#!/bin/bash

for cluster in all atrial cfmature cfwk8 epiicardium myocardium smcwk19 smcwk8 valvelate ventricular arteries capillary cfwk6 endocardium lymphec neuralcrest pericytes smcwk6 valveearly veins
do
    echo $cluster
    bedtools intersect -wa -wb -f 1.0 -a /mnt/lab_data3/soumyak/cardiogenesis/motif_ism/motif.hits.noheader.bed -b /mnt/lab_data3/soumyak/cardiogenesis/peaks/peaks_${cluster}_bpnet.csv > /mnt/lab_data3/soumyak/cardiogenesis/motifs_in_peaks/$cluster.motifs.in_peaks.bed
done

