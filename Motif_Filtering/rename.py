import os
import sys
import shutil

name_dict = {'all': 'pseudobulk',
             'arteries': 'aec',
             'atrial': 'acm',
             'capillary': 'cap',
             'cfmature': 'CF',
             'cfwk6': 'CFP',
             'cfwk8': 'preCF',
             'endocardium': 'Endo1_2',
             'epiicardium': 'epc',
             'lymphec': 'lec',
             'myocardium': 'ecm',
             'neuralcrest': 'nc',
             'pericytes': 'pc',
             'smcwk19': 'smc',
             'smcwk6': 'oft',
             'smcwk8': 'preSMC',
             'valveearly': 'fb1',
             'valvelate': 'fb2',
             'veins': 'vec',
             'ventricular': 'vcm'
            }

intersect_indir = '/oak/stanford/groups/akundaje/laks/illuminafiles/cardiomyogenesis/BPNET_deepshaps/motifs_in_peaks'
intersect_outdir = '/oak/stanford/groups/akundaje/projects/cardiogenesis/for_anshul/motifs_in_peaks'

ism_indir = '/oak/stanford/groups/akundaje/laks/illuminafiles/cardiomyogenesis/BPNET_deepshaps/ism_cleaned_motifs'
ism_outdir = '/oak/stanford/groups/akundaje/projects/cardiogenesis/for_anshul/ism_cleaned_motifs'

shap_indir = '/oak/stanford/groups/akundaje/laks/illuminafiles/cardiomyogenesis/BPNET_deepshaps/shap_cleaned_motifs'
shap_outdir = '/oak/stanford/groups/akundaje/projects/cardiogenesis/for_anshul/shap_cleaned_motifs'

all_indir = '/oak/stanford/groups/akundaje/laks/illuminafiles/cardiomyogenesis/BPNET_deepshaps/all_cleaned_motifs'
all_outdir = '/oak/stanford/groups/akundaje/projects/cardiogenesis/for_anshul/all_cleaned_motifs'

for celltype in name_dict:
    print(celltype)
    shutil.copyfile(intersect_indir + '/' + celltype + '.motifs.in_peaks.bed',
                    intersect_outdir + '/' + name_dict[celltype] + '.motifs.in_peaks.bed')

    shutil.copyfile(ism_indir + '/' + celltype + '.active_motifs.bed',
                    ism_outdir + '/' + name_dict[celltype] + '.active_motifs.bed')
    shutil.copyfile(ism_indir + '/' + celltype + '.active_motifs.tsv',
                    ism_outdir + '/' + name_dict[celltype] + '.active_motifs.tsv')

    shutil.copyfile(shap_indir + '/' + celltype + '.active_motifs.bed',
                    shap_outdir + '/' + name_dict[celltype] + '.active_motifs.bed')
    shutil.copyfile(shap_indir + '/' + celltype + '.active_motifs.tsv',
                    shap_outdir + '/' + name_dict[celltype] + '.active_motifs.tsv')

    shutil.copyfile(all_indir + '/' + celltype + '.active_motifs.tsv',
                    all_outdir + '/' + name_dict[celltype] + '.active_motifs.tsv')

    #print('tail -n +2 ' + all_outdir + '/' + name_dict[celltype] + '.active_motifs.tsv '                                                               +  all_outdir + '/' + name_dict[celltype] + '.active_motifs.bed')

    os.system('tail -n +2 ' + all_outdir + '/' + name_dict[celltype] + '.active_motifs.tsv > '
              +  all_outdir + '/' + name_dict[celltype] + '.active_motifs.bed')
