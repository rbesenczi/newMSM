#!/bin/bash -l

###########################################################
## Set the following lines according to your settings and cohort.
dataset="HCP"
input_folder=/scratch/users/k2258483/affined_${dataset}_sulc_curv
workdir=/scratch/users/k2258483/${dataset}_to_template
template=$workdir/templates/MSMStrain.L.sulc.curv.ico6.shape.gii
mesh=$workdir/templates/sunet.ico-6.sphere.surf.gii
config_file=$workdir/configs/${dataset}_to_template_config.txt
sublist=$workdir/all_subjects.txt
outdir=$workdir/output/
logdir=$workdir/logs/
###########################################################

num_subjects=$(cat $sublist | wc -l)

sbatch --array=1-${num_subjects} --output=$logdir/registration_%a.txt --export=ALL newMSM_HCP_to_template_v2.sh
