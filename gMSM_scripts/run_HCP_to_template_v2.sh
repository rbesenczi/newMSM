#!/bin/bash -l

###########################################################
## Set the following lines according to your settings and cohort.
dataset="HCP"
input_folder=/scratch/prj/cortical_imaging/Renato/affined_${dataset}_sulc
workdir=/scratch/prj/cortical_imaging/Renato/${dataset}_to_template
template=$workdir/MSMSulc.L.sulc.ico6.shape.gii
mesh=$workdir/sunet.ico-6.sphere.surf.gii
config_file=$workdir/${dataset}_to_template_config.txt
sublist=$workdir/subjects_in_study.txt
outdir=$workdir/output/
logdir=$workdir/logs/
###########################################################

num_subjects=$(cat $sublist | wc -l)

mkdir $outdir
mkdir $logdir

sbatch --array=1-${num_subjects} --output=$logdir/registration_%{a}.txt --export=ALL newMSM_HCP_to_template_v2.sh
