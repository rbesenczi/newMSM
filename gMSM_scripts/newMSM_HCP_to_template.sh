#!/bin/bash -l

###########################################################
## The following few lines are for running on CREATE cluster. Ignore it otherwise.
#SBATCH --partition=cpu
#SBATCH --time=0-12:00
#SBATCH --nodes=1
#SBATCH --mem=8192
#SBATCH --ntasks=8
#SBATCH --job-name=newMSM
#SBATCH --output=./newMSM_HCP_to_tmp_subs_NODE2214_39.txt
## how to run on create: $ sbatch newMSM_HCP_to_template.sh
###########################################################

###########################################################
## Set the following lines according to your settings and cohort.
workdir=$HOME/newMSM_HCP_to_tmp
input_folder=/scratch/prj/cortical_imaging/Renato/affined_sulc
template=/scratch/prj/cortical_imaging/Renato/groupwise/HCP/results/NODE2214_HPC/groupwise.NODE2214_HPC.mean.sulc.affine.dedrifted.ico6.shape.gii
mesh=$workdir/template/sunet.ico-6.template.surf.gii
config_file=$workdir/configs/newMSM_HCP_to_template_config_NODE2214.txt
sublist=/scratch/prj/cortical_imaging/Renato/groupwise/HCP/group_lists/NODE2214_HPC.csv
outdir=/scratch/prj/cortical_imaging/Renato/HCP_to_template/2
bunch=39
###########################################################

from=$(( bunch * 20 + 1 ))
to=$(( from + 19 ))

while IFS=',' read -r subject group
do
  echo "Registering subject $subject to the template..."

  time $HOME/fsldev/bin/newmsm \
  --inmesh=$mesh \
  --refmesh=$mesh \
  --indata="$input_folder/$subject.sulc.affine.ico6.shape.gii" \
  --refdata=$template \
  --conf=$config_file \
  --out=$outdir/$subject.MSMSulc.ico6. \
  --verbose

done < <(sed -n "${from},${to}p" $sublist)
