#!/bin/bash -l

#SBATCH --partition=cpu
#SBATCH --time=0-1:00
#SBATCH --nodes=1
#SBATCH --mem=2048
#SBATCH --ntasks=2
#SBATCH --job-name=MSM

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

while IFS=',' read -r subject
do
  echo "Registering subject $subject to the template..."

  time $HOME/fsldev/bin/newmsm \
    --inmesh=$mesh \
    --refmesh=$mesh \
    --indata="$input_folder/$subject.sulc.curv.affine.ico6.shape.gii" \
    --refdata=$template \
    --conf=$config_file \
    --out=$outdir/$subject.MSMSulc.ico6. \
    --verbose

  wb_command -surface-distortion $mesh $outdir/$subject.MSMSulc.ico6.sphere.reg.surf.gii $outdir/$subject.MSMSulc.ico6.sphere.distortion.func.gii -local-affine-method -log2

done < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $sublist)
